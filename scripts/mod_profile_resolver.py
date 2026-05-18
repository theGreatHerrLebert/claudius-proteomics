#!/usr/bin/env python3
"""
Modification Profile Resolver for San Jose.

Auto-detects the appropriate modification profile (standard, phospho, hla)
from paper/PRIDE metadata. Allows manual override via config.

Resolution order:
1. Manual override in config["dataset_metadata"][accession]["mod_profile"]
2. Auto-detect from paper_extraction.yaml (study_type, keywords)
3. Auto-detect from .d file names (classI/classII/HLA/MHC tokens)
4. Auto-detect from pride_metadata.yaml (title, sample prep)
5. Default: "standard"

Steps 3-4 are PDF-independent — they catch HLA immunopeptidomics datasets
even when no paper was found, which is the common case in bulk reprocessing.
"""

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml


# Keywords that trigger profile auto-detection
_HLA_KEYWORDS = {"hla", "mhc", "immunopeptid", "immunopeptidomics", "mhc-i", "mhc-ii"}
_PHOSPHO_KEYWORDS = {"phospho", "phosphoproteomics", "phosphorylation", "enrichment"}

# Study type -> profile mapping
_STUDY_TYPE_PROFILES = {
    "phosphoproteomics": "phospho",
    "immunopeptidomics": "hla",
}

# Immunopeptidomics tokens in .d run/file names. Matched against whole tokens
# (filename split on non-alphanumerics) so "classic" etc. cannot false-trigger.
# The class indicator (roman numeral or digit) is optional.
_HLA_TOKEN_RE = re.compile(
    r"^(?:hla|mhc)(?:[12]|i{1,2})?$"   # hla, mhc, hla2, mhci, hlaii, mhcii, ...
    r"|^class(?:i{1,2}|[12])$"         # classi, classii, class1, class2
    r"|^immunopep"                     # immunopeptidome(s) / immunopeptidomics
)


def resolve_mod_profile(
    accession: str,
    config: Dict[str, Any],
    metadata_dir: Path,
    run_names: Optional[List[str]] = None,
) -> Tuple[str, str, Optional[Dict[str, Any]]]:
    """Resolve modification profile for a dataset.

    Args:
        accession: PRIDE accession
        config: Pipeline configuration dict
        metadata_dir: Path to metadata directory for this dataset
        run_names: .d run/file names for this dataset, used as a PDF-independent
          signal for HLA immunopeptidomics detection

    Returns:
        (profile_name, source, paper_search_settings)
        - profile_name: "standard", "phospho", or "hla"
        - source: how the profile was determined ("manual", "auto:study_type",
          "auto:keyword", "auto:filename", "auto:pride_metadata", "default")
        - paper_search_settings: extracted search_settings dict from paper, or None
    """
    paper_settings = _load_paper_search_settings(metadata_dir)

    # 1. Manual override from config
    dataset_meta = config.get("dataset_metadata", {}).get(accession, {})
    manual_profile = dataset_meta.get("mod_profile")
    if manual_profile:
        available = config.get("mod_profiles", {})
        if manual_profile in available:
            return manual_profile, "manual", paper_settings
        else:
            print(f"  Warning: mod_profile '{manual_profile}' not in config, ignoring")

    # 2. Auto-detect from paper_extraction.yaml
    paper_data = _load_paper_extraction(metadata_dir)
    if paper_data:
        profile, source = _detect_from_paper(paper_data)
        if profile:
            return profile, source, paper_settings

    # 3. Auto-detect from .d file names — strong, PDF-independent HLA signal
    profile, source = _detect_from_filenames(run_names)
    if profile:
        return profile, source, paper_settings

    # 4. Auto-detect from PRIDE submission metadata (title + sample prep)
    profile, source = _detect_from_pride_metadata(metadata_dir)
    if profile:
        return profile, source, paper_settings

    # 5. Default
    return "standard", "default", paper_settings


def _detect_from_filenames(
    run_names: Optional[List[str]],
) -> Tuple[Optional[str], str]:
    """Detect an HLA immunopeptidomics dataset from .d run/file names.

    Immunopeptidomics runs almost always carry classI/classII/HLA/MHC tokens
    in their file names. This works even when no paper PDF is available.
    """
    for name in run_names or []:
        for token in re.split(r"[^a-z0-9]+", str(name).lower()):
            if token and _HLA_TOKEN_RE.match(token):
                return "hla", f"auto:filename={token}"
    return None, ""


def _detect_from_pride_metadata(metadata_dir: Path) -> Tuple[Optional[str], str]:
    """Detect profile from PRIDE submission metadata (title + sample prep).

    pride_metadata.yaml is populated from the PRIDE API and is available even
    when no paper PDF was found, so it is a useful independent fallback.
    """
    path = Path(metadata_dir) / "pride_metadata.yaml"
    if not path.exists():
        return None, ""
    try:
        with open(path) as f:
            data = yaml.safe_load(f) or {}
    except Exception:
        return None, ""

    texts: List[str] = []
    title = data.get("title", {})
    if isinstance(title, dict) and isinstance(title.get("value"), str):
        texts.append(title["value"].lower())
    sample_prep = data.get("fields", {}).get("sample_prep", {})
    if isinstance(sample_prep, dict) and isinstance(sample_prep.get("value"), str):
        texts.append(sample_prep["value"].lower())
    combined = " ".join(texts)

    for kw in _HLA_KEYWORDS:
        if kw in combined:
            return "hla", f"auto:pride_metadata={kw}"
    if "phospho" in combined:
        return "phospho", "auto:pride_metadata=phospho"
    return None, ""


def _load_paper_extraction(metadata_dir: Path) -> Optional[Dict[str, Any]]:
    """Load paper_extraction.yaml if it exists."""
    path = Path(metadata_dir) / "paper_extraction.yaml"
    if not path.exists():
        return None
    try:
        with open(path) as f:
            return yaml.safe_load(f)
    except Exception:
        return None


def _load_paper_search_settings(metadata_dir: Path) -> Optional[Dict[str, Any]]:
    """Load search_settings from paper_extraction.yaml."""
    data = _load_paper_extraction(metadata_dir)
    if not data:
        return None
    fields = data.get("fields", {})
    search_settings = fields.get("search_settings", {})
    value = search_settings.get("value") if isinstance(search_settings, dict) else None
    return value if isinstance(value, dict) else None


def _detect_from_paper(paper_data: Dict[str, Any]) -> Tuple[Optional[str], str]:
    """Auto-detect profile from paper extraction data.

    Returns:
        (profile_name or None, source_description)
    """
    fields = paper_data.get("fields", {})

    # Check study_type
    study_type_field = fields.get("study_type", {})
    study_type = study_type_field.get("value", "") if isinstance(study_type_field, dict) else ""
    if isinstance(study_type, str):
        study_type_lower = study_type.lower().strip()
        if study_type_lower in _STUDY_TYPE_PROFILES:
            return _STUDY_TYPE_PROFILES[study_type_lower], f"auto:study_type={study_type_lower}"

    # Check for HLA keywords across multiple fields
    text_fields = []
    for field_name in ("study_type", "sample_prep", "cell_lines_tissues"):
        field_data = fields.get(field_name, {})
        value = field_data.get("value", "") if isinstance(field_data, dict) else ""
        if isinstance(value, str):
            text_fields.append(value.lower())
        elif isinstance(value, list):
            text_fields.extend(str(v).lower() for v in value)

    combined_text = " ".join(text_fields)
    for kw in _HLA_KEYWORDS:
        if kw in combined_text:
            return "hla", f"auto:keyword={kw}"

    # Check for phospho keywords in variable_modifications
    search_settings = fields.get("search_settings", {})
    ss_value = search_settings.get("value") if isinstance(search_settings, dict) else None
    if isinstance(ss_value, dict):
        var_mods = ss_value.get("variable_modifications", [])
        if isinstance(var_mods, list):
            for mod in var_mods:
                mod_lower = str(mod).lower()
                if "phospho" in mod_lower:
                    return "phospho", "auto:variable_modifications=phospho"

    # Check for phospho keywords in text fields
    for kw in _PHOSPHO_KEYWORDS:
        if kw in combined_text:
            return "phospho", f"auto:keyword={kw}"

    return None, ""


def compare_profile_vs_paper(
    profile_name: str,
    profile_config: Dict[str, Any],
    paper_search_settings: Optional[Dict[str, Any]],
) -> List[str]:
    """Compare resolved profile against paper-extracted search settings.

    Returns human-readable diff strings for logging.

    Args:
        profile_name: Name of the resolved profile
        profile_config: The profile dict from config["mod_profiles"][name]
        paper_search_settings: search_settings extracted from paper, or None
    """
    if not paper_search_settings:
        return []

    diffs = []

    # Compare fixed modifications
    paper_fixed = paper_search_settings.get("fixed_modifications", [])
    if isinstance(paper_fixed, list) and paper_fixed:
        profile_fixed_names = {
            fm["name"] for fm in profile_config.get("fixed_modifications", [])
        }
        paper_fixed_names = set()
        for mod in paper_fixed:
            # Paper mods are like "Carbamidomethyl (C)"
            name = str(mod).split("(")[0].strip()
            paper_fixed_names.add(name)
        if profile_fixed_names != paper_fixed_names:
            diffs.append(
                f"Fixed mods: profile={sorted(profile_fixed_names)}, "
                f"paper={sorted(paper_fixed_names)}"
            )

    # Compare variable modifications
    paper_var = paper_search_settings.get("variable_modifications", [])
    if isinstance(paper_var, list) and paper_var:
        profile_var_names = {
            vm["name"] for vm in profile_config.get("variable_modifications", [])
        }
        paper_var_names = set()
        for mod in paper_var:
            name = str(mod).split("(")[0].strip()
            paper_var_names.add(name)
        if profile_var_names != paper_var_names:
            diffs.append(
                f"Variable mods: profile={sorted(profile_var_names)}, "
                f"paper={sorted(paper_var_names)}"
            )

    # Compare max_variable_mods
    paper_max_var = paper_search_settings.get("max_variable_mods")
    profile_max_var = profile_config.get("max_variable_mods")
    if paper_max_var is not None and profile_max_var is not None:
        if int(paper_max_var) != int(profile_max_var):
            diffs.append(
                f"Max variable mods: profile={profile_max_var}, paper={paper_max_var}"
            )

    # Compare peptide length range
    paper_min_len = paper_search_settings.get("min_peptide_length")
    paper_max_len = paper_search_settings.get("max_peptide_length")
    profile_min_len = profile_config.get("min_peptide_length")
    profile_max_len = profile_config.get("max_peptide_length")
    if paper_min_len is not None and profile_min_len is not None:
        if int(paper_min_len) != int(profile_min_len):
            diffs.append(
                f"Min peptide length: profile={profile_min_len}, paper={paper_min_len}"
            )
    if paper_max_len is not None and profile_max_len is not None:
        if int(paper_max_len) != int(profile_max_len):
            diffs.append(
                f"Max peptide length: profile={profile_max_len}, paper={paper_max_len}"
            )

    return diffs
