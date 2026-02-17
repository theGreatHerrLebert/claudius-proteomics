#!/usr/bin/env python3
"""
Modification Profile Resolver for San Jose.

Auto-detects the appropriate modification profile (standard, phospho, hla)
from paper/PRIDE metadata. Allows manual override via config.

Resolution order:
1. Manual override in config["dataset_metadata"][accession]["mod_profile"]
2. Auto-detect from paper_extraction.yaml (study_type, keywords)
3. Default: "standard"
"""

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


def resolve_mod_profile(
    accession: str,
    config: Dict[str, Any],
    metadata_dir: Path,
) -> Tuple[str, str, Optional[Dict[str, Any]]]:
    """Resolve modification profile for a dataset.

    Args:
        accession: PRIDE accession
        config: Pipeline configuration dict
        metadata_dir: Path to metadata directory for this dataset

    Returns:
        (profile_name, source, paper_search_settings)
        - profile_name: "standard", "phospho", or "hla"
        - source: how the profile was determined ("manual", "auto:study_type",
          "auto:keywords", "default")
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

    # 3. Default
    return "standard", "default", paper_settings


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
