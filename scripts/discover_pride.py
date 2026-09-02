#!/usr/bin/env python3
"""
Automatic timsTOF Dataset Discovery from PRIDE.

Queries the PRIDE REST API, classifies each dataset by instrument, acquisition
mode, organism, and lab, then scores for diversity and outputs a ranked catalog.

Three stages:
  1. Paginated search → candidate projects (includes metadata)
  2. Instrument validation → keep only real timsTOF projects
  3. File manifest check → confirm Bruker .d files exist

Usage:
    python scripts/discover_pride.py                          # Full discovery
    python scripts/discover_pride.py --cached-only            # Re-score from cache
    python scripts/discover_pride.py --top 20                 # Show top 20
    python scripts/discover_pride.py --mode DDA --top 50      # Filter to DDA
    python scripts/discover_pride.py --skip-file-check        # Faster (no Stage 3)
    python scripts/discover_pride.py --list-top 5             # Accessions only
    python scripts/discover_pride.py --update-config config/config.yaml
    python scripts/discover_pride.py --incremental            # Only new since last run
"""

import argparse
import csv
import fnmatch
import json
import math
import re
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import requests
import yaml

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PRIDE_API_BASE = "https://www.ebi.ac.uk/pride/ws/archive/v3"
# v2 is an alias for v3 on /projects and /search/projects, but its
# files/byProject endpoint was retired and now answers 200 with an empty
# body — see FileCheckError below.
FILE_PAGE_SIZE = 1000

DEFAULT_CACHE_DIR = Path("data/discovery/cache")
DEFAULT_OUTPUT_DIR = Path("data/discovery")

# Known timsTOF instrument model variants
TIMSTOF_MODELS = {
    "timstof": "timsTOF",
    "timstof pro": "timsTOF Pro",
    "timstof pro 2": "timsTOF Pro 2",
    "timstof flex": "timsTOF fleX",
    "timstof ht": "timsTOF HT",
    "timstof scp": "timsTOF SCP",
    "timstof ultra": "timsTOF Ultra",
}

# Bruker raw data file patterns (from download_pride.py)
BRUKER_FILE_PATTERNS = ["*.d.zip", "*.d.tar", "*.d"]

# Software that hints at acquisition mode
DIA_SOFTWARE = {"dia-nn", "diann", "spectronaut", "openswath", "encyclopedia"}
DDA_SOFTWARE = {"maxquant", "fragpipe", "msfragger", "proteome discoverer", "mascot"}

# Organism aliases (from sample_group_resolver.py)
ORGANISM_ALIASES: Dict[str, str] = {
    "homo sapiens": "human",
    "mus musculus": "mouse",
    "rattus norvegicus": "rat",
    "saccharomyces cerevisiae": "yeast",
    "escherichia coli": "ecoli",
    "drosophila melanogaster": "drosophila",
    "caenorhabditis elegans": "c_elegans",
    "arabidopsis thaliana": "arabidopsis",
    "danio rerio": "zebrafish",
    "sus scrofa": "pig",
    "bos taurus": "bovine",
    "gallus gallus": "chicken",
}

# Diversity scoring weights
DEFAULT_WEIGHTS = {
    "instrument_model": 0.25,
    "organism": 0.20,
    "acquisition_mode": 0.20,
    "lab": 0.15,
    "tissue": 0.10,
    "dataset_size": 0.10,
}


# ---------------------------------------------------------------------------
# Acquisition mode classification (extends pride_metadata._parse_protocol_text)
# ---------------------------------------------------------------------------

def _parse_protocol_text_for_mode(text: str) -> Tuple[Optional[str], float, str]:
    """Parse protocol/description text for acquisition mode using regex.

    Returns (mode, confidence, source_field).
    """
    if not text:
        return None, 0.0, ""

    patterns = [
        (r"MALDI[- ]?imaging|MALDI[- ]?MSI|MALDI[- ]?IMS", "MALDI", 0.95),
        (r"single[- ]?cell|SCoPE[- ]?MS|nanoPOTS|cellenONE", "single-cell", 0.90),
        (r"diaPASEF|dia-?PASEF", "diaPASEF", 0.95),
        (r"(?<!dia)PASEF", "PASEF", 0.90),
        (r"prm-?PASEF|PRM[- ]?PASEF", "PRM", 0.90),
        (r"data[- ]?independent\s+acquisition|(?<!\w)DIA(?!\w)", "DIA", 0.85),
        (r"data[- ]?dependent\s+acquisition|(?<!\w)DDA(?!\w)", "DDA", 0.85),
        (r"parallel\s+reaction\s+monitoring|(?<!\w)PRM(?!\w)", "PRM", 0.85),
    ]
    for pattern, mode, confidence in patterns:
        if re.search(pattern, text, re.IGNORECASE):
            return mode, confidence, "protocol_text"

    return None, 0.0, ""


def classify_acquisition_mode(
    project: Dict[str, Any],
) -> Tuple[str, float, str]:
    """Classify acquisition mode using layered heuristics.

    Priority:
      1. experimentTypes (structured field)
      2. Protocol/description text (regex)
      3. Keywords
      4. Software inference
      5. Title (last resort)

    Returns (mode, confidence, source).
    """
    # 1. experimentTypes (structured)
    exp_type_mode = None
    exp_types = project.get("experimentTypes", [])
    for et in exp_types:
        name = (et.get("name", "") if isinstance(et, dict) else str(et)).lower()
        if "data-independent" in name or "dia" == name:
            exp_type_mode = ("DIA", 0.95, "experimentTypes")
        elif "shotgun" in name or "data-dependent" in name:
            exp_type_mode = ("DDA", 0.90, "experimentTypes")
        elif "targeted" in name or "prm" in name or "srm" in name or "mrm" in name:
            exp_type_mode = ("PRM", 0.90, "experimentTypes")

    # 2. Protocol/description text
    # Check text first — it may give a more specific timsTOF mode (PASEF, diaPASEF)
    # that refines the generic experimentTypes classification
    text_fields = [
        project.get("projectDescription", ""),
        project.get("sampleProcessingProtocol", ""),
        project.get("dataProcessingProtocol", ""),
    ]
    combined_text = "\n".join(t or "" for t in text_fields)
    mode, conf, _ = _parse_protocol_text_for_mode(combined_text)
    if mode:
        return mode, conf, "protocol_text"

    # 3. Keywords
    keywords = project.get("keywords", [])
    kw_text = " ".join(str(k) for k in keywords).lower() if keywords else ""
    if kw_text:
        kw_mode, kw_conf, _ = _parse_protocol_text_for_mode(kw_text)
        if kw_mode:
            return kw_mode, min(kw_conf, 0.85), "keywords"

    # 4. Software inference from dataProcessingProtocol + software list
    software_text = (project.get("dataProcessingProtocol", "") or "").lower()
    software_list = project.get("softwares", [])
    for sw in software_list:
        name = sw.get("name", "") if isinstance(sw, dict) else str(sw)
        software_text += " " + name.lower()

    for sw_name in DIA_SOFTWARE:
        if sw_name in software_text:
            return "DIA", 0.65, f"software_inference:{sw_name}"
    for sw_name in DDA_SOFTWARE:
        if sw_name in software_text:
            return "DDA", 0.60, f"software_inference:{sw_name}"

    # 5. Title (last resort)
    title = project.get("title", "") or ""
    title_mode, title_conf, _ = _parse_protocol_text_for_mode(title)
    if title_mode:
        return title_mode, min(title_conf, 0.60), "title"

    # Fall back to experimentTypes if nothing else matched
    if exp_type_mode:
        return exp_type_mode

    return "unknown", 0.0, "none"


# ---------------------------------------------------------------------------
# Instrument classification
# ---------------------------------------------------------------------------

def classify_instrument(
    project: Dict[str, Any],
) -> Tuple[Optional[str], bool]:
    """Match instruments[] against known timsTOF variants.

    Returns (timstof_model, is_pure_timstof).
    is_pure_timstof is False when the project mixes timsTOF with other instruments.
    """
    instruments = project.get("instruments", [])
    if not instruments:
        return None, False

    timstof_model = None
    best_pattern_len = 0
    has_other = False

    for instr in instruments:
        name = instr.get("name", "") if isinstance(instr, dict) else str(instr)
        name_lower = name.lower().strip()

        matched = False
        # Check all patterns to find the most specific (longest) match
        for pattern, model in TIMSTOF_MODELS.items():
            if pattern in name_lower and len(pattern) > best_pattern_len:
                timstof_model = model
                best_pattern_len = len(pattern)
                matched = True

        if not matched:
            has_other = True

    is_pure = timstof_model is not None and not has_other
    return timstof_model, is_pure


# ---------------------------------------------------------------------------
# File manifest check
# ---------------------------------------------------------------------------

class FileCheckError(Exception):
    """The file listing could not be retrieved. NOT the same as "no files".

    The retired v2 ``files/byProject`` endpoint answers HTTP 200 with an empty
    body, which the previous version of this check read as "this dataset ships
    no Bruker data" and wrote into the catalog as a measurement of 0 GB. A
    failed fetch must never be representable as a size.
    """


def _fetch_file_listing(
    accession: str,
    cache_dir: Path,
    rate_limit: float,
) -> List[Dict[str, Any]]:
    """Every file record for a project, following v3 pagination.

    Raises FileCheckError if any page fails — a short read must not look like
    the end of the listing.
    """
    cache_file = cache_dir / f"files_v3_{accession}.json"
    cached = _load_cache(cache_file)
    if cached is not None:
        return cached

    files: List[Dict[str, Any]] = []
    page = 0
    while True:
        url = (
            f"{PRIDE_API_BASE}/projects/{accession}/files"
            f"?pageSize={FILE_PAGE_SIZE}&page={page}"
        )
        chunk = _api_get_with_retry(url)
        if chunk is None:
            raise FileCheckError(f"page {page} fetch failed")
        if not isinstance(chunk, list):
            raise FileCheckError(
                f"page {page}: expected a list, got {type(chunk).__name__}"
            )
        files.extend(chunk)
        time.sleep(rate_limit)
        if len(chunk) < FILE_PAGE_SIZE:
            break
        page += 1

    if not files:
        # v3 answers 200 with [] for an accession that is not in the database,
        # so an empty listing is ambiguous: genuinely no files, or a withdrawn
        # or mistyped accession. Confirm the project record before believing it.
        if _api_get_with_retry(f"{PRIDE_API_BASE}/projects/{accession}") is None:
            raise FileCheckError("empty listing and no project record")

    _save_cache(cache_file, files)
    return files


def validate_bruker_files(
    accession: str,
    cache_dir: Path,
    rate_limit: float = 0.5,
) -> Tuple[bool, int, float, float, str]:
    """Check if a project has Bruker raw data on PRIDE.

    Returns (has_files, n_files, size_gb, raw_size_gb, size_source); raises
    FileCheckError if the listing could not be fetched.

    Two views of the raw data, because neither alone is right: filename globs
    (``*.d.zip`` and friends) miss archive-wrapped raw — ``raw.zip``,
    ``*.tar.gz``, ``Topology.rar`` — while PRIDE's own RAW ``fileCategory``
    catches those. By this stage the catalog is timsTOF-only, so RAW here means
    Bruker .d in some wrapper. ``n_files`` is the union of both views and
    ``size_gb`` the better-supported estimate, so a glob miss no longer sinks a
    dataset at the ``--min-files`` gate or bins it as "small" when scoring.
    Same convention as ``audit_pool_size_license.py``'s ``size_source`` column.
    """
    files_data = _fetch_file_listing(accession, cache_dir, rate_limit)

    n_files = 0
    bruker_bytes = raw_bytes = 0
    for f in files_data:
        name = (f.get("fileName") or "").lower()
        size = f.get("fileSizeBytes", 0)
        is_glob = any(
            fnmatch.fnmatch(name, pat.lower()) for pat in BRUKER_FILE_PATTERNS
        )
        is_raw = ((f.get("fileCategory") or {}).get("value") or "").upper() == "RAW"
        if is_glob:
            bruker_bytes += size
        if is_raw:
            raw_bytes += size
        if is_glob or is_raw:
            n_files += 1

    size_gb = round(bruker_bytes / (1024**3), 2)
    raw_gb = round(raw_bytes / (1024**3), 2)
    if size_gb > 0:
        best_gb, source = size_gb, "filename_glob"
    elif raw_gb > 0:
        best_gb, source = raw_gb, "raw_category"
    else:
        best_gb, source = 0.0, "none"

    return bool(n_files), n_files, best_gb, raw_gb, source


# ---------------------------------------------------------------------------
# API + caching helpers
# ---------------------------------------------------------------------------

def _api_get_with_retry(
    url: str,
    timeout: int = 60,
    max_retries: int = 3,
) -> Optional[Any]:
    """GET with exponential backoff on 429/503."""
    for attempt in range(max_retries):
        try:
            resp = requests.get(url, timeout=timeout)
            if resp.status_code == 200:
                return resp.json()
            if resp.status_code in (429, 503):
                wait = 2 ** (attempt + 1)
                print(f"  Rate limited ({resp.status_code}), waiting {wait}s...")
                time.sleep(wait)
                continue
            resp.raise_for_status()
        except requests.RequestException as e:
            if attempt < max_retries - 1:
                wait = 2 ** (attempt + 1)
                print(f"  Request error: {e}, retrying in {wait}s...")
                time.sleep(wait)
            else:
                print(f"  Failed after {max_retries} attempts: {e}")
                return None
    return None


def _load_cache(path: Path, max_age_hours: float = 24.0) -> Optional[Any]:
    """Load cached JSON if it exists and is fresh enough."""
    if not path.exists():
        return None
    age = datetime.now() - datetime.fromtimestamp(path.stat().st_mtime)
    if age > timedelta(hours=max_age_hours):
        return None
    with open(path) as f:
        return json.load(f)


def _save_cache(path: Path, data: Any) -> None:
    """Save data as JSON cache file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f)


# ---------------------------------------------------------------------------
# Stage 1: Paginated PRIDE search
# ---------------------------------------------------------------------------

def discover_all_timstof_projects(
    keyword: str = "timsTOF",
    cache_dir: Path = DEFAULT_CACHE_DIR,
    rate_limit: float = 1.0,
    refresh: bool = False,
    page_size: int = 100,
) -> List[Dict[str, Any]]:
    """Paginate PRIDE search API for timsTOF projects.

    Returns list of raw project dicts from the API.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    all_projects = []
    page = 0

    while True:
        cache_file = cache_dir / f"search_page_{page}.json"

        page_data = None
        if not refresh:
            page_data = _load_cache(cache_file)

        if page_data is None:
            url = (
                f"{PRIDE_API_BASE}/search/projects"
                f"?keyword={keyword}&pageSize={page_size}&page={page}"
            )
            print(f"  Fetching page {page}...")
            page_data = _api_get_with_retry(url)

            if page_data is None:
                print(f"  Failed to fetch page {page}, stopping.")
                break

            _save_cache(cache_file, page_data)
            time.sleep(rate_limit)

        # PRIDE search API returns a list of project objects directly,
        # or a dict with _embedded.compactprojects
        projects = []
        if isinstance(page_data, list):
            projects = page_data
        elif isinstance(page_data, dict):
            # Try _embedded.compactprojects (HAL format)
            embedded = page_data.get("_embedded", {})
            projects = embedded.get("compactprojects", [])
            if not projects:
                # Maybe it's just a list under a different key
                projects = page_data.get("projects", [])

        if not projects:
            break

        all_projects.extend(projects)
        print(f"  Page {page}: {len(projects)} projects (total: {len(all_projects)})")

        # Check if there are more pages
        if len(projects) < page_size:
            break

        page += 1

    print(f"  Stage 1 complete: {len(all_projects)} candidate projects")
    return all_projects


# ---------------------------------------------------------------------------
# Stage 2: Instrument validation + classification
# ---------------------------------------------------------------------------

def _normalize_organism(name: str) -> str:
    """Normalize an organism name to a short key."""
    name_lower = name.lower().strip()
    for full_name, key in ORGANISM_ALIASES.items():
        if full_name in name_lower:
            return key
    # Return cleaned lowercase if no alias match
    return name_lower.replace(" ", "_")[:30]


def _extract_tissue(project: Dict[str, Any]) -> Optional[str]:
    """Try to extract tissue/cell type from project metadata."""
    # Check organism_parts field (search results use "organismsPart")
    parts = project.get("organismsPart", project.get("organism_parts", project.get("organismParts", [])))
    if parts:
        part = parts[0]
        name = part.get("name", "") if isinstance(part, dict) else str(part)
        if name:
            return name.lower()

    # Check diseases field as a proxy
    diseases = project.get("diseases", [])
    if diseases:
        d = diseases[0]
        name = d.get("name", "") if isinstance(d, dict) else str(d)
        if name and name.lower() not in ("not applicable", "not available"):
            return name.lower()

    return None


def classify_projects(
    projects: List[Dict[str, Any]],
    cache_dir: Path = DEFAULT_CACHE_DIR,
    rate_limit: float = 0.3,
    refresh: bool = False,
) -> List[Dict[str, Any]]:
    """Validate and classify each project.

    For compact search results, fetch full project details when needed.
    Returns list of classified dataset dicts (only timsTOF projects).
    """
    classified = []

    for i, proj in enumerate(projects):
        accession = proj.get("accession", "")
        if not accession:
            continue

        # The search API may return compact results; fetch full project if needed
        if "instruments" not in proj:
            full_proj = _fetch_full_project(accession, cache_dir, rate_limit, refresh)
            if full_proj is None:
                continue
            proj = full_proj

        # Instrument check
        timstof_model, is_pure = classify_instrument(proj)
        if timstof_model is None:
            continue

        # Acquisition mode
        acq_mode, acq_confidence, acq_source = classify_acquisition_mode(proj)

        # Organism
        organisms = proj.get("organisms", [])
        organism_names = []
        organism_keys = []
        for org in organisms:
            name = org.get("name", "") if isinstance(org, dict) else str(org)
            if name:
                organism_names.append(name)
                organism_keys.append(_normalize_organism(name))

        # Lab — search results use affiliations/labPIs, full projects use contacts
        lab = None
        submitter = None
        affiliations = proj.get("affiliations", [])
        if affiliations:
            lab = affiliations[0] if isinstance(affiliations[0], str) else affiliations[0].get("name")
        lab_pis = proj.get("labPIs", [])
        submitters = proj.get("submitters", [])
        if lab_pis:
            submitter = lab_pis[0] if isinstance(lab_pis[0], str) else lab_pis[0].get("name")
        elif submitters:
            submitter = submitters[0] if isinstance(submitters[0], str) else submitters[0].get("name")
        # Fallback to contacts (full project format)
        if not lab:
            contacts = proj.get("contacts", [])
            if contacts:
                lab = contacts[0].get("affiliation")
                if not submitter:
                    submitter = contacts[0].get("name")

        # Tissue
        tissue = _extract_tissue(proj)

        # Title & dates
        title = proj.get("title", "")
        pub_date = proj.get("publicationDate", proj.get("submissionDate", ""))

        # Warnings
        warnings = []
        if not is_pure:
            warnings.append("mixed_instruments")
        if len(organisms) > 1:
            warnings.append("multi_organism")
        if acq_mode == "unknown":
            warnings.append("unknown_acquisition_mode")

        entry = {
            "accession": accession,
            "title": title,
            "instrument_model": timstof_model,
            "is_pure_timstof": is_pure,
            "acquisition_mode": acq_mode,
            "acq_confidence": acq_confidence,
            "acq_source": acq_source,
            "organisms": organism_names,
            "organism_keys": organism_keys,
            "lab": lab,
            "submitter": submitter,
            "tissue": tissue,
            "publication_date": pub_date,
            "warnings": warnings,
            # Placeholders for Stage 3
            "has_bruker_files": None,
            "n_bruker_files": None,
            "bruker_size_gb": None,
            "raw_size_gb": None,
            "size_source": None,
        }
        classified.append(entry)

        if (i + 1) % 50 == 0:
            print(f"  Classified {i + 1}/{len(projects)}...")

    print(f"  Stage 2 complete: {len(classified)} timsTOF projects")
    return classified


def _fetch_full_project(
    accession: str,
    cache_dir: Path,
    rate_limit: float = 0.3,
    refresh: bool = False,
) -> Optional[Dict[str, Any]]:
    """Fetch full project details, with caching."""
    cache_file = cache_dir / f"project_{accession}.json"

    if not refresh:
        cached = _load_cache(cache_file)
        if cached is not None:
            return cached

    url = f"{PRIDE_API_BASE}/projects/{accession}"
    data = _api_get_with_retry(url)
    if data is not None:
        _save_cache(cache_file, data)
    time.sleep(rate_limit)
    return data


# ---------------------------------------------------------------------------
# Stage 3: File manifest check
# ---------------------------------------------------------------------------

def validate_all_files(
    datasets: List[Dict[str, Any]],
    cache_dir: Path = DEFAULT_CACHE_DIR,
    rate_limit: float = 0.5,
    min_files: int = 0,
) -> List[Dict[str, Any]]:
    """Check each dataset for Bruker .d files on PRIDE.

    Modifies datasets in-place and returns filtered list.
    """
    valid = []
    failed = []
    for i, ds in enumerate(datasets):
        acc = ds["accession"]
        try:
            has_files, n_files, size_gb, raw_gb, size_src = validate_bruker_files(
                acc, cache_dir, rate_limit
            )
        except FileCheckError as e:
            # Fields stay None so a later --cached-only re-score cannot read a
            # throttled fetch as "0 GB, no Bruker files". The dataset is held
            # back from `valid` because we genuinely do not know, and the count
            # is reported below rather than swallowed.
            ds["has_bruker_files"] = None
            ds["n_bruker_files"] = None
            ds["bruker_size_gb"] = None
            ds["raw_size_gb"] = None
            ds["size_source"] = None
            ds.setdefault("warnings", []).append(f"file_check_failed: {e}")
            failed.append(acc)
            continue

        ds["has_bruker_files"] = has_files
        ds["n_bruker_files"] = n_files
        ds["bruker_size_gb"] = size_gb
        ds["raw_size_gb"] = raw_gb
        ds["size_source"] = size_src

        if has_files and n_files >= min_files:
            valid.append(ds)

        if (i + 1) % 20 == 0:
            print(f"  Checked files for {i + 1}/{len(datasets)}...")

    print(f"  Stage 3 complete: {len(valid)}/{len(datasets)} have Bruker files")
    if failed:
        print(
            f"  WARNING: {len(failed)} file checks could not be completed and "
            f"were left unmeasured (not counted as empty): "
            f"{', '.join(failed[:5])}{' ...' if len(failed) > 5 else ''}"
        )
    return valid


# ---------------------------------------------------------------------------
# Diversity scoring (greedy set-cover)
# ---------------------------------------------------------------------------

def score_and_rank(
    datasets: List[Dict[str, Any]],
    weights: Optional[Dict[str, float]] = None,
) -> List[Dict[str, Any]]:
    """Greedy set-cover: iteratively pick dataset adding most diversity.

    Dimensions scored:
      - instrument_model, organism, acquisition_mode, lab, tissue, dataset_size

    Returns datasets sorted by priority_rank (1 = most diverse addition).
    """
    if not datasets:
        return []

    if weights is None:
        weights = DEFAULT_WEIGHTS

    # Track how often each dimension value has been seen
    seen_counts: Dict[str, Dict[str, int]] = {
        dim: {} for dim in weights
    }

    remaining = list(range(len(datasets)))
    ranked = []

    while remaining:
        best_idx = None
        best_score = -1.0

        for idx in remaining:
            ds = datasets[idx]
            score = _marginal_gain(ds, seen_counts, weights)
            if score > best_score:
                best_score = score
                best_idx = idx

        if best_idx is None:
            break

        # Select this dataset
        ds = datasets[best_idx]
        ranked.append(best_idx)
        remaining.remove(best_idx)

        # Update seen counts
        for dim in weights:
            val = _get_dimension_value(ds, dim)
            if val:
                seen_counts[dim][val] = seen_counts[dim].get(val, 0) + 1

    # Assign ranks
    for rank, idx in enumerate(ranked, 1):
        datasets[idx]["priority_rank"] = rank

    # Sort by rank
    return sorted(datasets, key=lambda d: d.get("priority_rank", 9999))


def _get_dimension_value(ds: Dict[str, Any], dim: str) -> Optional[str]:
    """Extract a single representative value for a diversity dimension."""
    if dim == "instrument_model":
        return ds.get("instrument_model")
    if dim == "organism":
        keys = ds.get("organism_keys", [])
        return keys[0] if keys else None
    if dim == "acquisition_mode":
        return ds.get("acquisition_mode")
    if dim == "lab":
        return ds.get("lab")
    if dim == "tissue":
        return ds.get("tissue")
    if dim == "dataset_size":
        n = ds.get("n_bruker_files")
        if n is None:
            return "unknown"
        if n <= 5:
            return "small"
        if n <= 20:
            return "medium"
        if n <= 100:
            return "large"
        return "very_large"
    return None


def _marginal_gain(
    ds: Dict[str, Any],
    seen_counts: Dict[str, Dict[str, int]],
    weights: Dict[str, float],
) -> float:
    """Compute marginal diversity gain for a candidate dataset."""
    total = 0.0
    for dim, weight in weights.items():
        val = _get_dimension_value(ds, dim)
        if val is None:
            continue
        count = seen_counts[dim].get(val, 0)
        if count == 0:
            gain = 1.0  # New value
        else:
            gain = 1.0 / math.sqrt(count + 1)  # Diminishing returns
        total += weight * gain
    return total


# ---------------------------------------------------------------------------
# Output generation
# ---------------------------------------------------------------------------

def generate_catalog(
    datasets: List[Dict[str, Any]],
    output_dir: Path,
) -> Path:
    """Write full catalog YAML."""
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "timstof_catalog.yaml"

    catalog = {
        "generated_at": datetime.now().isoformat(),
        "total_datasets": len(datasets),
        "datasets": datasets,
    }

    with open(path, "w") as f:
        yaml.dump(catalog, f, default_flow_style=False, sort_keys=False, width=120)

    print(f"Wrote catalog: {path} ({len(datasets)} datasets)")
    return path


def generate_csv(
    datasets: List[Dict[str, Any]],
    output_dir: Path,
) -> Path:
    """Write exhaustive CSV table of all candidate datasets."""
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "timstof_catalog.csv"

    fieldnames = [
        "priority_rank",
        "accession",
        "title",
        "instrument_model",
        "is_pure_timstof",
        "acquisition_mode",
        "acq_confidence",
        "acq_source",
        "organisms",
        "organism_keys",
        "lab",
        "submitter",
        "tissue",
        "publication_date",
        "has_bruker_files",
        "n_bruker_files",
        "bruker_size_gb",
        "raw_size_gb",
        "size_source",
        "warnings",
    ]

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for ds in datasets:
            row = dict(ds)
            # Flatten list fields to semicolon-separated strings
            row["organisms"] = "; ".join(row.get("organisms") or [])
            row["organism_keys"] = "; ".join(row.get("organism_keys") or [])
            row["warnings"] = "; ".join(row.get("warnings") or [])
            writer.writerow(row)

    print(f"Wrote CSV: {path} ({len(datasets)} rows)")
    return path


def generate_recommended(
    datasets: List[Dict[str, Any]],
    output_dir: Path,
    top_n: int = 20,
) -> Path:
    """Write config-compatible recommended datasets YAML."""
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "recommended_datasets.yaml"

    top = datasets[:top_n]

    # Build config-compatible format
    ds_list = []
    ds_metadata = {}

    for ds in top:
        acc = ds["accession"]
        ds_list.append(acc)

        org_keys = ds.get("organism_keys", [])
        organism = org_keys[0] if org_keys else "unknown"

        ds_metadata[acc] = {
            "description": ds.get("title", ""),
            "instrument": ds.get("instrument_model", ""),
            "acquisition_mode": ds.get("acquisition_mode", "unknown"),
            "organism": organism,
            "lab_id": ds.get("lab", ""),
            "priority_rank": ds.get("priority_rank"),
            "n_bruker_files": ds.get("n_bruker_files"),
            "bruker_size_gb": ds.get("bruker_size_gb"),
        }

        if ds.get("warnings"):
            ds_metadata[acc]["warnings"] = ds["warnings"]

    recommended = {
        "generated_at": datetime.now().isoformat(),
        "description": f"Top {len(top)} diversity-ranked timsTOF datasets from PRIDE",
        "datasets": ds_list,
        "dataset_metadata": ds_metadata,
    }

    with open(path, "w") as f:
        yaml.dump(recommended, f, default_flow_style=False, sort_keys=False, width=120)

    print(f"Wrote recommended: {path} ({len(top)} datasets)")
    return path


def generate_report(
    datasets: List[Dict[str, Any]],
    output_dir: Path,
    top_n: int = 20,
) -> Path:
    """Write human-readable discovery report."""
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "discovery_report.md"

    lines = [
        "# PRIDE timsTOF Discovery Report",
        f"",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"",
        f"## Summary",
        f"",
        f"- **Total timsTOF datasets:** {len(datasets)}",
    ]

    # Distribution counts
    instruments = {}
    modes = {}
    organisms = {}
    labs = {}
    for ds in datasets:
        m = ds.get("instrument_model", "unknown")
        instruments[m] = instruments.get(m, 0) + 1

        mode = ds.get("acquisition_mode", "unknown")
        modes[mode] = modes.get(mode, 0) + 1

        for ok in ds.get("organism_keys", []):
            organisms[ok] = organisms.get(ok, 0) + 1

        lab = ds.get("lab", "unknown")
        if lab:
            labs[lab] = labs.get(lab, 0) + 1

    lines.append(f"- **Instruments:** {len(instruments)} models")
    lines.append(f"- **Acquisition modes:** {len(modes)} types")
    lines.append(f"- **Organisms:** {len(organisms)} species")
    lines.append(f"- **Labs:** {len(labs)} institutions")

    # Instrument breakdown
    lines.extend(["", "## Instruments", ""])
    lines.append("| Model | Count |")
    lines.append("|-------|-------|")
    for m, c in sorted(instruments.items(), key=lambda x: -x[1]):
        lines.append(f"| {m} | {c} |")

    # Acquisition mode breakdown
    lines.extend(["", "## Acquisition Modes", ""])
    lines.append("| Mode | Count |")
    lines.append("|------|-------|")
    for m, c in sorted(modes.items(), key=lambda x: -x[1]):
        lines.append(f"| {m} | {c} |")

    # Organism breakdown
    lines.extend(["", "## Organisms", ""])
    lines.append("| Organism | Count |")
    lines.append("|----------|-------|")
    for o, c in sorted(organisms.items(), key=lambda x: -x[1]):
        lines.append(f"| {o} | {c} |")

    # Top-N recommended
    top = datasets[:top_n]
    lines.extend(["", f"## Top {len(top)} Recommended Datasets", ""])
    lines.append("| Rank | Accession | Instrument | Mode | Organism | Files | Size (GB) |")
    lines.append("|------|-----------|------------|------|----------|-------|-----------|")
    for ds in top:
        rank = ds.get("priority_rank", "-")
        acc = ds.get("accession", "")
        instr = ds.get("instrument_model", "")
        mode = ds.get("acquisition_mode", "unknown")
        orgs = ", ".join(ds.get("organism_keys", [])[:2])
        n = ds.get("n_bruker_files", "?")
        sz = ds.get("bruker_size_gb", "?")
        lines.append(f"| {rank} | {acc} | {instr} | {mode} | {orgs} | {n} | {sz} |")

    # Coverage analysis for top-N
    if top:
        top_instruments = set(d.get("instrument_model") for d in top)
        top_modes = set(d.get("acquisition_mode") for d in top)
        top_organisms = set()
        for d in top:
            top_organisms.update(d.get("organism_keys", []))

        lines.extend(["", "## Coverage Analysis", ""])
        lines.append(
            f"Top {len(top)} covers: "
            f"{len(top_instruments)}/{len(instruments)} instruments, "
            f"{len(top_modes)}/{len(modes)} modes, "
            f"{len(top_organisms)}/{len(organisms)} organisms"
        )

    # Datasets needing review
    review = [d for d in datasets if d.get("warnings")]
    if review:
        lines.extend(["", "## Datasets Needing Review", ""])
        lines.append("| Accession | Warnings |")
        lines.append("|-----------|----------|")
        for ds in review[:30]:
            acc = ds["accession"]
            warns = ", ".join(ds["warnings"])
            lines.append(f"| {acc} | {warns} |")
        if len(review) > 30:
            lines.append(f"| ... | {len(review) - 30} more |")

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"Wrote report: {path}")
    return path


# ---------------------------------------------------------------------------
# Config update
# ---------------------------------------------------------------------------

def update_config(
    config_path: Path,
    datasets: List[Dict[str, Any]],
    top_n: int = 10,
) -> None:
    """Append top-N datasets to an existing config YAML."""
    with open(config_path) as f:
        config = yaml.safe_load(f) or {}

    existing = set(config.get("datasets", []))
    metadata = config.get("dataset_metadata", {})

    added = 0
    for ds in datasets[:top_n]:
        acc = ds["accession"]
        if acc in existing:
            continue

        config.setdefault("datasets", []).append(acc)
        org_keys = ds.get("organism_keys", [])
        organism = org_keys[0] if org_keys else "unknown"

        metadata[acc] = {
            "description": ds.get("title", "")[:120],
            "instrument": ds.get("instrument_model", ""),
            "acquisition_mode": ds.get("acquisition_mode", "unknown"),
            "lab_id": ds.get("lab", ""),
        }
        added += 1

    config["dataset_metadata"] = metadata

    with open(config_path, "w") as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False, width=120)

    print(f"Updated {config_path}: added {added} new datasets")


# ---------------------------------------------------------------------------
# Blacklist
# ---------------------------------------------------------------------------

def load_blacklist(path: Optional[Path] = None) -> set:
    """Load blacklisted accessions."""
    if path is None:
        path = Path("config/blacklist.yaml")
    if not path.exists():
        return set()
    with open(path) as f:
        data = yaml.safe_load(f) or {}
    bl = data.get("blacklist", {})
    if bl is None:
        return set()
    return set(bl.keys())


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Discover and rank timsTOF datasets on PRIDE for San José",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/discover_pride.py                          # Full discovery
  python scripts/discover_pride.py --cached-only            # Re-score from cache
  python scripts/discover_pride.py --top 20                 # Show top 20
  python scripts/discover_pride.py --mode DDA --top 50      # Filter to DDA
  python scripts/discover_pride.py --skip-file-check        # Skip Stage 3
  python scripts/discover_pride.py --list-top 5             # Accessions only
  python scripts/discover_pride.py --update-config config/config.yaml
  python scripts/discover_pride.py --incremental            # Only new since last run
""",
    )

    parser.add_argument(
        "--output", "-o", type=Path, default=DEFAULT_OUTPUT_DIR,
        help="Output directory (default: data/discovery/)",
    )
    parser.add_argument(
        "--keyword", default="timsTOF",
        help="PRIDE search keyword (default: timsTOF)",
    )
    parser.add_argument(
        "--mode", choices=["DDA", "DIA", "PASEF", "diaPASEF", "PRM", "MALDI", "single-cell"],
        help="Filter results to specific acquisition mode",
    )
    parser.add_argument(
        "--organism", help="Filter results to specific organism key (e.g. human, yeast)",
    )
    parser.add_argument(
        "--instrument", help="Filter results to specific instrument model substring",
    )
    parser.add_argument(
        "--min-files", type=int, default=1,
        help="Minimum number of Bruker files required (default: 1)",
    )
    parser.add_argument(
        "--top", type=int, default=20,
        help="Number of top-ranked datasets to include in recommendations (default: 20)",
    )
    parser.add_argument(
        "--list-top", type=int, metavar="N",
        help="Print top N accessions and exit (scriptable output)",
    )
    parser.add_argument(
        "--exclude-processed", nargs="*", metavar="ACC",
        help="Exclude already-processed accessions from ranking",
    )
    parser.add_argument(
        "--cached-only", action="store_true",
        help="Only use cached API responses; do not make new requests",
    )
    parser.add_argument(
        "--refresh", action="store_true",
        help="Force re-fetch all cached data",
    )
    parser.add_argument(
        "--skip-file-check", action="store_true",
        help="Skip Stage 3 (Bruker file validation) for faster results",
    )
    parser.add_argument(
        "--rate-limit", type=float, default=1.0,
        help="Seconds between search page requests (default: 1.0)",
    )
    parser.add_argument(
        "--update-config", type=Path, metavar="PATH",
        help="Auto-append top-N datasets to this config YAML",
    )
    parser.add_argument(
        "--csv", action="store_true",
        help="Export exhaustive CSV table of all candidates",
    )
    parser.add_argument(
        "--incremental", action="store_true",
        help="Only discover datasets newer than the last run",
    )
    parser.add_argument(
        "--blacklist", type=Path, default=None,
        help="Path to blacklist YAML (default: config/blacklist.yaml)",
    )

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    cache_dir = args.output / "cache"
    output_dir = args.output

    # Load blacklist
    blacklisted = load_blacklist(args.blacklist)
    if blacklisted:
        print(f"Loaded {len(blacklisted)} blacklisted accessions")

    # Load previous catalog for incremental mode
    prev_accessions: set = set()
    if args.incremental:
        prev_catalog = output_dir / "timstof_catalog.yaml"
        if prev_catalog.exists():
            with open(prev_catalog) as f:
                prev = yaml.safe_load(f) or {}
            prev_accessions = {d["accession"] for d in prev.get("datasets", [])}
            print(f"Incremental mode: {len(prev_accessions)} previously discovered")

    # ---- Stage 1: Search ----
    print("\n=== Stage 1: PRIDE Search ===")
    if args.cached_only:
        # Load from existing cache files
        projects = []
        page = 0
        while True:
            cache_file = cache_dir / f"search_page_{page}.json"
            if not cache_file.exists():
                break
            with open(cache_file) as f:
                page_data = json.load(f)
            if isinstance(page_data, list):
                projects.extend(page_data)
            elif isinstance(page_data, dict):
                embedded = page_data.get("_embedded", {})
                p = embedded.get("compactprojects", [])
                if not p:
                    p = page_data.get("projects", [])
                projects.extend(p)
            page += 1
        print(f"  Loaded {len(projects)} projects from cache ({page} pages)")
    else:
        projects = discover_all_timstof_projects(
            keyword=args.keyword,
            cache_dir=cache_dir,
            rate_limit=args.rate_limit,
            refresh=args.refresh,
        )

    if not projects:
        print("No projects found. Check your network connection or try --cached-only.")
        return 1

    # ---- Stage 2: Classify ----
    print("\n=== Stage 2: Instrument Validation & Classification ===")
    datasets = classify_projects(
        projects, cache_dir=cache_dir,
        rate_limit=args.rate_limit * 0.3,
        refresh=args.refresh,
    )

    # Remove blacklisted
    before = len(datasets)
    datasets = [d for d in datasets if d["accession"] not in blacklisted]
    if before != len(datasets):
        print(f"  Removed {before - len(datasets)} blacklisted datasets")

    # Remove already-known in incremental mode
    if prev_accessions:
        before = len(datasets)
        datasets = [d for d in datasets if d["accession"] not in prev_accessions]
        print(f"  Incremental: {before - len(datasets)} already known, {len(datasets)} new")

    # Remove already-processed if specified
    if args.exclude_processed:
        excl = set(args.exclude_processed)
        datasets = [d for d in datasets if d["accession"] not in excl]

    # ---- Stage 3: File check (optional) ----
    if not args.skip_file_check and not args.cached_only:
        print("\n=== Stage 3: Bruker File Validation ===")
        datasets = validate_all_files(
            datasets, cache_dir=cache_dir,
            rate_limit=args.rate_limit * 0.5,
            min_files=args.min_files,
        )

    # ---- Apply filters ----
    if args.mode:
        datasets = [d for d in datasets if d.get("acquisition_mode") == args.mode]
        print(f"  Filter --mode={args.mode}: {len(datasets)} datasets")

    if args.organism:
        datasets = [
            d for d in datasets
            if args.organism in d.get("organism_keys", [])
        ]
        print(f"  Filter --organism={args.organism}: {len(datasets)} datasets")

    if args.instrument:
        datasets = [
            d for d in datasets
            if args.instrument.lower() in (d.get("instrument_model") or "").lower()
        ]
        print(f"  Filter --instrument={args.instrument}: {len(datasets)} datasets")

    if not datasets:
        print("\nNo datasets match the specified filters.")
        return 1

    # ---- Rank ----
    print(f"\n=== Diversity Ranking ({len(datasets)} datasets) ===")
    datasets = score_and_rank(datasets)

    # ---- Output ----
    if args.csv and not args.list_top:
        generate_csv(datasets, output_dir)
        return 0

    if args.list_top:
        for ds in datasets[: args.list_top]:
            print(ds["accession"])
        return 0

    print(f"\n=== Generating Outputs ===")
    generate_catalog(datasets, output_dir)
    generate_csv(datasets, output_dir)
    generate_recommended(datasets, output_dir, top_n=args.top)
    generate_report(datasets, output_dir, top_n=args.top)

    if args.update_config:
        update_config(args.update_config, datasets, top_n=args.top)

    # ---- Summary ----
    print(f"\n=== Top {min(args.top, len(datasets))} Datasets ===")
    print(f"{'Rank':<5} {'Accession':<12} {'Instrument':<16} {'Mode':<12} {'Organism':<12} {'Files':>5} {'GB':>7}")
    print("-" * 75)
    for ds in datasets[: args.top]:
        rank = ds.get("priority_rank", "-")
        acc = ds["accession"]
        instr = (ds.get("instrument_model") or "")[:15]
        mode = (ds.get("acquisition_mode") or "unknown")[:11]
        orgs = ", ".join(ds.get("organism_keys", [])[:1]) or "?"
        n = ds.get("n_bruker_files", "?")
        sz = ds.get("bruker_size_gb", "?")
        print(f"{rank:<5} {acc:<12} {instr:<16} {mode:<12} {orgs:<12} {str(n):>5} {str(sz):>7}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
