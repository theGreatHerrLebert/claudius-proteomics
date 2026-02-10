#!/usr/bin/env python3
"""
Sample Group Resolver for Multi-Organism PRIDE Datasets

Resolves per-file organism + enzyme assignments from PRIDE archive file names.
The PRIDE REST API provides no per-file species annotation — organisms are listed
only at project level. File-to-organism mapping must be inferred from file names.

Vocabulary:
  - Project: PRIDE submission (PXD accession)
  - Sample Group: Collection of runs sharing (organism, enzyme). Unit of search execution.
  - Run: Single .d folder (one MS acquisition).

Usage:
    python scripts/sample_group_resolver.py PXD019086 --raw-dir data/raw/PXD019086
"""

import json
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml

# ---------------------------------------------------------------------------
# Lookup tables
# ---------------------------------------------------------------------------

ORGANISM_ALIASES: Dict[str, str] = {
    "hela": "human",
    "human": "human",
    "sapiens": "human",
    "ecoli": "ecoli",
    "e.coli": "ecoli",
    "coli": "ecoli",
    "yeast": "yeast",
    "cerevisiae": "yeast",
    "drosophila": "drosophila",
    "melanogaster": "drosophila",
    "celegans": "c_elegans",
    "c.elegans": "c_elegans",
    "elegans": "c_elegans",
    "mouse": "mouse",
    "proteometools": "human",
}

ENZYME_ALIASES: Dict[str, str] = {
    "trp": "trypsin",
    "trypsin": "trypsin",
    "tryp": "trypsin",
    "lysc": "lysc",
    "lysn": "lysn",
}

TAXON_TO_ORGANISM: Dict[int, str] = {
    9606: "human",
    562: "ecoli",
    83333: "ecoli",
    4932: "yeast",
    559292: "yeast",
    7227: "drosophila",
    6239: "c_elegans",
    10090: "mouse",
}

ORGANISM_NAMES: Dict[str, str] = {
    "human": "Homo sapiens",
    "ecoli": "Escherichia coli",
    "yeast": "Saccharomyces cerevisiae",
    "drosophila": "Drosophila melanogaster",
    "c_elegans": "Caenorhabditis elegans",
    "mouse": "Mus musculus",
}

ORGANISM_TAXONS: Dict[str, int] = {
    "human": 9606,
    "ecoli": 562,
    "yeast": 559292,
    "drosophila": 7227,
    "c_elegans": 6239,
    "mouse": 10090,
}


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class SampleGroup:
    """A collection of runs sharing the same organism and enzyme."""

    group_id: str           # e.g. "human_trypsin"
    organism_key: str       # config key: "human"
    organism_name: str      # "Homo sapiens"
    taxon_id: int           # 9606
    enzyme: str             # "trypsin" | "lysc" | "lysn"
    sample_type: str        # "biological" | "synthetic"
    source_archive: str     # "Raw_HeLa_Trp.zip"
    runs: List[str] = field(default_factory=list)  # .d folder names

    @property
    def n_runs(self) -> int:
        return len(self.runs)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "group_id": self.group_id,
            "organism_key": self.organism_key,
            "organism_name": self.organism_name,
            "taxon_id": self.taxon_id,
            "enzyme": self.enzyme,
            "sample_type": self.sample_type,
            "source_archive": self.source_archive,
            "n_runs": self.n_runs,
            "runs": self.runs,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "SampleGroup":
        return cls(
            group_id=d["group_id"],
            organism_key=d["organism_key"],
            organism_name=d.get("organism_name", ""),
            taxon_id=d.get("taxon_id", 0),
            enzyme=d["enzyme"],
            sample_type=d.get("sample_type", "biological"),
            source_archive=d.get("source_archive", ""),
            runs=d.get("runs", []),
        )


@dataclass
class SampleGroupManifest:
    """Complete sample group resolution for a PRIDE accession."""

    accession: str
    groups: List[SampleGroup] = field(default_factory=list)
    unassigned_runs: List[str] = field(default_factory=list)
    generated_at: str = field(default_factory=lambda: datetime.now().isoformat())

    @property
    def is_multi_organism(self) -> bool:
        organism_keys = {g.organism_key for g in self.groups}
        return len(organism_keys) > 1

    @property
    def run_to_group(self) -> Dict[str, str]:
        """Map .d folder name -> group_id."""
        mapping = {}
        for group in self.groups:
            for run in group.runs:
                mapping[run] = group.group_id
        return mapping

    @property
    def organisms(self) -> List[Dict[str, Any]]:
        """Unique organisms across all groups."""
        seen = set()
        result = []
        for g in self.groups:
            if g.organism_key not in seen:
                seen.add(g.organism_key)
                result.append({
                    "name": g.organism_name,
                    "taxon_id": g.taxon_id,
                    "key": g.organism_key,
                })
        return result

    def to_dict(self) -> Dict[str, Any]:
        return {
            "accession": self.accession,
            "generated_at": self.generated_at,
            "is_multi_organism": self.is_multi_organism,
            "n_groups": len(self.groups),
            "organisms": self.organisms,
            "groups": [g.to_dict() for g in self.groups],
            "unassigned_runs": self.unassigned_runs,
        }

    def to_yaml(self, path: Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False, sort_keys=False)

    @classmethod
    def from_yaml(cls, path: Path) -> "SampleGroupManifest":
        with open(path) as f:
            data = yaml.safe_load(f)
        return cls.from_dict(data)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "SampleGroupManifest":
        manifest = cls(
            accession=data["accession"],
            generated_at=data.get("generated_at", ""),
            unassigned_runs=data.get("unassigned_runs", []),
        )
        for gd in data.get("groups", []):
            manifest.groups.append(SampleGroup.from_dict(gd))
        return manifest


# ---------------------------------------------------------------------------
# File name heuristic
# ---------------------------------------------------------------------------

def _parse_archive_name(filename: str) -> Tuple[Optional[str], Optional[str], str]:
    """Extract organism_key and enzyme from an archive file name.

    Returns:
        (organism_key, enzyme, sample_type)
    """
    # Normalize: lowercase, remove extension, replace separators
    stem = Path(filename).stem
    if stem.endswith(".d"):
        stem = stem[:-2]
    normalized = stem.lower().replace("-", "_").replace(".", "_")
    tokens = re.split(r"[_\s]+", normalized)

    organism_key = None
    enzyme = None
    sample_type = "biological"

    # Check for synthetic libraries (ProteomeTools)
    if "proteometools" in normalized or "pt_" in normalized:
        sample_type = "synthetic"
        organism_key = "human"
        enzyme = "trypsin"  # ProteomeTools peptides are searched with trypsin rules

    # Scan tokens for organism
    for token in tokens:
        if token in ORGANISM_ALIASES:
            organism_key = ORGANISM_ALIASES[token]
            break

    # If no hit from tokens, try substring matching on the full normalized name
    if organism_key is None:
        for alias, key in sorted(ORGANISM_ALIASES.items(), key=lambda x: -len(x[0])):
            if alias in normalized:
                organism_key = key
                break

    # Scan tokens for enzyme
    for token in tokens:
        if token in ENZYME_ALIASES:
            enzyme = ENZYME_ALIASES[token]
            break

    # Substring fallback for enzyme
    if enzyme is None:
        for alias, enz in sorted(ENZYME_ALIASES.items(), key=lambda x: -len(x[0])):
            if alias in normalized:
                enzyme = enz
                break

    return organism_key, enzyme, sample_type


def _make_group_id(organism_key: str, enzyme: str, sample_type: str) -> str:
    """Build a unique group_id string."""
    if sample_type == "synthetic":
        return f"{organism_key}_{enzyme}_synthetic"
    return f"{organism_key}_{enzyme}"


# ---------------------------------------------------------------------------
# PRIDE API helpers
# ---------------------------------------------------------------------------

PRIDE_API_BASE = "https://www.ebi.ac.uk/pride/ws/archive/v2"


def _fetch_pride_file_list(accession: str) -> List[Dict[str, Any]]:
    """Fetch the file list from PRIDE REST API."""
    import requests

    url = f"{PRIDE_API_BASE}/files/byProject?accession={accession}"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    return response.json()


def _fetch_pride_organisms(accession: str) -> List[Dict[str, Any]]:
    """Fetch organisms from PRIDE project metadata."""
    import requests

    url = f"{PRIDE_API_BASE}/projects/{accession}"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    project = response.json()
    return project.get("organisms", [])


def _get_raw_archives(file_list: List[Dict[str, Any]]) -> List[str]:
    """Filter PRIDE file list for RAW category zip archives."""
    raw_files = []
    for f in file_list:
        category = f.get("fileCategory", {})
        cat_value = category.get("value", "") if isinstance(category, dict) else str(category)
        filename = f.get("fileName", "")
        if cat_value.upper() == "RAW" and filename.endswith(".zip"):
            raw_files.append(filename)
    return sorted(raw_files)


# ---------------------------------------------------------------------------
# Core resolution
# ---------------------------------------------------------------------------

def resolve_from_file_list(
    accession: str,
    archive_names: List[str],
    project_organisms: Optional[List[Dict[str, Any]]] = None,
) -> SampleGroupManifest:
    """Resolve sample groups from a list of archive file names.

    This does NOT require downloaded data — only the archive names (from PRIDE API
    or from directory listing).

    Args:
        accession: PRIDE accession
        archive_names: List of zip file names (e.g. "Raw_HeLa_Trp.zip")
        project_organisms: Optional list of project-level organisms from PRIDE API

    Returns:
        SampleGroupManifest with groups (runs will be empty until populated)
    """
    # Build project-level organism set for cross-validation
    project_organism_keys = set()
    if project_organisms:
        for org in project_organisms:
            accession_str = org.get("accession", "")
            if accession_str:
                try:
                    taxon_id = int(accession_str)
                    if taxon_id in TAXON_TO_ORGANISM:
                        project_organism_keys.add(TAXON_TO_ORGANISM[taxon_id])
                except ValueError:
                    pass
            name = org.get("name", "").lower()
            for alias, key in ORGANISM_ALIASES.items():
                if alias in name:
                    project_organism_keys.add(key)
                    break

    # Parse each archive
    groups_by_id: Dict[str, SampleGroup] = {}
    for archive in archive_names:
        organism_key, enzyme, sample_type = _parse_archive_name(archive)

        if organism_key is None or enzyme is None:
            continue

        group_id = _make_group_id(organism_key, enzyme, sample_type)
        if group_id not in groups_by_id:
            groups_by_id[group_id] = SampleGroup(
                group_id=group_id,
                organism_key=organism_key,
                organism_name=ORGANISM_NAMES.get(organism_key, organism_key),
                taxon_id=ORGANISM_TAXONS.get(organism_key, 0),
                enzyme=enzyme,
                sample_type=sample_type,
                source_archive=archive,
            )
        else:
            # Multiple archives can map to same group; keep first as source_archive
            pass

    # Cross-validate if project organisms available
    if project_organism_keys:
        resolved_keys = {g.organism_key for g in groups_by_id.values()}
        extra = resolved_keys - project_organism_keys
        if extra:
            print(f"  Warning: resolved organisms {extra} not in project metadata")

    manifest = SampleGroupManifest(
        accession=accession,
        groups=sorted(groups_by_id.values(), key=lambda g: g.group_id),
    )
    return manifest


def populate_runs(manifest: SampleGroupManifest, raw_dir: Path) -> None:
    """Scan raw_dir for .d folders and assign each to a sample group.

    Modifies the manifest in-place: populates group.runs and manifest.unassigned_runs.
    """
    if not raw_dir.exists():
        return

    d_folders = sorted(p.name for p in raw_dir.iterdir() if p.is_dir() and p.name.endswith(".d"))
    assigned = set()

    for d_name in d_folders:
        # Try to match .d name to a group via the same heuristic
        organism_key, enzyme, sample_type = _parse_archive_name(d_name)
        if organism_key is not None and enzyme is not None:
            group_id = _make_group_id(organism_key, enzyme, sample_type)
            for group in manifest.groups:
                if group.group_id == group_id:
                    group.runs.append(d_name)
                    assigned.add(d_name)
                    break
            else:
                # No group with this ID — check if organism matches any group
                for group in manifest.groups:
                    if group.organism_key == organism_key and group.enzyme == enzyme:
                        group.runs.append(d_name)
                        assigned.add(d_name)
                        break

    manifest.unassigned_runs = [d for d in d_folders if d not in assigned]


def resolve_sample_groups(
    accession: str,
    raw_dir: Path,
    metadata_dir: Path,
    config: Dict[str, Any],
) -> SampleGroupManifest:
    """Full resolution: PRIDE API + file name heuristic + .d folder scan.

    Args:
        accession: PRIDE accession
        raw_dir: Directory containing .d folders
        metadata_dir: Directory for metadata files
        config: Pipeline configuration

    Returns:
        SampleGroupManifest with populated groups
    """
    # 1. Try to get archive names from PRIDE API
    archive_names: List[str] = []
    project_organisms: List[Dict[str, Any]] = []
    try:
        file_list = _fetch_pride_file_list(accession)
        archive_names = _get_raw_archives(file_list)
        project_organisms = _fetch_pride_organisms(accession)
        print(f"  PRIDE API: {len(archive_names)} raw archives, {len(project_organisms)} organisms")
    except Exception as e:
        print(f"  Warning: Could not fetch PRIDE file list: {e}")

    # 2. If no PRIDE data, fall back to scanning raw_dir for .d folders directly
    if not archive_names and raw_dir.exists():
        d_folders = sorted(
            p.name for p in raw_dir.iterdir() if p.is_dir() and p.name.endswith(".d")
        )
        archive_names = d_folders
        print(f"  Fallback: inferring groups from {len(d_folders)} .d folders")

    # 3. Resolve groups from archive names
    manifest = resolve_from_file_list(accession, archive_names, project_organisms)

    # 4. If we got no groups but have a single-organism config, create a default group
    if not manifest.groups:
        organism_key = config.get("dataset_metadata", {}).get(accession, {}).get("organism")
        if organism_key:
            default_group = SampleGroup(
                group_id=f"{organism_key}_trypsin",
                organism_key=organism_key,
                organism_name=ORGANISM_NAMES.get(organism_key, organism_key),
                taxon_id=ORGANISM_TAXONS.get(organism_key, 0),
                enzyme="trypsin",
                sample_type="biological",
                source_archive="",
            )
            manifest.groups.append(default_group)
            print(f"  Fallback: created default group {default_group.group_id}")

    # 5. Populate runs from .d folders
    populate_runs(manifest, raw_dir)

    # 6. If some .d files are unassigned and there's only one group, put them there
    if manifest.unassigned_runs and len(manifest.groups) == 1:
        manifest.groups[0].runs.extend(manifest.unassigned_runs)
        manifest.unassigned_runs = []

    return manifest


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Resolve sample groups for a multi-organism PRIDE dataset"
    )
    parser.add_argument("accession", help="PRIDE accession (e.g. PXD019086)")
    parser.add_argument("--raw-dir", type=Path, help="Directory containing .d folders")
    parser.add_argument("--output", "-o", type=Path, help="Output YAML path")
    parser.add_argument(
        "--config", "-c", default="config/config.yaml", help="Config file"
    )

    args = parser.parse_args()

    # Load config
    config: Dict[str, Any] = {}
    config_path = Path(args.config)
    if config_path.exists():
        with open(config_path) as f:
            config = yaml.safe_load(f) or {}

    raw_dir = args.raw_dir or Path(f"data/raw/{args.accession}")
    metadata_dir = Path(f"data/metadata/{args.accession}")

    print(f"Resolving sample groups for {args.accession}")
    manifest = resolve_sample_groups(args.accession, raw_dir, metadata_dir, config)

    # Output
    output_path = args.output or metadata_dir / "sample_groups.yaml"
    manifest.to_yaml(output_path)
    print(f"\nWrote {output_path}")

    # Summary
    print(f"\nSample Group Summary:")
    print(f"  Accession: {manifest.accession}")
    print(f"  Multi-organism: {manifest.is_multi_organism}")
    print(f"  Groups: {len(manifest.groups)}")
    for g in manifest.groups:
        print(f"    {g.group_id}: {g.organism_name} / {g.enzyme} "
              f"({g.sample_type}) — {g.n_runs} runs [{g.source_archive}]")
    if manifest.unassigned_runs:
        print(f"  Unassigned runs: {len(manifest.unassigned_runs)}")
        for r in manifest.unassigned_runs[:5]:
            print(f"    {r}")
        if len(manifest.unassigned_runs) > 5:
            print(f"    ... and {len(manifest.unassigned_runs) - 5} more")


if __name__ == "__main__":
    main()
