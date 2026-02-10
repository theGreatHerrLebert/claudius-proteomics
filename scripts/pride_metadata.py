#!/usr/bin/env python3
"""
PRIDE API metadata extraction for San José.

Fetches dataset metadata from PRIDE REST API and parses protocol text
to extract bias-aware fields like gradient_length, column_type, acquisition_mode.
"""

import json
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import requests
import yaml

# PRIDE REST API v2 base URL
PRIDE_API_BASE = "https://www.ebi.ac.uk/pride/ws/archive/v2"


@dataclass
class MetadataField:
    """
    A single metadata field with provenance tracking.

    Attributes:
        value: The field value (any type)
        status: "auto" | "inferred" | "manual" | "missing"
        source: Where the value came from (e.g., "pride_api.organisms[0].name")
        confidence: Confidence score for inferred values (0.0-1.0)
        hint: Suggestion for manual entry if missing
    """

    value: Any
    status: str  # "auto" | "inferred" | "manual" | "missing"
    source: Optional[str] = None
    confidence: float = 1.0
    hint: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for YAML serialization."""
        result = {"value": self.value, "status": self.status}
        if self.source:
            result["source"] = self.source
        if self.status == "inferred" and self.confidence < 1.0:
            result["confidence"] = self.confidence
        if self.hint:
            result["hint"] = self.hint
        return result

    @classmethod
    def auto(cls, value: Any, source: str) -> "MetadataField":
        """Create an auto-extracted field."""
        return cls(value=value, status="auto", source=source)

    @classmethod
    def inferred(
        cls, value: Any, source: str, confidence: float = 0.8
    ) -> "MetadataField":
        """Create an inferred field (parsed from text)."""
        return cls(
            value=value, status="inferred", source=source, confidence=confidence
        )

    @classmethod
    def manual(cls, value: Any) -> "MetadataField":
        """Create a manually-entered field."""
        return cls(value=value, status="manual", source="manual_entry")

    @classmethod
    def missing(cls, hint: Optional[str] = None) -> "MetadataField":
        """Create a missing field placeholder."""
        return cls(value=None, status="missing", hint=hint)


@dataclass
class DatasetMetadata:
    """
    Complete metadata for a PRIDE dataset with provenance tracking.

    All fields use MetadataField to track extraction status and source.
    """

    accession: str

    # Core identifiers
    dataset_id: MetadataField = field(default_factory=lambda: MetadataField.missing())
    lab_id: MetadataField = field(default_factory=lambda: MetadataField.missing())

    # Scientific context
    organism: MetadataField = field(default_factory=lambda: MetadataField.missing())
    organisms_all: MetadataField = field(default_factory=lambda: MetadataField.missing())
    is_multi_organism: bool = False
    instrument: MetadataField = field(default_factory=lambda: MetadataField.missing())
    title: MetadataField = field(default_factory=lambda: MetadataField.missing())
    publication_doi: MetadataField = field(default_factory=lambda: MetadataField.missing())

    # Bias-aware fields (often need manual entry or raw data extraction)
    gradient_length: MetadataField = field(
        default_factory=lambda: MetadataField.missing(
            hint="Check publication methods section or raw data LC method"
        )
    )
    column_type: MetadataField = field(
        default_factory=lambda: MetadataField.missing(
            hint="Check publication methods section for column specifications"
        )
    )
    acquisition_mode: MetadataField = field(
        default_factory=lambda: MetadataField.missing(
            hint="Check if DDA, DIA, PASEF, diaPASEF, PRM, etc."
        )
    )
    lc_system: MetadataField = field(
        default_factory=lambda: MetadataField.missing(
            hint="Check raw data chromatography info or publication"
        )
    )

    # File information
    num_files: MetadataField = field(default_factory=lambda: MetadataField.missing())
    file_types: MetadataField = field(default_factory=lambda: MetadataField.missing())
    total_size_gb: MetadataField = field(default_factory=lambda: MetadataField.missing())

    # Generation metadata
    generated_at: str = field(default_factory=lambda: datetime.now().isoformat())
    version: str = "1.0"

    # Raw API response (for debugging)
    _raw_response: Optional[Dict] = field(default=None, repr=False)

    def is_complete(self) -> bool:
        """Check if all required fields are populated (not missing)."""
        required_fields = [
            self.dataset_id,
            self.lab_id,
            self.organism,
            self.instrument,
            self.gradient_length,
            self.acquisition_mode,
        ]
        return all(f.status != "missing" for f in required_fields)

    def validation_summary(self) -> Dict[str, Any]:
        """Get a summary of field validation status."""
        all_fields = {
            "dataset_id": self.dataset_id,
            "lab_id": self.lab_id,
            "organism": self.organism,
            "instrument": self.instrument,
            "title": self.title,
            "publication_doi": self.publication_doi,
            "gradient_length": self.gradient_length,
            "column_type": self.column_type,
            "acquisition_mode": self.acquisition_mode,
            "lc_system": self.lc_system,
        }

        auto_count = sum(1 for f in all_fields.values() if f.status == "auto")
        inferred_count = sum(1 for f in all_fields.values() if f.status == "inferred")
        manual_count = sum(1 for f in all_fields.values() if f.status == "manual")
        missing_count = sum(1 for f in all_fields.values() if f.status == "missing")
        missing_fields = [k for k, f in all_fields.items() if f.status == "missing"]

        return {
            "complete": self.is_complete(),
            "auto_fields": auto_count,
            "inferred_fields": inferred_count,
            "manual_fields": manual_count,
            "missing_fields": missing_count,
            "missing_field_names": missing_fields,
            "requires_review": missing_count > 0 or inferred_count > 0,
        }

    def to_yaml_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for YAML serialization."""
        return {
            "accession": self.accession,
            "version": self.version,
            "generated_at": self.generated_at,
            "title": self.title.to_dict() if self.title.value else None,
            "publication_doi": self.publication_doi.to_dict()
            if self.publication_doi.value
            else None,
            "is_multi_organism": self.is_multi_organism,
            "fields": {
                "dataset_id": self.dataset_id.to_dict(),
                "lab_id": self.lab_id.to_dict(),
                "organism": self.organism.to_dict(),
                "organisms_all": self.organisms_all.to_dict()
                if self.organisms_all.value
                else None,
                "instrument": self.instrument.to_dict(),
                "gradient_length": self.gradient_length.to_dict(),
                "column_type": self.column_type.to_dict(),
                "acquisition_mode": self.acquisition_mode.to_dict(),
                "lc_system": self.lc_system.to_dict(),
            },
            "files": {
                "num_files": self.num_files.to_dict() if self.num_files.value else None,
                "file_types": self.file_types.to_dict()
                if self.file_types.value
                else None,
                "total_size_gb": self.total_size_gb.to_dict()
                if self.total_size_gb.value
                else None,
            },
            "validation": self.validation_summary(),
        }

    def to_yaml(self, path: Path) -> None:
        """Write metadata to YAML file."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        yaml_dict = self.to_yaml_dict()

        # Custom YAML dumper for nice formatting
        def str_representer(dumper, data):
            if "\n" in data:
                return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")
            return dumper.represent_scalar("tag:yaml.org,2002:str", data)

        yaml.add_representer(str, str_representer)

        with open(path, "w") as f:
            yaml.dump(yaml_dict, f, default_flow_style=False, sort_keys=False)

    @classmethod
    def from_pride_api(cls, accession: str) -> "DatasetMetadata":
        """
        Fetch metadata from PRIDE REST API.

        Args:
            accession: PRIDE accession (e.g., "PXD019086")

        Returns:
            DatasetMetadata with auto-extracted fields
        """
        metadata = cls(accession=accession)

        # Fetch project metadata
        project_url = f"{PRIDE_API_BASE}/projects/{accession}"
        try:
            response = requests.get(project_url, timeout=30)
            response.raise_for_status()
            project = response.json()
            metadata._raw_response = project
        except requests.RequestException as e:
            print(f"Warning: Failed to fetch PRIDE metadata: {e}")
            return metadata

        # Extract core fields
        metadata.dataset_id = MetadataField.auto(accession, "pride_api.accession")
        metadata.title = MetadataField.auto(
            project.get("title"), "pride_api.title"
        )

        # Lab ID from contacts
        contacts = project.get("contacts", [])
        if contacts:
            affiliation = contacts[0].get("affiliation")
            if affiliation:
                metadata.lab_id = MetadataField.auto(
                    affiliation, "pride_api.contacts[0].affiliation"
                )

        # Organism from organisms list
        organisms = project.get("organisms", [])
        if organisms:
            # Store all organisms with name and taxon_id
            all_orgs = []
            for org in organisms:
                org_entry = {"name": org.get("name", "")}
                accession_str = org.get("accession", "")
                if accession_str:
                    try:
                        org_entry["taxon_id"] = int(accession_str)
                    except ValueError:
                        pass
                all_orgs.append(org_entry)

            metadata.organisms_all = MetadataField.auto(
                all_orgs, "pride_api.organisms"
            )
            metadata.is_multi_organism = len(organisms) > 1

            # Keep backward-compatible single organism (first in list)
            org_name = organisms[0].get("name")
            if org_name:
                metadata.organism = MetadataField.auto(
                    org_name, "pride_api.organisms[0].name"
                )

        # Instrument from instruments list
        instruments = project.get("instruments", [])
        if instruments:
            instr_name = instruments[0].get("name")
            if instr_name:
                metadata.instrument = MetadataField.auto(
                    instr_name, "pride_api.instruments[0].name"
                )

        # DOI from references
        references = project.get("references", [])
        for ref in references:
            doi = ref.get("doi")
            if doi:
                metadata.publication_doi = MetadataField.auto(
                    doi, "pride_api.references[0].doi"
                )
                break

        # Parse protocol text for additional fields
        protocol = project.get("projectDescription", "") or ""
        protocol += "\n" + (project.get("sampleProcessingProtocol", "") or "")
        protocol += "\n" + (project.get("dataProcessingProtocol", "") or "")

        parsed = _parse_protocol_text(protocol)
        if parsed.get("gradient_length"):
            metadata.gradient_length = MetadataField.inferred(
                parsed["gradient_length"],
                "pride_api.projectDescription",
                confidence=parsed.get("gradient_confidence", 0.7),
            )
        if parsed.get("column_type"):
            metadata.column_type = MetadataField.inferred(
                parsed["column_type"],
                "pride_api.projectDescription",
                confidence=parsed.get("column_confidence", 0.7),
            )
        if parsed.get("acquisition_mode"):
            metadata.acquisition_mode = MetadataField.inferred(
                parsed["acquisition_mode"],
                "pride_api.projectDescription",
                confidence=parsed.get("mode_confidence", 0.7),
            )

        # Fetch file information
        files_url = f"{PRIDE_API_BASE}/files/byProject?accession={accession}"
        try:
            files_response = requests.get(files_url, timeout=30)
            files_response.raise_for_status()
            files_data = files_response.json()

            if files_data:
                metadata.num_files = MetadataField.auto(
                    len(files_data), "pride_api.files.count"
                )

                # Aggregate file types
                file_types = {}
                total_size = 0
                for f in files_data:
                    ext = Path(f.get("fileName", "")).suffix.lower()
                    file_types[ext] = file_types.get(ext, 0) + 1
                    total_size += f.get("fileSize", 0)

                metadata.file_types = MetadataField.auto(
                    file_types, "pride_api.files"
                )
                metadata.total_size_gb = MetadataField.auto(
                    round(total_size / (1024**3), 2), "pride_api.files"
                )

        except requests.RequestException:
            pass  # File info is optional

        return metadata

    @classmethod
    def from_yaml(cls, path: Path) -> "DatasetMetadata":
        """Load metadata from YAML file."""
        with open(path) as f:
            data = yaml.safe_load(f)

        metadata = cls(accession=data.get("accession", ""))
        metadata.version = data.get("version", "1.0")
        metadata.generated_at = data.get("generated_at", "")
        metadata.is_multi_organism = data.get("is_multi_organism", False)

        # Load title and DOI
        if data.get("title"):
            metadata.title = _dict_to_field(data["title"])
        if data.get("publication_doi"):
            metadata.publication_doi = _dict_to_field(data["publication_doi"])

        # Load fields
        fields = data.get("fields", {})
        for field_name, field_data in fields.items():
            if hasattr(metadata, field_name) and field_data:
                setattr(metadata, field_name, _dict_to_field(field_data))

        return metadata


def _dict_to_field(d: Dict[str, Any]) -> MetadataField:
    """Convert dictionary to MetadataField."""
    return MetadataField(
        value=d.get("value"),
        status=d.get("status", "missing"),
        source=d.get("source"),
        confidence=d.get("confidence", 1.0),
        hint=d.get("hint"),
    )


def _parse_protocol_text(protocol: str) -> Dict[str, Any]:
    """
    Parse protocol text to extract gradient_length, column_type, acquisition_mode.

    This uses regex patterns to find common descriptions of LC-MS parameters.
    """
    result = {}

    if not protocol:
        return result

    protocol_lower = protocol.lower()

    # Gradient length patterns
    gradient_patterns = [
        r"(\d+)\s*[-–]?\s*min(?:ute)?(?:s)?\s+gradient",
        r"gradient\s+(?:of\s+)?(\d+)\s*min",
        r"(\d+)\s*min\s+(?:linear\s+)?gradient",
        r"(\d+)\s*[-–]?\s*min(?:ute)?(?:s)?\s+(?:LC|liquid chromatography)",
        r"elution\s+(?:over|for)\s+(\d+)\s*min",
    ]
    for pattern in gradient_patterns:
        match = re.search(pattern, protocol_lower)
        if match:
            result["gradient_length"] = int(match.group(1))
            result["gradient_confidence"] = 0.8
            break

    # Column type patterns
    column_patterns = [
        r"(C18|C8|HILIC|RP|reverse[d]?\s*phase)[^\n,]*column",
        r"column[^\n,]*(C18|C8|HILIC|RP)",
        r"((?:\d+\s*[µu]m|\d+\s*cm)[^\n,]*(?:C18|C8|column))",
        r"(PepMap|Acclaim|BEH|HSS|CSH)[^\n,]*column",
    ]
    for pattern in column_patterns:
        match = re.search(pattern, protocol, re.IGNORECASE)
        if match:
            result["column_type"] = match.group(1).strip()
            result["column_confidence"] = 0.7
            break

    # Acquisition mode patterns
    mode_patterns = [
        (r"diaPASEF|dia-?PASEF", "diaPASEF", 0.95),
        (r"(?<!dia)PASEF(?!.*DIA)", "PASEF", 0.9),
        (r"data[- ]?independent|DIA", "DIA", 0.85),
        (r"data[- ]?dependent|DDA", "DDA", 0.85),
        (r"parallel reaction monitoring|PRM", "PRM", 0.9),
        (r"targeted|MRM|SRM", "targeted", 0.8),
    ]
    for pattern, mode, confidence in mode_patterns:
        if re.search(pattern, protocol, re.IGNORECASE):
            result["acquisition_mode"] = mode
            result["mode_confidence"] = confidence
            break

    return result


def fetch_pride_metadata(accession: str) -> DatasetMetadata:
    """
    Convenience function to fetch metadata from PRIDE.

    Args:
        accession: PRIDE accession (e.g., "PXD019086")

    Returns:
        DatasetMetadata with auto-extracted fields
    """
    return DatasetMetadata.from_pride_api(accession)


def save_raw_response(accession: str, output_dir: Path) -> Path:
    """
    Save raw PRIDE API response for debugging.

    Args:
        accession: PRIDE accession
        output_dir: Directory to save response

    Returns:
        Path to saved JSON file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    project_url = f"{PRIDE_API_BASE}/projects/{accession}"
    response = requests.get(project_url, timeout=30)
    response.raise_for_status()

    output_file = output_dir / "raw_api_response.json"
    with open(output_file, "w") as f:
        json.dump(response.json(), f, indent=2)

    return output_file


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Fetch metadata from PRIDE REST API"
    )
    parser.add_argument(
        "accession",
        help="PRIDE accession (e.g., PXD019086)"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        help="Output YAML file (default: stdout)"
    )
    parser.add_argument(
        "--save-raw",
        type=Path,
        help="Directory to save raw API response"
    )

    args = parser.parse_args()

    # Fetch metadata
    metadata = fetch_pride_metadata(args.accession)

    # Save raw response if requested
    if args.save_raw:
        raw_file = save_raw_response(args.accession, args.save_raw)
        print(f"Saved raw API response to {raw_file}")

    # Output
    if args.output:
        metadata.to_yaml(args.output)
        print(f"Wrote metadata to {args.output}")

        # Print validation summary
        summary = metadata.validation_summary()
        print(f"\nValidation summary:")
        print(f"  Complete: {summary['complete']}")
        print(f"  Auto fields: {summary['auto_fields']}")
        print(f"  Inferred fields: {summary['inferred_fields']}")
        print(f"  Missing fields: {summary['missing_fields']}")
        if summary["missing_field_names"]:
            print(f"  Missing: {', '.join(summary['missing_field_names'])}")
    else:
        # Print to stdout
        yaml_dict = metadata.to_yaml_dict()
        yaml.dump(yaml_dict, default_flow_style=False, sort_keys=False)
