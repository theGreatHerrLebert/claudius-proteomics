#!/usr/bin/env python3
"""
Extract metadata from timsTOF .d folders.

Sources:
- analysis.tdf/GlobalMetadata: instrument, operator, method, datetime
- SampleInfo.xml: LC/MS method names (parse gradient length, mode)
- Frames table: actual run duration
- chromatography-data.sqlite: LC system type
"""

import re
import sqlite3
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional


# Patterns for extracting gradient length from LC method names
GRADIENT_PATTERNS = [
    r"(\d+)min",           # "120min" → 120
    r"(\d+)_min",          # "120_min" → 120
    r"_(\d+)m_",           # "_60m_" → 60
    r"(\d+)m[in]*_",       # "60min_" or "60m_" → 60
    r"gradient[_-]?(\d+)", # "gradient_60" → 60
    r"(\d+)[_-]?gradient", # "60_gradient" → 60
]

# Patterns for extracting acquisition mode from MS method names
MODE_PATTERNS = [
    (r"diaPASEF", "diaPASEF"),
    (r"dia[_-]?PASEF", "diaPASEF"),
    (r"PASEF", "PASEF"),
    (r"pasef", "PASEF"),
    (r"DDA", "DDA"),
    (r"dda", "DDA"),
    (r"DIA", "DIA"),
    (r"dia", "DIA"),
    (r"PRM", "PRM"),
    (r"prm", "PRM"),
    (r"MRM", "MRM"),
    (r"SRM", "SRM"),
]


@dataclass
class RawMetadata:
    """Metadata extracted from a single .d folder."""

    # Source file
    d_path: str
    d_name: str

    # From GlobalMetadata
    instrument_name: Optional[str] = None
    instrument_serial: Optional[str] = None
    operator_name: Optional[str] = None
    acquisition_datetime: Optional[str] = None
    ms_method_name: Optional[str] = None

    # From SampleInfo.xml
    lc_method_name: Optional[str] = None
    ms_method_name_hystar: Optional[str] = None
    sample_name: Optional[str] = None

    # Derived from method names
    gradient_length_minutes: Optional[int] = None
    acquisition_mode: Optional[str] = None

    # From Frames table
    run_duration_seconds: Optional[float] = None
    run_duration_minutes: Optional[float] = None
    num_frames: Optional[int] = None
    num_ms1_frames: Optional[int] = None
    num_ms2_frames: Optional[int] = None

    # From chromatography-data.sqlite
    lc_system: Optional[str] = None

    # Extraction metadata
    extraction_timestamp: str = field(default_factory=lambda: datetime.now().isoformat())
    extraction_errors: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "d_path": self.d_path,
            "d_name": self.d_name,
            "instrument_name": self.instrument_name,
            "instrument_serial": self.instrument_serial,
            "operator_name": self.operator_name,
            "acquisition_datetime": self.acquisition_datetime,
            "ms_method_name": self.ms_method_name,
            "lc_method_name": self.lc_method_name,
            "ms_method_name_hystar": self.ms_method_name_hystar,
            "sample_name": self.sample_name,
            "gradient_length_minutes": self.gradient_length_minutes,
            "acquisition_mode": self.acquisition_mode,
            "run_duration_seconds": self.run_duration_seconds,
            "run_duration_minutes": self.run_duration_minutes,
            "num_frames": self.num_frames,
            "num_ms1_frames": self.num_ms1_frames,
            "num_ms2_frames": self.num_ms2_frames,
            "lc_system": self.lc_system,
            "extraction_timestamp": self.extraction_timestamp,
            "extraction_errors": self.extraction_errors,
        }


def parse_gradient_from_method_name(method_name: str) -> Optional[int]:
    """
    Extract gradient length from LC method name.

    Examples:
    - "20181107_120min_UpstairsCopy.m" → 120
    - "60min_gradient.m" → 60
    - "gradient_90.m" → 90
    """
    if not method_name:
        return None

    for pattern in GRADIENT_PATTERNS:
        match = re.search(pattern, method_name, re.IGNORECASE)
        if match:
            try:
                return int(match.group(1))
            except (ValueError, IndexError):
                continue

    return None


def parse_acquisition_mode(ms_method_name: str) -> Optional[str]:
    """
    Extract acquisition mode from MS method name.

    Examples:
    - "PASEF100ms_linearCE.m" → "PASEF"
    - "DDA_standard.m" → "DDA"
    - "diaPASEF.m" → "diaPASEF"
    """
    if not ms_method_name:
        return None

    for pattern, mode in MODE_PATTERNS:
        if re.search(pattern, ms_method_name):
            return mode

    return None


def _extract_global_metadata(tdf_path: Path) -> Dict[str, Any]:
    """Extract metadata from analysis.tdf GlobalMetadata table."""
    result = {}

    try:
        conn = sqlite3.connect(str(tdf_path))
        cursor = conn.cursor()

        # GlobalMetadata is a key-value table
        cursor.execute("SELECT Key, Value FROM GlobalMetadata")
        rows = cursor.fetchall()

        metadata = {row[0]: row[1] for row in rows}

        result["instrument_name"] = metadata.get("InstrumentName")
        result["instrument_serial"] = metadata.get("InstrumentSerialNumber")
        result["operator_name"] = metadata.get("OperatorName")
        result["acquisition_datetime"] = metadata.get("AcquisitionDateTime")
        result["ms_method_name"] = metadata.get("MethodName")

        conn.close()
    except Exception as e:
        result["_error"] = str(e)

    return result


def _extract_frames_info(tdf_path: Path) -> Dict[str, Any]:
    """Extract frame statistics and run duration from Frames table."""
    result = {}

    try:
        conn = sqlite3.connect(str(tdf_path))
        cursor = conn.cursor()

        # Get frame count and types
        cursor.execute("SELECT COUNT(*) FROM Frames")
        result["num_frames"] = cursor.fetchone()[0]

        # MsMsType: 0 = MS1, 8/9 = PASEF MS2
        cursor.execute("SELECT COUNT(*) FROM Frames WHERE MsMsType = 0")
        result["num_ms1_frames"] = cursor.fetchone()[0]

        cursor.execute("SELECT COUNT(*) FROM Frames WHERE MsMsType != 0")
        result["num_ms2_frames"] = cursor.fetchone()[0]

        # Get run duration from Time column (seconds)
        cursor.execute("SELECT MIN(Time), MAX(Time) FROM Frames")
        time_range = cursor.fetchone()
        if time_range[0] is not None and time_range[1] is not None:
            duration_sec = time_range[1] - time_range[0]
            result["run_duration_seconds"] = round(duration_sec, 2)
            result["run_duration_minutes"] = round(duration_sec / 60, 2)

        conn.close()
    except Exception as e:
        result["_error"] = str(e)

    return result


def _extract_sample_info(sample_info_path: Path) -> Dict[str, Any]:
    """Extract metadata from SampleInfo.xml."""
    result = {}

    try:
        tree = ET.parse(str(sample_info_path))
        root = tree.getroot()

        # Look for relevant elements - structure varies by HyStar version
        # Common elements: HyStar_LC_Method_Name, HyStar_MS_Method_Name, SampleName

        def find_text(tag: str) -> Optional[str]:
            elem = root.find(f".//{tag}")
            if elem is not None and elem.text:
                return elem.text.strip()
            return None

        result["lc_method_name"] = find_text("HyStar_LC_Method_Name")
        result["ms_method_name_hystar"] = find_text("HyStar_MS_Method_Name")
        result["sample_name"] = find_text("SampleName")

        # Alternative element names
        if not result["lc_method_name"]:
            result["lc_method_name"] = find_text("LC_Method_Name")
        if not result["ms_method_name_hystar"]:
            result["ms_method_name_hystar"] = find_text("MS_Method_Name")

    except Exception as e:
        result["_error"] = str(e)

    return result


def _extract_chromatography_info(chrom_db_path: Path) -> Dict[str, Any]:
    """Extract LC system info from chromatography-data.sqlite."""
    result = {}

    try:
        conn = sqlite3.connect(str(chrom_db_path))
        cursor = conn.cursor()

        # TraceSources table contains instrument info
        cursor.execute("""
            SELECT DISTINCT Instrument FROM TraceSources
            WHERE Instrument IS NOT NULL AND Instrument != ''
        """)
        instruments = [row[0] for row in cursor.fetchall()]

        if instruments:
            result["lc_system"] = instruments[0]  # Take first if multiple

        conn.close()
    except Exception as e:
        result["_error"] = str(e)

    return result


def extract_raw_metadata(d_path: Path) -> RawMetadata:
    """
    Extract metadata from a .d folder.

    Args:
        d_path: Path to the .d folder

    Returns:
        RawMetadata with extracted fields
    """
    d_path = Path(d_path)

    metadata = RawMetadata(
        d_path=str(d_path.absolute()),
        d_name=d_path.name,
    )

    errors = []

    # Extract from analysis.tdf
    tdf_path = d_path / "analysis.tdf"
    if tdf_path.exists():
        global_meta = _extract_global_metadata(tdf_path)
        if "_error" in global_meta:
            errors.append(f"GlobalMetadata: {global_meta['_error']}")
        else:
            metadata.instrument_name = global_meta.get("instrument_name")
            metadata.instrument_serial = global_meta.get("instrument_serial")
            metadata.operator_name = global_meta.get("operator_name")
            metadata.acquisition_datetime = global_meta.get("acquisition_datetime")
            metadata.ms_method_name = global_meta.get("ms_method_name")

        frames_info = _extract_frames_info(tdf_path)
        if "_error" in frames_info:
            errors.append(f"Frames: {frames_info['_error']}")
        else:
            metadata.run_duration_seconds = frames_info.get("run_duration_seconds")
            metadata.run_duration_minutes = frames_info.get("run_duration_minutes")
            metadata.num_frames = frames_info.get("num_frames")
            metadata.num_ms1_frames = frames_info.get("num_ms1_frames")
            metadata.num_ms2_frames = frames_info.get("num_ms2_frames")
    else:
        errors.append("analysis.tdf not found")

    # Extract from SampleInfo.xml
    sample_info_path = d_path / "SampleInfo.xml"
    if sample_info_path.exists():
        sample_info = _extract_sample_info(sample_info_path)
        if "_error" in sample_info:
            errors.append(f"SampleInfo.xml: {sample_info['_error']}")
        else:
            metadata.lc_method_name = sample_info.get("lc_method_name")
            metadata.ms_method_name_hystar = sample_info.get("ms_method_name_hystar")
            metadata.sample_name = sample_info.get("sample_name")
    else:
        errors.append("SampleInfo.xml not found")

    # Extract from chromatography-data.sqlite
    chrom_db_path = d_path / "chromatography-data.sqlite"
    if chrom_db_path.exists():
        chrom_info = _extract_chromatography_info(chrom_db_path)
        if "_error" in chrom_info:
            errors.append(f"chromatography-data.sqlite: {chrom_info['_error']}")
        else:
            metadata.lc_system = chrom_info.get("lc_system")
    # Note: chromatography-data.sqlite may not exist in all .d folders

    # Derive gradient length from LC method name
    lc_method = metadata.lc_method_name or metadata.ms_method_name
    if lc_method:
        metadata.gradient_length_minutes = parse_gradient_from_method_name(lc_method)

    # Fallback: use run duration if gradient not parsed
    if metadata.gradient_length_minutes is None and metadata.run_duration_minutes:
        # Round to nearest common gradient length
        duration = metadata.run_duration_minutes
        common_gradients = [15, 30, 45, 60, 90, 120, 180, 240]
        for cg in common_gradients:
            if abs(duration - cg) < 5:  # Within 5 minutes
                metadata.gradient_length_minutes = cg
                break

    # Derive acquisition mode from MS method name
    ms_method = metadata.ms_method_name_hystar or metadata.ms_method_name
    if ms_method:
        metadata.acquisition_mode = parse_acquisition_mode(ms_method)

    metadata.extraction_errors = errors

    return metadata


def get_run_duration_minutes(d_path: Path) -> Optional[float]:
    """
    Get actual run duration from Frames table.

    This is a convenience function for quick duration lookup.
    """
    tdf_path = d_path / "analysis.tdf"
    if not tdf_path.exists():
        return None

    frames_info = _extract_frames_info(tdf_path)
    return frames_info.get("run_duration_minutes")


def extract_all_raw_metadata(raw_dir: Path) -> List[RawMetadata]:
    """
    Extract metadata from all .d folders in a directory.

    Args:
        raw_dir: Directory containing .d folders

    Returns:
        List of RawMetadata for each .d folder
    """
    raw_dir = Path(raw_dir)
    d_folders = sorted(raw_dir.glob("*.d"))

    results = []
    for d_folder in d_folders:
        if d_folder.is_dir():
            results.append(extract_raw_metadata(d_folder))

    return results


def aggregate_metadata(metadata_list: List[RawMetadata]) -> Dict[str, Any]:
    """
    Aggregate metadata from multiple .d folders into a single summary.

    Checks consistency across files and returns consensus values.
    """
    if not metadata_list:
        return {}

    # Collect values for consistency checking
    instruments = [m.instrument_name for m in metadata_list if m.instrument_name]
    lc_systems = [m.lc_system for m in metadata_list if m.lc_system]
    gradients = [m.gradient_length_minutes for m in metadata_list if m.gradient_length_minutes]
    modes = [m.acquisition_mode for m in metadata_list if m.acquisition_mode]
    durations = [m.run_duration_minutes for m in metadata_list if m.run_duration_minutes]

    def most_common(values: List[Any]) -> Optional[Any]:
        if not values:
            return None
        from collections import Counter
        counter = Counter(values)
        return counter.most_common(1)[0][0]

    def check_consistency(values: List[Any], name: str) -> Dict[str, Any]:
        if not values:
            return {"value": None, "consistent": True, "unique_values": []}
        unique = list(set(values))
        return {
            "value": most_common(values),
            "consistent": len(unique) == 1,
            "unique_values": unique if len(unique) > 1 else [],
        }

    return {
        "num_files": len(metadata_list),
        "instrument": check_consistency(instruments, "instrument"),
        "lc_system": check_consistency(lc_systems, "lc_system"),
        "gradient_length": check_consistency(gradients, "gradient_length"),
        "acquisition_mode": check_consistency(modes, "acquisition_mode"),
        "run_duration": {
            "min": min(durations) if durations else None,
            "max": max(durations) if durations else None,
            "mean": sum(durations) / len(durations) if durations else None,
        },
        "files": [m.d_name for m in metadata_list],
    }


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(
        description="Extract metadata from timsTOF .d folders"
    )
    parser.add_argument(
        "d_path",
        type=Path,
        help="Path to a .d folder or directory containing .d folders"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        help="Output JSON file (default: stdout)"
    )
    parser.add_argument(
        "--aggregate",
        action="store_true",
        help="Aggregate metadata from multiple .d folders"
    )

    args = parser.parse_args()

    if args.d_path.suffix == ".d" and args.d_path.is_dir():
        # Single .d folder
        metadata = extract_raw_metadata(args.d_path)
        result = metadata.to_dict()
    elif args.d_path.is_dir():
        # Directory with multiple .d folders
        metadata_list = extract_all_raw_metadata(args.d_path)
        if args.aggregate:
            result = aggregate_metadata(metadata_list)
        else:
            result = [m.to_dict() for m in metadata_list]
    else:
        print(f"Error: {args.d_path} is not a valid .d folder or directory")
        exit(1)

    output_json = json.dumps(result, indent=2)

    if args.output:
        args.output.write_text(output_json)
        print(f"Wrote metadata to {args.output}")
    else:
        print(output_json)
