#!/usr/bin/env python3
"""
Paper-based metadata extraction for San José.

Downloads open-access PDFs via Unpaywall, extracts text, and uses Claude API
to pull structured LC-MS parameters from Materials & Methods sections.

Graceful degradation: every failure is a skip, never a block.
"""

import json
import os
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml


def download_paper_pdf(
    doi: str,
    papers_dir: Path,
    unpaywall_email: Optional[str] = None,
) -> Optional[Path]:
    """
    Download open-access PDF via Unpaywall API.

    Args:
        doi: Publication DOI (e.g., "10.1038/s41467-021-21352-8")
        papers_dir: Directory to save PDFs
        unpaywall_email: Email for Unpaywall API (required by their TOS)

    Returns:
        Path to downloaded PDF, or None if unavailable/paywalled
    """
    papers_dir = Path(papers_dir)
    papers_dir.mkdir(parents=True, exist_ok=True)

    # Sanitize DOI for filename
    sanitized = doi.replace("/", "_").replace("\\", "_")
    pdf_path = papers_dir / f"{sanitized}.pdf"

    # Skip if already downloaded
    if pdf_path.exists() and pdf_path.stat().st_size > 0:
        print(f"    PDF already exists: {pdf_path.name}")
        return pdf_path

    if not unpaywall_email:
        print("    No unpaywall_email configured, skipping PDF download")
        print("    Hint: Set paper_extraction.unpaywall_email in config.yaml")
        return None

    try:
        import requests
    except ImportError:
        print("    requests not installed, skipping PDF download")
        return None

    # Query Unpaywall API
    api_url = f"https://api.unpaywall.org/v2/{doi}?email={unpaywall_email}"
    try:
        resp = requests.get(api_url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except Exception as e:
        print(f"    Unpaywall API error: {e}")
        return None

    # Find best OA PDF URL
    pdf_url = None
    best_oa = data.get("best_oa_location")
    if best_oa:
        pdf_url = best_oa.get("url_for_pdf") or best_oa.get("url")

    if not pdf_url:
        # Check all OA locations
        for loc in data.get("oa_locations", []):
            url = loc.get("url_for_pdf") or loc.get("url")
            if url and url.endswith(".pdf"):
                pdf_url = url
                break

    if not pdf_url:
        print(f"    Paper is paywalled or no OA PDF found for DOI: {doi}")
        print(f"    You can manually place the PDF in: {papers_dir}/")
        return None

    # Download PDF
    try:
        print(f"    Downloading PDF from: {pdf_url[:80]}...")
        pdf_resp = requests.get(pdf_url, timeout=120, stream=True)
        pdf_resp.raise_for_status()

        # Write to file
        with open(pdf_path, "wb") as f:
            for chunk in pdf_resp.iter_content(chunk_size=8192):
                f.write(chunk)

        # Validate PDF magic bytes
        with open(pdf_path, "rb") as f:
            magic = f.read(5)
        if magic != b"%PDF-":
            print(f"    Downloaded file is not a valid PDF (got {magic!r}), removing")
            pdf_path.unlink()
            return None

        size_kb = pdf_path.stat().st_size / 1024
        print(f"    Downloaded: {pdf_path.name} ({size_kb:.0f} KB)")
        return pdf_path

    except Exception as e:
        print(f"    PDF download failed: {e}")
        if pdf_path.exists():
            pdf_path.unlink()
        return None


def extract_pdf_text(pdf_path: Path, max_chars: int = 80_000) -> Optional[str]:
    """
    Extract text from PDF, prioritizing methods-related pages.

    Tries pdfplumber first, falls back to PyPDF2.

    Args:
        pdf_path: Path to PDF file
        max_chars: Maximum characters to return

    Returns:
        Extracted text, or None on failure
    """
    pdf_path = Path(pdf_path)
    if not pdf_path.exists():
        return None

    text = None

    # Try pdfplumber first
    try:
        import pdfplumber

        pages_text = []
        with pdfplumber.open(pdf_path) as pdf:
            for page in pdf.pages:
                page_text = page.extract_text()
                if page_text:
                    pages_text.append(page_text)
        if pages_text:
            text = "\n\n".join(pages_text)
    except ImportError:
        pass
    except Exception as e:
        print(f"    pdfplumber failed: {e}, trying PyPDF2...")

    # Fallback to PyPDF2
    if text is None:
        try:
            from PyPDF2 import PdfReader

            reader = PdfReader(str(pdf_path))
            pages_text = []
            for page in reader.pages:
                page_text = page.extract_text()
                if page_text:
                    pages_text.append(page_text)
            if pages_text:
                text = "\n\n".join(pages_text)
        except ImportError:
            print("    Neither pdfplumber nor PyPDF2 installed, cannot extract text")
            print("    Install with: pip install pdfplumber")
            return None
        except Exception as e:
            print(f"    PyPDF2 failed: {e}")
            return None

    if not text:
        print(f"    No text extracted from {pdf_path.name}")
        return None

    # If text exceeds limit, prioritize methods-related pages
    if len(text) > max_chars:
        text = _prioritize_methods_text(text, max_chars)

    return text


def _prioritize_methods_text(full_text: str, max_chars: int) -> str:
    """
    If text is too long, prioritize pages with methods-related keywords.
    Split by double-newline (page breaks), score each section, take highest-scoring.
    """
    methods_keywords = [
        "methods", "materials", "gradient", "column", "chromatograph",
        "lc-ms", "lc/ms", "nano-lc", "nanolc", "nanoelute", "timstof",
        "mobile phase", "flow rate", "acetonitrile", "formic acid",
        "sample prep", "digestion", "trypsin", "pasef", "acquisition",
    ]

    sections = full_text.split("\n\n")

    # Score each section by keyword density
    scored = []
    for section in sections:
        lower = section.lower()
        score = sum(1 for kw in methods_keywords if kw in lower)
        scored.append((score, section))

    # Sort by score descending, take highest-scoring sections up to limit
    scored.sort(key=lambda x: x[0], reverse=True)

    result_parts = []
    total_chars = 0
    for score, section in scored:
        if total_chars + len(section) > max_chars:
            # Add as much as we can
            remaining = max_chars - total_chars
            if remaining > 200:
                result_parts.append(section[:remaining])
            break
        result_parts.append(section)
        total_chars += len(section)

    return "\n\n".join(result_parts)


def extract_metadata_with_llm(
    paper_text: str,
    accession: str,
    api_key: str,
    model: str = "claude-sonnet-4-20250514",
) -> Optional[Dict[str, Any]]:
    """
    Use Claude API to extract structured LC-MS metadata from paper text.

    Args:
        paper_text: Extracted text from the publication PDF
        accession: PRIDE accession for context
        api_key: Anthropic API key
        model: Claude model to use

    Returns:
        Dict with extracted fields, or None on failure
    """
    try:
        import anthropic
    except ImportError:
        print("    anthropic package not installed, skipping LLM extraction")
        print("    Install with: pip install anthropic")
        return None

    prompt = f"""You are extracting LC-MS/MS experimental parameters from a proteomics publication.
The data is deposited in PRIDE as {accession}.

Extract the following fields from the paper text. For each field, provide:
- "value": the extracted value
- "confidence": float 0.0-1.0 indicating how certain you are
- "evidence": the exact quote from the paper supporting this value

Fields to extract:
1. "gradient_length" - LC gradient length in minutes (integer). Look for statements like "X min gradient" or gradient time descriptions.
2. "column_type" - Full column specification (e.g., "25 cm x 75 um IonOpticks Aurora C18 1.6 um"). Include length, ID, packing material, particle size.
3. "lc_system" - LC system name (e.g., "nanoElute 2", "EASY-nLC 1200", "Vanquish Neo")
4. "acquisition_mode" - One of: DDA, DIA, PASEF, diaPASEF, PRM, targeted
5. "sample_prep" - Brief sample preparation description (e.g., "SP3 bead cleanup, trypsin digestion")
6. "lab_id" - Corresponding author name and institution
7. "additional_lc_params" - Object with any of: flow_rate, mobile_phase_a, mobile_phase_b, column_temp, injection_volume

Return ONLY valid JSON (no markdown code fences) with this structure:
{{
  "gradient_length": {{"value": 120, "confidence": 0.95, "evidence": "...exact quote..."}},
  "column_type": {{"value": "...", "confidence": 0.9, "evidence": "..."}},
  "lc_system": {{"value": "...", "confidence": 0.9, "evidence": "..."}},
  "acquisition_mode": {{"value": "PASEF", "confidence": 0.95, "evidence": "..."}},
  "sample_prep": {{"value": "...", "confidence": 0.8, "evidence": "..."}},
  "lab_id": {{"value": "...", "confidence": 0.9, "evidence": "..."}},
  "additional_lc_params": {{"value": {{...}}, "confidence": 0.7, "evidence": "..."}}
}}

If a field cannot be determined from the text, set value to null and confidence to 0.0.

Paper text:
{paper_text}"""

    try:
        client = anthropic.Anthropic(api_key=api_key)
        response = client.messages.create(
            model=model,
            max_tokens=2000,
            messages=[{"role": "user", "content": prompt}],
        )

        response_text = response.content[0].text

        # Try to parse JSON directly
        try:
            result = json.loads(response_text)
        except json.JSONDecodeError:
            # Try to extract JSON block from response
            result = _extract_json_from_text(response_text)

        if result is None:
            print("    LLM returned invalid JSON, skipping")
            return None

        # Attach token usage for logging
        result["_usage"] = {
            "model": model,
            "input_tokens": response.usage.input_tokens,
            "output_tokens": response.usage.output_tokens,
        }

        return result

    except Exception as e:
        print(f"    LLM extraction failed: {e}")
        return None


def _extract_json_from_text(text: str) -> Optional[Dict]:
    """Try to extract a JSON object from text that may contain markdown fences."""
    # Try code fence extraction
    patterns = [
        r"```json\s*\n(.*?)\n```",
        r"```\s*\n(.*?)\n```",
        r"\{[\s\S]*\}",
    ]
    for pattern in patterns:
        match = re.search(pattern, text, re.DOTALL)
        if match:
            try:
                candidate = match.group(1) if match.lastindex else match.group(0)
                return json.loads(candidate)
            except (json.JSONDecodeError, IndexError):
                continue
    return None


def run_paper_extraction(
    accession: str,
    metadata_dir: Path,
    config: Dict[str, Any],
) -> Optional[Dict[str, Any]]:
    """
    Orchestrate paper-based metadata extraction.

    1. Load DOI from pride_metadata.yaml
    2. Try Unpaywall download, check for manual PDFs
    3. Extract text from PDFs
    4. Call LLM extraction if API key available
    5. Save paper_extraction.yaml
    6. Merge into pride_metadata.yaml

    Args:
        accession: PRIDE accession
        metadata_dir: Path to metadata directory for this dataset
        config: Pipeline configuration dict

    Returns:
        Dict with extraction results, or None if skipped entirely
    """
    metadata_dir = Path(metadata_dir)
    paper_config = config.get("paper_extraction", {})

    # Check if disabled
    if not paper_config.get("enabled", True):
        print("    Paper extraction disabled in config")
        return None

    # Check skip_if_exists
    extraction_path = metadata_dir / "paper_extraction.yaml"
    if paper_config.get("skip_if_exists", True) and extraction_path.exists():
        print(f"    Paper extraction already exists: {extraction_path}")
        return _load_yaml(extraction_path)

    # Load DOI from pride_metadata.yaml
    pride_metadata_path = metadata_dir / "pride_metadata.yaml"
    doi = None
    if pride_metadata_path.exists():
        pride_data = _load_yaml(pride_metadata_path)
        if pride_data:
            doi_field = pride_data.get("publication_doi")
            if isinstance(doi_field, dict):
                doi = doi_field.get("value")
            elif isinstance(doi_field, str):
                doi = doi_field

    papers_dir = metadata_dir / "papers"
    papers_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Try to download PDF
    pdf_sources = []
    if doi:
        unpaywall_email = paper_config.get("unpaywall_email")
        pdf_path = download_paper_pdf(doi, papers_dir, unpaywall_email)
        if pdf_path:
            pdf_sources.append({
                "filename": pdf_path.name,
                "source": "unpaywall",
            })
    else:
        print("    No DOI found in PRIDE metadata, skipping download")
        print(f"    You can manually place PDFs in: {papers_dir}/")

    # Step 2: Find all PDFs (auto-downloaded + manually placed)
    all_pdfs = sorted(papers_dir.glob("*.pdf"))
    if not all_pdfs:
        print(f"    No PDFs found in {papers_dir}/")
        # Save minimal result
        result = {
            "accession": accession,
            "extracted_at": datetime.now().isoformat(),
            "pdf_sources": [],
            "fields": {},
            "status": "no_pdf",
        }
        _save_yaml(result, extraction_path)
        return result

    # Track manually-placed PDFs
    auto_filenames = {s["filename"] for s in pdf_sources}
    for pdf in all_pdfs:
        if pdf.name not in auto_filenames:
            pdf_sources.append({
                "filename": pdf.name,
                "source": "manual",
            })

    # Step 3: Extract text from all PDFs
    max_chars = paper_config.get("max_pdf_chars", 80_000)
    combined_text = []
    for pdf in all_pdfs:
        print(f"    Extracting text from: {pdf.name}")
        text = extract_pdf_text(pdf, max_chars=max_chars)
        if text:
            # Track chars per source
            for src in pdf_sources:
                if src["filename"] == pdf.name:
                    src["chars_extracted"] = len(text)
            combined_text.append(text)

    if not combined_text:
        print("    No text extracted from any PDF")
        result = {
            "accession": accession,
            "extracted_at": datetime.now().isoformat(),
            "pdf_sources": pdf_sources,
            "fields": {},
            "status": "no_text",
        }
        _save_yaml(result, extraction_path)
        return result

    full_text = "\n\n---\n\n".join(combined_text)

    # Step 4: LLM extraction
    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        print("    ANTHROPIC_API_KEY not set, skipping LLM extraction")
        print("    Set the env var and re-run with --enrich-metadata to extract")
        result = {
            "accession": accession,
            "extracted_at": datetime.now().isoformat(),
            "pdf_sources": pdf_sources,
            "fields": {},
            "status": "no_api_key",
        }
        _save_yaml(result, extraction_path)
        return result

    model = paper_config.get("model", "claude-sonnet-4-20250514")
    print(f"    Calling {model} for metadata extraction...")
    llm_result = extract_metadata_with_llm(full_text, accession, api_key, model)

    if not llm_result:
        result = {
            "accession": accession,
            "extracted_at": datetime.now().isoformat(),
            "model": model,
            "pdf_sources": pdf_sources,
            "fields": {},
            "status": "llm_failed",
        }
        _save_yaml(result, extraction_path)
        return result

    # Step 5: Build structured result
    usage = llm_result.pop("_usage", {})
    fields = {}
    for field_name, field_data in llm_result.items():
        if not isinstance(field_data, dict):
            continue
        value = field_data.get("value")
        if value is None:
            continue
        # Find which PDF this came from (use first with text)
        source_pdf = next(
            (s["filename"] for s in pdf_sources if s.get("chars_extracted", 0) > 0),
            "unknown",
        )
        fields[field_name] = {
            "value": value,
            "status": "inferred",
            "source": f"paper:{source_pdf}",
            "confidence": field_data.get("confidence", 0.8),
            "evidence": field_data.get("evidence", ""),
        }

    n_fields = len(fields)
    print(f"    Extracted {n_fields} fields from paper")

    result = {
        "accession": accession,
        "extracted_at": datetime.now().isoformat(),
        "model": model,
        "pdf_sources": pdf_sources,
        "fields": fields,
        "raw_llm_response": usage,
        "status": "success",
    }

    # Step 6: Save and merge
    _save_yaml(result, extraction_path)
    print(f"    Saved: {extraction_path}")

    # Merge into pride_metadata.yaml
    if pride_metadata_path.exists():
        try:
            try:
                from scripts.pride_metadata import DatasetMetadata
            except ImportError:
                import sys
                sys.path.insert(0, str(Path(__file__).parent.parent))
                from scripts.pride_metadata import DatasetMetadata

            metadata = DatasetMetadata.from_yaml(pride_metadata_path)
            metadata.merge_paper_extraction(extraction_path)
            metadata.to_yaml(pride_metadata_path)
            print(f"    Merged paper fields into {pride_metadata_path.name}")
        except Exception as e:
            print(f"    Could not merge into pride_metadata.yaml: {e}")

    return result


def _load_yaml(path: Path) -> Optional[Dict]:
    """Load YAML file, return None on error."""
    try:
        with open(path) as f:
            return yaml.safe_load(f)
    except Exception:
        return None


def _save_yaml(data: Dict, path: Path) -> None:
    """Save dict to YAML file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Extract LC-MS metadata from publication PDFs"
    )
    parser.add_argument(
        "accession",
        help="PRIDE accession (e.g., PXD019086)",
    )
    parser.add_argument(
        "--metadata-dir",
        type=Path,
        help="Metadata directory (default: data/metadata/{accession})",
    )
    parser.add_argument(
        "--config", "-c",
        default="config/config.yaml",
        help="Config file (default: config/config.yaml)",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Skip PDF download (use manually-placed PDFs only)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run even if paper_extraction.yaml exists",
    )

    args = parser.parse_args()

    # Determine metadata dir
    metadata_dir = args.metadata_dir or Path(f"data/metadata/{args.accession}")

    # Load config
    config_path = Path(args.config)
    if config_path.exists():
        with open(config_path) as f:
            config = yaml.safe_load(f) or {}
    else:
        config = {}

    # Apply CLI overrides
    if args.skip_download:
        config.setdefault("paper_extraction", {})["unpaywall_email"] = None
    if args.force:
        config.setdefault("paper_extraction", {})["skip_if_exists"] = False

    # Run
    print(f"Paper metadata extraction for {args.accession}")
    print(f"  Metadata dir: {metadata_dir}")
    result = run_paper_extraction(args.accession, metadata_dir, config)

    if result:
        status = result.get("status", "unknown")
        n_fields = len(result.get("fields", {}))
        print(f"\nResult: status={status}, fields={n_fields}")
        if result.get("fields"):
            for name, field in result["fields"].items():
                val = field.get("value", "?")
                conf = field.get("confidence", "?")
                print(f"  {name}: {val} (confidence={conf})")
    else:
        print("\nResult: skipped")
