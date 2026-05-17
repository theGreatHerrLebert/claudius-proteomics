#!/usr/bin/env python3
"""
Download raw data from PRIDE Archive with metadata extraction.

This script:
1. Fetches metadata from PRIDE REST API
2. Parses protocol text for gradient/column/mode hints
3. Downloads raw files over HTTPS from the ftp.ebi.ac.uk mirror
4. Creates metadata YAML with validation status

Usage:
    # Metadata only (no download)
    python download_pride.py --accession PXD019086 --metadata-only

    # Full download
    python download_pride.py --accession PXD019086 --output data/raw/PXD019086
"""

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Optional

import requests

from pride_metadata import (
    DatasetMetadata,
    MetadataField,
    _parse_protocol_text,
    fetch_pride_metadata,
    save_raw_response,
)

# PRIDE REST API v2 base URL
PRIDE_API_BASE = "https://www.ebi.ac.uk/pride/ws/archive/v2"

# ProteomeXchange API for cross-repository resolution
PX_API_BASE = "https://proteomecentral.proteomexchange.org/cgi/GetDataset"

# PRIDE raw files are served from ftp.pride.ebi.ac.uk over FTP/Aspera. Those
# protocols (and that host) are unreachable from networks that only permit
# HTTPS through an allow-list proxy — e.g. the Mogon NHR HPC cluster. The
# general EBI server ftp.ebi.ac.uk mirrors the identical /pride/data/archive
# tree over HTTPS, so every PRIDE file URL is rewritten to that host before
# download.
EBI_HTTPS_BASE = "https://ftp.ebi.ac.uk"


def to_ebi_https_url(url: str) -> str:
    """Rewrite a PRIDE file URL to the proxy-reachable ftp.ebi.ac.uk HTTPS mirror.

    ``ftp://ftp.pride.ebi.ac.uk/pride/data/archive/...`` (and the http[s]
    variants) become ``https://ftp.ebi.ac.uk/pride/data/archive/...``. The
    archive path is identical between the two servers. Any other URL — and
    Aspera ``prd_ascp@...`` locations, which cannot be rewritten — is returned
    unchanged.
    """
    for host_prefix in (
        "ftp://ftp.pride.ebi.ac.uk",
        "https://ftp.pride.ebi.ac.uk",
        "http://ftp.pride.ebi.ac.uk",
    ):
        if url.startswith(host_prefix):
            return EBI_HTTPS_BASE + url[len(host_prefix):]
    return url


def _resolve_repository(accession: str) -> tuple:
    """
    Resolve which repository hosts a PXD accession.

    Tries PRIDE first, then falls back to ProteomeXchange to discover
    the actual repository (jPOST, MassIVE, etc.).

    Args:
        accession: ProteomeXchange accession (e.g., PXD019746)

    Returns:
        Tuple of (repo_name, repo_info) where repo_name is "pride", "jpost",
        "massive", etc. and repo_info contains FTP URLs and metadata.
    """
    # Try PRIDE first
    proj_url = f"{PRIDE_API_BASE}/projects/{accession}"
    try:
        resp = requests.get(proj_url, timeout=30)
        if resp.status_code == 200:
            return ("pride", resp.json())
    except requests.RequestException:
        pass

    # PRIDE returned 404 or failed — query ProteomeXchange
    print(f"   PRIDE API returned 404 for {accession}, querying ProteomeXchange...")
    px_url = f"{PX_API_BASE}?ID={accession}&outputMode=json"
    try:
        resp = requests.get(px_url, timeout=60)
        resp.raise_for_status()

        # ProteomeXchange may return JSONP; strip callback wrapper
        text = resp.text.strip()
        if text.startswith("merged(") and text.endswith(");"):
            text = text[len("merged("):-len(");")]
        elif text.startswith("(") and text.endswith(")"):
            text = text[1:-1]
        px_data = json.loads(text)
    except Exception as e:
        print(f"   ProteomeXchange query failed: {e}")
        return ("unknown", {})

    # Find FTP links in fullDatasetLinks
    ftp_base = None
    repo_name = "unknown"
    repo_accession = None

    for link in px_data.get("fullDatasetLinks", []):
        name = link.get("name", "")
        value = link.get("value", "")

        # Detect repository from dataset URI or FTP links
        if "jpost" in name.lower() or "jpost" in value.lower():
            repo_name = "jpost"
        elif "massive" in name.lower() or "massive" in value.lower():
            repo_name = "massive"

        # Capture FTP location (separate from dataset URI)
        if value.lower().startswith("ftp://"):
            ftp_base = value

    # Extract repository-specific accession from identifiers
    for ident in px_data.get("identifiers", []):
        name = ident.get("name", "")
        value = ident.get("value", "")
        if "jpost" in name.lower() and value:
            repo_accession = value  # e.g., "JPST000848"
            repo_name = "jpost"
        elif "massive" in name.lower() and value:
            repo_accession = value
            repo_name = "massive"

    # Fallback: extract repo accession from FTP URL
    if not repo_accession and ftp_base:
        match = re.search(r"(JPST\d+|MSV\d+)", ftp_base)
        if match:
            repo_accession = match.group(1)

    repo_info = {
        "ftp_base": ftp_base,
        "px_accession": accession,
        "repo_accession": repo_accession,
        "px_data": px_data,
    }

    print(f"   Resolved to {repo_name}: {repo_accession or 'unknown'}")
    if ftp_base:
        print(f"   FTP: {ftp_base}")

    return (repo_name, repo_info)


def _get_file_list_jpost(ftp_url: str) -> List[dict]:
    """
    List files from a jPOST FTP directory.

    Args:
        ftp_url: FTP URL (e.g., ftp://ftp.jpostdb.org/JPST000848/)

    Returns:
        List of file info dicts with fileName, fileSize, ftpLink
        (same format as PRIDE API file list).
    """
    import ftplib
    from urllib.parse import urlparse

    parsed = urlparse(ftp_url)
    host = parsed.hostname
    ftp_dir = parsed.path.rstrip("/")

    print(f"   Listing jPOST FTP: {host}{ftp_dir}")

    files = []
    try:
        ftp = ftplib.FTP(host, timeout=60)
        ftp.login()

        # Try to list the directory; jPOST may have subdirectories
        dirs_to_scan = [ftp_dir]

        # Some jPOST deposits put data in subdirectories
        try:
            entries = []
            ftp.dir(ftp_dir, entries.append)

            # Check for subdirectories that might contain raw data
            for entry in entries:
                parts = entry.split()
                if len(parts) < 9:
                    continue
                if entry.startswith("d"):  # directory
                    subdir_name = " ".join(parts[8:])
                    dirs_to_scan.append(f"{ftp_dir}/{subdir_name}")
        except Exception:
            pass

        for scan_dir in dirs_to_scan:
            try:
                entries = []
                ftp.dir(scan_dir, entries.append)
            except Exception:
                continue

            for entry in entries:
                parts = entry.split()
                if len(parts) < 9:
                    continue
                if entry.startswith("d"):
                    continue  # Skip directories
                size = int(parts[4])
                name = " ".join(parts[8:])
                ftp_link = f"ftp://{host}{scan_dir}/{name}"
                files.append({
                    "fileName": name,
                    "fileSize": size,
                    "ftpLink": ftp_link,
                })

        ftp.quit()
        print(f"   Found {len(files)} files via jPOST FTP")
    except Exception as e:
        print(f"   jPOST FTP listing failed: {e}")

    return files


def _reassemble_tdf_to_d(output_dir: Path) -> int:
    """
    Reassemble flat TDF files into Bruker .d folder structure.

    jPOST stores Bruker timsTOF data as flat .tdf + .tdf_bin file pairs.
    Search engines require .d folders containing analysis.tdf + analysis.tdf_bin.

    Args:
        output_dir: Directory containing downloaded TDF files

    Returns:
        Number of .d folders created
    """
    output_dir = Path(output_dir)

    # Find all .tdf files (not already inside .d folders)
    tdf_files = [
        f for f in output_dir.glob("*.tdf")
        if not any(p.suffix == ".d" for p in f.parents)
    ]

    if not tdf_files:
        return 0

    import shutil

    reassembled = 0
    for tdf_file in tdf_files:
        base_name = tdf_file.stem  # e.g., "181115ko03ProRP_P_dP_Fr02_Slot1-01_1_896"
        tdf_bin = output_dir / f"{base_name}.tdf_bin"

        if not tdf_bin.exists():
            print(f"   Warning: No .tdf_bin companion for {tdf_file.name}, skipping")
            continue

        # Create .d directory
        d_dir = output_dir / f"{base_name}.d"
        d_dir.mkdir(exist_ok=True)

        # Move files into .d as analysis.tdf / analysis.tdf_bin
        shutil.move(str(tdf_file), str(d_dir / "analysis.tdf"))
        shutil.move(str(tdf_bin), str(d_dir / "analysis.tdf_bin"))

        reassembled += 1

    if reassembled:
        print(f"   Reassembled {reassembled} .d folders from TDF files")

    return reassembled


def _create_px_metadata(accession: str, px_data: dict, metadata_dir: Path) -> None:
    """
    Create pride_metadata.yaml from ProteomeXchange data.

    Used when PRIDE API returns 404 and metadata comes from ProteomeXchange.

    Args:
        accession: PXD accession
        px_data: ProteomeXchange JSON response
        metadata_dir: Directory to write metadata
    """
    metadata_dir = Path(metadata_dir)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    metadata = DatasetMetadata(accession=accession)
    metadata.dataset_id = MetadataField.auto(accession, "proteomexchange")

    # Title
    title = px_data.get("title", "")
    if title:
        metadata.title = MetadataField.auto(title, "proteomexchange.title")

    # Description / abstract
    description = px_data.get("description", "")

    # Species from species list
    species_list = px_data.get("species", [])
    if species_list:
        org_names = []
        for sp in species_list:
            name = sp.get("name", "")
            if name:
                org_names.append(name)
        if org_names:
            metadata.organism = MetadataField.auto(
                org_names[0], "proteomexchange.species[0].name"
            )
            if len(org_names) > 1:
                metadata.is_multi_organism = True
                metadata.organisms_all = MetadataField.auto(
                    [{"name": n} for n in org_names], "proteomexchange.species"
                )

    # Instruments
    instruments = px_data.get("instruments", [])
    if instruments:
        instr_name = instruments[0].get("name", "")
        if instr_name:
            metadata.instrument = MetadataField.auto(
                instr_name, "proteomexchange.instruments[0].name"
            )

    # Contacts / lab
    contacts = px_data.get("contacts", [])
    if contacts:
        affiliation = contacts[0].get("affiliation", "")
        if affiliation:
            metadata.lab_id = MetadataField.auto(
                affiliation, "proteomexchange.contacts[0].affiliation"
            )

    # Publication DOI
    publications = px_data.get("publications", [])
    for pub in publications:
        doi = pub.get("doi", "")
        if doi:
            metadata.publication_doi = MetadataField.auto(
                doi, "proteomexchange.publications[0].doi"
            )
            break

    # Parse description for LC-MS hints
    protocol_text = f"{title}\n{description}"
    parsed = _parse_protocol_text(protocol_text)
    if parsed.get("acquisition_mode"):
        metadata.acquisition_mode = MetadataField.inferred(
            parsed["acquisition_mode"],
            "proteomexchange.description",
            confidence=parsed.get("mode_confidence", 0.7),
        )

    # Save as pride_metadata.yaml (same format, different source)
    metadata.to_yaml(metadata_dir / "pride_metadata.yaml")

    # Also save raw PX response
    raw_file = metadata_dir / "raw_api_response.json"
    with open(raw_file, "w") as f:
        json.dump(px_data, f, indent=2)


def _get_file_list_pride_api(accession: str) -> List[dict]:
    """
    List all files for a PRIDE project via the v2 REST API.

    Uses the paginated ``/projects/{accession}/files`` endpoint and normalises
    each entry to a dict with ``fileName``, ``fileSize``, ``fileCategory`` and
    ``ftpLink`` — where ``ftpLink`` is the FTP-protocol location already
    rewritten to the ftp.ebi.ac.uk HTTPS mirror.
    """
    files: List[dict] = []
    page = 0
    page_size = 100

    while True:
        url = (
            f"{PRIDE_API_BASE}/projects/{accession}/files"
            f"?pageSize={page_size}&page={page}"
        )
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        batch = resp.json()
        if not batch:
            break

        for f in batch:
            ftp_url = None
            for loc in f.get("publicFileLocations", []):
                if loc.get("name") == "FTP Protocol":
                    ftp_url = loc.get("value")
                    break
            files.append({
                "fileName": f.get("fileName", ""),
                "fileSize": f.get("fileSizeBytes", 0),
                "fileCategory": f.get("fileCategory", {}).get("value", ""),
                "ftpLink": to_ebi_https_url(ftp_url) if ftp_url else None,
            })

        if len(batch) < page_size:
            break
        page += 1

    return files


def get_file_list(accession: str) -> tuple:
    """
    Get list of files for a PRIDE accession.

    Returns tuple of (file_list, repo_name) where file_list is a list of
    file info dicts with fileName, fileSize, ftpLink, etc.

    Falls back through: PRIDE API → EBI HTTPS index → ProteomeXchange → jPOST.
    """
    # Primary: PRIDE v2 REST API project files endpoint.
    try:
        api_files = _get_file_list_pride_api(accession)
        if api_files:
            return (api_files, "pride")
        print("   PRIDE API returned no files, trying HTTPS index fallback...")
    except requests.RequestException as e:
        print(f"   PRIDE file API failed: {e}")

    # Fallback: scrape the EBI HTTPS directory index for the archive path.
    ftp_files = _get_file_list_ftp(accession)
    if ftp_files:
        return (ftp_files, "pride")

    # PRIDE didn't work — resolve via ProteomeXchange
    repo_name, repo_info = _resolve_repository(accession)

    if repo_name == "jpost" and repo_info.get("ftp_base"):
        files = _get_file_list_jpost(repo_info["ftp_base"])
        return (files, "jpost")

    # For other repositories, try using the FTP base URL directly
    if repo_info.get("ftp_base"):
        files = _get_file_list_jpost(repo_info["ftp_base"])  # Generic FTP listing
        return (files, repo_name)

    return ([], repo_name)


def _get_file_list_ftp(accession: str) -> List[dict]:
    """
    Fallback: list files by scraping the EBI HTTPS directory index.

    Discovers the archive path from the PRIDE API publication date and parses
    the Apache autoindex page at ftp.ebi.ac.uk (the HTTPS mirror), so it works
    behind HTTPS-only proxies where FTP is unavailable.
    """
    import html

    # Get publication date to construct the archive path
    proj_url = f"{PRIDE_API_BASE}/projects/{accession}"
    proj_resp = requests.get(proj_url, timeout=60)
    proj_resp.raise_for_status()
    proj_data = proj_resp.json()

    pub_date = proj_data.get("publicationDate", "")
    if not pub_date:
        print("   Cannot determine archive path: no publication date")
        return []

    year, month = pub_date.split("-")[:2]
    archive_dir = f"/pride/data/archive/{year}/{month}/{accession}"
    index_url = f"{EBI_HTTPS_BASE}{archive_dir}/"

    print(f"   HTTPS index: {index_url}")

    files = []
    try:
        resp = requests.get(index_url, timeout=60)
        resp.raise_for_status()

        # Apache autoindex: <a href="filename">. Skip column-sort links
        # (href="?...") and the parent-directory / sub-directory links
        # (absolute paths or trailing slash).
        seen = set()
        for match in re.finditer(r'href="([^"?/][^"]*)"', resp.text):
            name = html.unescape(match.group(1))
            if name.endswith("/") or name in seen:
                continue
            seen.add(name)
            files.append({
                "fileName": name,
                "fileSize": 0,  # not reliably available from the index page
                "ftpLink": f"{index_url}{name}",
            })

        print(f"   Found {len(files)} files via HTTPS index")
    except Exception as e:
        print(f"   HTTPS index listing failed: {e}")

    return files


def filter_raw_files(
    files: List[dict],
    patterns: Optional[List[str]] = None,
    max_files: int = 0,
) -> List[dict]:
    """
    Filter files to download based on patterns.

    Args:
        files: List of file info dicts from PRIDE API
        patterns: Glob patterns to match (e.g., ["*.d.zip", "*.d.tar"])
        max_files: Maximum number of files (0 = unlimited)

    Returns:
        Filtered list of files
    """
    if patterns is None:
        # Bruker timsTOF raw data: .d folders (usually archived) and loose TDF.
        patterns = ["*.d.zip", "*.d.tar", "*.d.tar.*", "*.d.rar", "*.d.7z",
                    "*.d", "*.tdf", "*.tdf_bin"]

    # Non-Bruker raw / derived formats. San José processes Bruker .d data only,
    # but some timsTOF-classified PRIDE datasets are multi-instrument and ship
    # foreign raw files — reject these even when tagged fileCategory "RAW".
    non_bruker = (".raw", ".wiff", ".wiff2", ".wiff.scan", ".lcd",
                  ".mzml", ".mzxml", ".mgf")

    import fnmatch

    filtered = []
    for f in files:
        filename = f.get("fileName", "")
        low = filename.lower()
        if low.endswith(non_bruker):
            continue
        # PRIDE tags raw-data files with fileCategory "RAW"; trust that over
        # filename globs so Bruker bundle archives (e.g. original_data.zip,
        # Raw_HeLa_Trp.zip) are still selected — the non-Bruker extension check
        # above has already excluded foreign raw files.
        if f.get("fileCategory", "") == "RAW":
            filtered.append(f)
            continue
        for pattern in patterns:
            if fnmatch.fnmatch(low, pattern.lower()):
                filtered.append(f)
                break

    # Sort by filename for consistent ordering
    filtered.sort(key=lambda x: x.get("fileName", ""))

    # Apply max files limit
    # For TDF pairs (.tdf + .tdf_bin), both files form one raw acquisition,
    # so limit by unique base names to get max_files acquisitions
    if max_files > 0:
        has_tdf = any(f.get("fileName", "").endswith(".tdf") for f in filtered)
        if has_tdf:
            # Group by base name and take max_files groups
            seen_bases = []
            tdf_filtered = []
            for f in filtered:
                name = f.get("fileName", "")
                if name.endswith(".tdf_bin"):
                    base = name[:-len(".tdf_bin")]
                elif name.endswith(".tdf"):
                    base = name[:-len(".tdf")]
                else:
                    base = name
                if base not in seen_bases:
                    if len(seen_bases) >= max_files:
                        break
                    seen_bases.append(base)
                tdf_filtered.append(f)
            filtered = tdf_filtered
        else:
            filtered = filtered[:max_files]

    return filtered


def download_with_wget(
    files: List[dict],
    output_dir: Path,
    retry_count: int = 3,
) -> bool:
    """
    Fallback download using wget for individual files.

    Args:
        files: List of file info dicts with ftpLink
        output_dir: Directory to save files
        retry_count: Number of retries per file

    Returns:
        True if all downloads succeeded
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    success_count = 0
    fail_count = 0

    for f in files:
        url = f.get("ftpLink")
        if not url:
            # Prefer the FTP-protocol location; Aspera locations cannot be
            # fetched over HTTPS. Order in publicFileLocations is not fixed.
            for loc in f.get("publicFileLocations", []):
                if loc.get("name") == "FTP Protocol":
                    url = loc.get("value")
                    break
        url = to_ebi_https_url(url) if url else None
        filename = f.get("fileName", "")

        if not url:
            print(f"Warning: No download URL for {filename}")
            fail_count += 1
            continue

        output_path = output_dir / filename

        # Skip if already downloaded
        if output_path.exists():
            expected_size = f.get("fileSize", 0)
            actual_size = output_path.stat().st_size
            if actual_size == expected_size:
                print(f"Skipping {filename} (already downloaded)")
                success_count += 1
                continue

        print(f"Downloading {filename}...")

        cmd = [
            "wget",
            "-c",   # Continue partial downloads
            "-nv",  # One concise line per file — no dotted progress spam in logs
            "-O", str(output_path),
            url,
        ]

        for attempt in range(retry_count):
            try:
                result = subprocess.run(cmd, timeout=3600)
                if result.returncode == 0:
                    success_count += 1
                    break
            except subprocess.TimeoutExpired:
                print(f"Download timed out for {filename}")
            except Exception as e:
                print(f"Download error for {filename}: {e}")

            if attempt == retry_count - 1:
                fail_count += 1

    print(f"Downloaded {success_count}/{len(files)} files ({fail_count} failed)")
    return fail_count == 0


def _flatten_d_wrappers(output_dir: Path) -> None:
    """Hoist Bruker .d folders nested one level inside a wrapper directory.

    Some PRIDE archives extract to a wrapper dir (e.g. ``MS_run_files.d/``) that
    itself contains the real ``.d`` acquisition folders. The runner globs
    ``output_dir/*.d`` flat, so nested .d folders are moved up. A genuine .d
    folder contains ``analysis.tdf`` directly; a wrapper does not.
    """
    import shutil

    def is_real_d(p: Path) -> bool:
        return p.is_dir() and (
            (p / "analysis.tdf").exists() or (p / "analysis.tdf_bin").exists()
        )

    for sub in list(output_dir.iterdir()):
        if not sub.is_dir() or is_real_d(sub):
            continue
        inner = [c for c in sub.iterdir() if c.name.endswith(".d") and is_real_d(c)]
        for d in inner:
            dest = output_dir / d.name
            if not dest.exists():
                shutil.move(str(d), str(dest))
                print(f"  flattened {sub.name}/{d.name} -> {d.name}")
        if inner and not any(c.is_dir() and c.name.endswith(".d") for c in sub.iterdir()):
            shutil.rmtree(sub, ignore_errors=True)


def extract_archives(output_dir: Path, cleanup: bool = True) -> None:
    """Extract downloaded Bruker .d archives (.zip / .tar* / .rar / .7z),
    flatten any nested .d wrappers, and delete archives after a clean extraction.
    """
    import shutil
    import zipfile
    import tarfile

    output_dir = Path(output_dir)
    extracted: List[Path] = []  # archives removed after a successful extraction

    # --- .zip ---
    for zip_file in output_dir.glob("*.zip"):
        print(f"Extracting {zip_file.name}...")
        try:
            with zipfile.ZipFile(zip_file, "r") as zf:
                names = [n for n in zf.namelist() if n]
                if zip_file.name.lower().endswith(".d.zip"):
                    d_name = zip_file.name[:-4]  # strip .zip, keep .d
                    has_parent = all(
                        n.replace("\\", "/").split("/")[0] == d_name for n in names
                    )
                    if has_parent:
                        zf.extractall(output_dir)
                    else:
                        (output_dir / d_name).mkdir(parents=True, exist_ok=True)
                        zf.extractall(output_dir / d_name)
                else:
                    zf.extractall(output_dir)
            extracted.append(zip_file)
        except Exception as e:
            print(f"  failed: {e}")

    # --- .tar / .tar.gz / .tar.bz2 / .tar.xz ---
    for tar_file in list(output_dir.glob("*.tar")) + list(output_dir.glob("*.tar.*")):
        print(f"Extracting {tar_file.name}...")
        try:
            with tarfile.open(tar_file, "r:*") as tf:
                tf.extractall(output_dir)
            extracted.append(tar_file)
        except Exception as e:
            print(f"  failed: {e}")

    # --- .rar ---
    rar_files = list(output_dir.glob("*.rar"))
    if rar_files:
        unrar_bin = shutil.which("unrar")
        if not unrar_bin:
            for cand in (Path.home() / "unrar" / "unrar",
                         Path.home() / ".local" / "bin" / "unrar"):
                if cand.exists():
                    unrar_bin = str(cand)
                    break
        if unrar_bin:
            for rar_file in rar_files:
                print(f"Extracting {rar_file.name}...")
                try:
                    r = subprocess.run(
                        [unrar_bin, "x", "-o+", str(rar_file), str(output_dir) + "/"],
                        capture_output=True, text=True, timeout=7200,
                    )
                    if r.returncode == 0:
                        extracted.append(rar_file)
                    else:
                        print(f"  unrar failed: {r.stderr[:200]}")
                except Exception as e:
                    print(f"  failed: {e}")
        else:
            print(f"Warning: unrar not found — skipping {len(rar_files)} .rar files")

    # --- .7z ---
    sz_files = list(output_dir.glob("*.7z"))
    if sz_files:
        sevenzip = next((shutil.which(n) for n in ("7z", "7za", "7zr")
                         if shutil.which(n)), None)
        for sz_file in sz_files:
            print(f"Extracting {sz_file.name}...")
            try:
                if sevenzip:
                    r = subprocess.run(
                        [sevenzip, "x", "-y", f"-o{output_dir}", str(sz_file)],
                        capture_output=True, text=True, timeout=7200,
                    )
                    if r.returncode == 0:
                        extracted.append(sz_file)
                    else:
                        print(f"  7z failed: {r.stderr[:200]}")
                else:
                    import py7zr
                    with py7zr.SevenZipFile(sz_file, "r") as z:
                        z.extractall(output_dir)
                    extracted.append(sz_file)
            except ImportError:
                print(f"  no 7z tool and no py7zr — skipping {sz_file.name}")
            except Exception as e:
                print(f"  failed: {e}")

    # --- normalise layout, then clean up archives that extracted cleanly ---
    _flatten_d_wrappers(output_dir)
    if cleanup:
        for arc in extracted:
            try:
                arc.unlink()
            except Exception:
                pass


def download_pride(
    accession: str,
    output_dir: Path,
    metadata_dir: Optional[Path] = None,
    retry_count: int = 3,
    max_files: int = 0,
    file_patterns: Optional[List[str]] = None,
    extract: bool = True,
) -> bool:
    """
    Download raw files from PRIDE Archive.

    Args:
        accession: PRIDE accession number (e.g., PXD019086)
        output_dir: Directory to save downloaded files
        metadata_dir: Directory for metadata files (default: data/metadata/{accession})
        retry_count: Number of retries for failed downloads
        max_files: Maximum number of files (0 = unlimited, for testing)
        file_patterns: Glob patterns for file selection
        extract: Whether to extract archives after download

    Returns:
        True if download succeeded
    """
    output_dir = Path(output_dir)
    if metadata_dir is None:
        metadata_dir = Path(f"data/metadata/{accession}")

    print(f"=" * 60)
    print(f"Download: {accession}")
    print(f"=" * 60)

    # Step 1: Fetch metadata
    print(f"\n1. Fetching metadata...")
    repo_name = "pride"  # Default, may be updated in Step 2

    try:
        metadata = fetch_pride_metadata(accession)
        # Save PRIDE metadata
        metadata_dir.mkdir(parents=True, exist_ok=True)
        pride_meta_file = metadata_dir / "pride_metadata.yaml"
        metadata.to_yaml(pride_meta_file)
        print(f"   Saved metadata to {pride_meta_file}")

        # Save raw API response
        try:
            raw_file = save_raw_response(accession, metadata_dir)
            print(f"   Saved raw API response to {raw_file}")
        except Exception as e:
            print(f"   Warning: Could not save raw API response: {e}")
    except requests.HTTPError as e:
        if e.response is not None and e.response.status_code == 404:
            print(f"   PRIDE API returned 404, will try ProteomeXchange...")
        else:
            print(f"   Warning: PRIDE metadata fetch failed: {e}")
    except Exception as e:
        print(f"   Warning: PRIDE metadata fetch failed: {e}")

    # Step 2: Get file list
    print(f"\n2. Getting file list...")
    try:
        all_files, repo_name = get_file_list(accession)
        print(f"   Found {len(all_files)} total files (source: {repo_name})")
    except Exception as e:
        print(f"   Error getting file list: {e}")
        return False

    # If resolved to non-PRIDE repo, create metadata from ProteomeXchange.
    # Overwrite if existing file is empty (PRIDE API returned 404 but
    # fetch_pride_metadata swallowed the error and wrote all-missing fields).
    if repo_name != "pride":
        px_metadata_file = metadata_dir / "pride_metadata.yaml"
        needs_px_metadata = not px_metadata_file.exists()
        if not needs_px_metadata:
            # Check if existing metadata is actually populated
            try:
                import yaml
                with open(px_metadata_file) as f:
                    existing = yaml.safe_load(f) or {}
                fields = existing.get("fields", {})
                has_data = any(
                    isinstance(v, dict) and v.get("status") != "missing"
                    for v in fields.values() if v
                )
                if not has_data:
                    needs_px_metadata = True
            except Exception:
                needs_px_metadata = True

        if needs_px_metadata:
            print(f"   Creating metadata from ProteomeXchange...")
            _, repo_info = _resolve_repository(accession)
            px_data = repo_info.get("px_data", {})
            if px_data:
                _create_px_metadata(accession, px_data, metadata_dir)
                print(f"   Saved ProteomeXchange metadata to {px_metadata_file}")

    # Filter files
    files = filter_raw_files(all_files, file_patterns, max_files)
    print(f"   Selected {len(files)} files for download")

    if not files:
        print("   No matching files found!")
        return False

    # Show files to download
    total_size = sum(f.get("fileSize", 0) for f in files)
    print(f"   Total size: {total_size / (1024**3):.2f} GB")
    for f in files[:5]:
        print(f"     - {f.get('fileName')}")
    if len(files) > 5:
        print(f"     ... and {len(files) - 5} more")

    # Step 3: Download files over HTTPS from the ftp.ebi.ac.uk mirror.
    # FTP, Aspera and pridepy are unreachable behind HTTPS-only allow-list
    # proxies (e.g. Mogon NHR), so every file is fetched via HTTPS instead.
    print(f"\n3. Downloading {len(files)} files via HTTPS...")
    output_dir.mkdir(parents=True, exist_ok=True)
    success = download_with_wget(files, output_dir, retry_count)

    if not success:
        print("   Download failed!")
        return False

    # Step 4: Extract archives
    if extract:
        print(f"\n4. Extracting archives...")
        extract_archives(output_dir)

    # Step 4b: Reassemble TDF files into .d folders (jPOST deposits)
    n_reassembled = _reassemble_tdf_to_d(output_dir)
    if n_reassembled > 0:
        print(f"   Reassembled {n_reassembled} .d folders from flat TDF files")

    # Step 5: Create download complete flag
    flag_file = output_dir / ".download_complete"
    flag_file.write_text(f"Downloaded at: {time.strftime('%Y-%m-%dT%H:%M:%S')}\n")
    print(f"\n5. Download complete!")
    print(f"   Files saved to: {output_dir}")
    print(f"   Metadata saved to: {metadata_dir}")

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Download raw data from PRIDE Archive"
    )
    parser.add_argument(
        "--accession", "-a",
        required=True,
        help="PRIDE accession number (e.g., PXD019086)"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        help="Output directory for raw files"
    )
    parser.add_argument(
        "--metadata-dir",
        type=Path,
        help="Directory for metadata files (default: data/metadata/{accession})"
    )
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="Only fetch metadata, don't download files"
    )
    parser.add_argument(
        "--retry",
        type=int,
        default=3,
        help="Number of download retries (default: 3)"
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=0,
        help="Maximum number of files to download (0 = all)"
    )
    parser.add_argument(
        "--no-extract",
        action="store_true",
        help="Don't extract archives after download"
    )

    args = parser.parse_args()

    # Set default output directory
    if args.output is None:
        args.output = Path(f"data/raw/{args.accession}")

    if args.metadata_dir is None:
        args.metadata_dir = Path(f"data/metadata/{args.accession}")

    if args.metadata_only:
        # Just fetch and save metadata
        print(f"Fetching metadata for {args.accession}...")
        metadata = fetch_pride_metadata(args.accession)

        args.metadata_dir.mkdir(parents=True, exist_ok=True)
        pride_meta_file = args.metadata_dir / "pride_metadata.yaml"
        metadata.to_yaml(pride_meta_file)
        print(f"Saved PRIDE metadata to {pride_meta_file}")

        # Also save raw API response
        try:
            raw_file = save_raw_response(args.accession, args.metadata_dir)
            print(f"Saved raw API response to {raw_file}")
        except Exception as e:
            print(f"Warning: Could not save raw API response: {e}")

        # Print summary
        summary = metadata.validation_summary()
        print(f"\nValidation summary:")
        print(f"  Auto fields: {summary['auto_fields']}")
        print(f"  Inferred fields: {summary['inferred_fields']}")
        print(f"  Missing fields: {summary['missing_fields']}")
        if summary["missing_field_names"]:
            print(f"  Missing: {', '.join(summary['missing_field_names'])}")

        return 0

    # Full download
    success = download_pride(
        accession=args.accession,
        output_dir=args.output,
        metadata_dir=args.metadata_dir,
        retry_count=args.retry,
        max_files=args.max_files,
        extract=not args.no_extract,
    )

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
