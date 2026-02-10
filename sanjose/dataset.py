"""Dataset: top-level entry point for San José data.

Auto-discovers precursor_store.parquet, extracted/ blobs, and manifest.json
from a dataset directory (either merged/{accession}/ or a flat directory).
"""

import json
from pathlib import Path
from typing import Optional, List, Iterator, Dict, Any

import pandas as pd

from .store import PrecursorStore
from .blob import BlobReader
from .models import Precursor


class Dataset:
    """San José dataset: query precursors, iterate for training, access raw signal.

    Usage:
        ds = Dataset("data/merged/PXD019086")
        df = ds.query(sequence="YRPGTVALR", charge=2)
        for p in ds.iter(min_engines=3, limit=10):
            print(p.signal.fragment_spectrum.n_peaks)
    """

    def __init__(self, path: str):
        """Open a dataset directory.

        Auto-discovers:
          - precursor_store.parquet (required)
          - extracted/{accession}/ or extracted/ (optional, for blob access)
          - manifest.json (optional, for metadata)

        Args:
            path: Path to merged dataset directory, or parent data/ directory.
                  Examples:
                    "data/merged/PXD019086"
                    "/scratch/claudius-proteomics/PXD019086"
        """
        self._path = Path(path).resolve()
        self._store_path: Optional[Path] = None
        self._extracted_dir: Optional[Path] = None
        self._manifest_path: Optional[Path] = None

        self._discover()

        # Initialize blob reader if extracted dir found
        blob_reader = None
        if self._extracted_dir is not None:
            blob_reader = BlobReader(str(self._extracted_dir))

        self._store = PrecursorStore(str(self._store_path), blob_reader)
        self._manifest: Optional[Dict] = None

    def _discover(self):
        """Auto-discover store, blobs, and manifest from path."""
        p = self._path

        # Find precursor_store.parquet
        candidates = [
            p / "precursor_store.parquet",
            p / "merged" / "precursor_store.parquet",
        ]
        for c in candidates:
            if c.exists():
                self._store_path = c
                break

        if self._store_path is None:
            # Try to find in any subdirectory matching the accession pattern
            for f in p.rglob("precursor_store.parquet"):
                self._store_path = f
                break

        if self._store_path is None:
            raise FileNotFoundError(
                f"precursor_store.parquet not found in {p}. "
                f"Expected at {p}/precursor_store.parquet"
            )

        # Find extracted/ directory (for blobs)
        store_parent = self._store_path.parent
        accession = self._guess_accession(store_parent)

        # Common layouts:
        #   data/merged/PXD019086/precursor_store.parquet → data/extracted/PXD019086/
        #   PXD019086/precursor_store.parquet → PXD019086/extracted/
        extracted_candidates = [
            store_parent.parent.parent / "extracted" / accession if accession else None,
            store_parent.parent / "extracted" / accession if accession else None,
            store_parent / "extracted",
            p / "extracted" / accession if accession else None,
            p / "extracted",
            p.parent / "extracted" / accession if accession else None,
        ]
        for c in extracted_candidates:
            if c is not None and c.exists() and c.is_dir():
                self._extracted_dir = c
                break

        # Find manifest.json
        manifest_candidates = [
            store_parent / "manifest.json",
            p / "manifest.json",
        ]
        for c in manifest_candidates:
            if c.exists():
                self._manifest_path = c
                break

    def _guess_accession(self, directory: Path) -> Optional[str]:
        """Guess PRIDE accession from directory name (PXDxxxxxx)."""
        name = directory.name
        if name.startswith("PXD") and len(name) >= 9:
            return name
        # Check parent
        parent_name = directory.parent.name
        if parent_name.startswith("PXD") and len(parent_name) >= 9:
            return parent_name
        return None

    @property
    def accession(self) -> Optional[str]:
        """PRIDE accession (from manifest or directory name)."""
        if self.manifest and 'accession' in self.manifest:
            return self.manifest['accession']
        return self._guess_accession(self._store_path.parent)

    @property
    def num_precursors(self) -> int:
        return self._store.num_precursors

    @property
    def manifest(self) -> Optional[Dict]:
        """Parsed manifest.json metadata."""
        if self._manifest is None and self._manifest_path is not None:
            with open(self._manifest_path) as f:
                self._manifest = json.load(f)
        return self._manifest

    @property
    def has_blobs(self) -> bool:
        return self._store.blob_reader is not None

    @property
    def store(self) -> PrecursorStore:
        """Direct access to the underlying PrecursorStore."""
        return self._store

    # Delegate query methods to store

    def get(
        self, precursor_id: int, raw_file: Optional[str] = None,
    ) -> Optional[Precursor]:
        """Get a single precursor by ID (+ raw_file for uniqueness)."""
        return self._store.get(precursor_id, raw_file)

    def query(
        self,
        *,
        columns: Optional[List[str]] = None,
        limit: Optional[int] = None,
        **filter_kwargs,
    ) -> pd.DataFrame:
        """Query precursors → DataFrame. See PrecursorStore.query()."""
        return self._store.query(columns=columns, limit=limit, **filter_kwargs)

    def iter(
        self,
        *,
        limit: Optional[int] = None,
        **filter_kwargs,
    ) -> Iterator[Precursor]:
        """Iterate precursors as Precursor objects. See PrecursorStore.iter()."""
        return self._store.iter(limit=limit, **filter_kwargs)

    def batches(
        self,
        batch_size: int = 1000,
        *,
        columns: Optional[List[str]] = None,
        include_signal: bool = False,
        **filter_kwargs,
    ) -> Iterator[pd.DataFrame]:
        """Streaming batch iteration. See PrecursorStore.batches()."""
        return self._store.batches(
            batch_size, columns=columns, include_signal=include_signal,
            **filter_kwargs,
        )

    def summary(self) -> Dict[str, Any]:
        """Dataset summary with manifest metadata."""
        s = self._store.summary()
        s['accession'] = self.accession
        s['path'] = str(self._path)
        if self.manifest:
            s['pipeline_version'] = self.manifest.get('pipeline_version')
            s['generated_at'] = self.manifest.get('generated_at')
        return s

    def __len__(self) -> int:
        return self.num_precursors

    def __repr__(self) -> str:
        acc = self.accession or "unknown"
        blobs = ", with blobs" if self.has_blobs else ""
        return f"<Dataset {acc}: {self.num_precursors:,} precursors{blobs}>"
