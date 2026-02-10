"""BlobReader: reads and decompresses raw 4D signal from blobs.bin files.

Blob format (per precursor):
  [zstd or gzip compressed payload]
  After decompression:
    [4-byte little-endian metadata_len]
    [metadata JSON (metadata_len bytes)]
    [NPZ arrays]

Arrays in NPZ:
  - frag_mz, frag_intensity, frag_mobility, frag_scan (fragment spectrum)
  - ms1_rt_coords, ms1_rt_intensities (XIC)
  - ms1_im_coords, ms1_im_intensities (mobilogram)
  - raw_rt, raw_mz, raw_mobility, raw_intensity (4D point cloud)
"""

import gzip
import io
import json
from functools import lru_cache
from pathlib import Path
from typing import Optional, List

import numpy as np

from .models import FragmentSpectrum, MS1Signal, RawPointCloud, RawSignal
from ._compat import HAS_ZSTD, get_zstd_decompressor


class BlobReader:
    """Reads raw 4D signal data from blobs.bin files."""

    def __init__(self, extracted_dir: str):
        """Initialize with path to extracted/ directory.

        Args:
            extracted_dir: Path to directory containing {raw_file}.d/blobs.bin
        """
        self.extracted_dir = Path(extracted_dir)

    @lru_cache(maxsize=64)
    def _resolve_blob_path(self, raw_file: str) -> Optional[Path]:
        """Find the blobs.bin file for a given raw file name."""
        raw_clean = raw_file.replace('.d', '')

        candidates = [
            self.extracted_dir / f"{raw_clean}.d" / "blobs.bin",
            self.extracted_dir / raw_file / "blobs.bin",
        ]
        for p in candidates:
            if p.exists():
                return p
        return None

    def _decompress(self, data: bytes) -> bytes:
        """Detect compression and decompress."""
        # zstd magic: 0x28 0xB5 0x2F 0xFD
        if data[:4] == b'\x28\xb5\x2f\xfd':
            dctx = get_zstd_decompressor()
            return dctx.decompress(data)
        # gzip magic: 0x1F 0x8B
        elif data[:2] == b'\x1f\x8b':
            return gzip.decompress(data)
        else:
            raise ValueError(f"Unknown compression format: {data[:4].hex()}")

    def _parse_blob(self, combined: bytes) -> RawSignal:
        """Parse decompressed blob into RawSignal."""
        metadata_len = int.from_bytes(combined[:4], "little")
        metadata = json.loads(combined[4:4 + metadata_len].decode("utf-8"))

        npz_bytes = combined[4 + metadata_len:]
        arrays = np.load(io.BytesIO(npz_bytes))

        def get(name: str) -> np.ndarray:
            return arrays[name] if name in arrays else np.array([], dtype=np.float32)

        fragment_spectrum = FragmentSpectrum(
            mz=get("frag_mz"),
            intensity=get("frag_intensity"),
            mobility=get("frag_mobility"),
            scan=get("frag_scan") if "frag_scan" in arrays else None,
        )

        ms1_signal = MS1Signal(
            xic_rt=get("ms1_rt_coords"),
            xic_intensity=get("ms1_rt_intensities"),
            mobilogram_im=get("ms1_im_coords"),
            mobilogram_intensity=get("ms1_im_intensities"),
            isotope_intensity=np.array(
                metadata.get("ms1_isotope_intensities", []), dtype=np.float64
            ),
        )

        raw_point_cloud = RawPointCloud(
            rt=get("raw_rt"),
            mz=get("raw_mz"),
            mobility=get("raw_mobility"),
            intensity=get("raw_intensity"),
        )

        return RawSignal(
            fragment_spectrum=fragment_spectrum,
            ms1_signal=ms1_signal,
            raw_point_cloud=raw_point_cloud,
            metadata=metadata,
        )

    def read(self, raw_file: str, offset: int, size: int) -> Optional[RawSignal]:
        """Read and decompress a single precursor blob.

        Args:
            raw_file: Raw file name (with or without .d suffix)
            offset: Byte offset into blobs.bin
            size: Number of bytes to read

        Returns:
            RawSignal or None if blob file not found
        """
        blob_path = self._resolve_blob_path(raw_file)
        if blob_path is None:
            return None

        with open(blob_path, 'rb') as f:
            f.seek(offset)
            compressed = f.read(size)

        combined = self._decompress(compressed)
        return self._parse_blob(combined)

    def read_batch(
        self, raw_file: str, offsets: List[int], sizes: List[int]
    ) -> List[Optional[RawSignal]]:
        """Read multiple blobs from the same raw file efficiently.

        Opens the blob file once and seeks to each offset.

        Args:
            raw_file: Raw file name
            offsets: List of byte offsets
            sizes: List of byte sizes

        Returns:
            List of RawSignal (or None for missing blobs)
        """
        blob_path = self._resolve_blob_path(raw_file)
        if blob_path is None:
            return [None] * len(offsets)

        results = []
        with open(blob_path, 'rb') as f:
            for offset, size in zip(offsets, sizes):
                try:
                    f.seek(offset)
                    compressed = f.read(size)
                    combined = self._decompress(compressed)
                    results.append(self._parse_blob(combined))
                except Exception:
                    results.append(None)

        return results

    def has_blobs(self, raw_file: str) -> bool:
        """Check if blob file exists for this raw file."""
        return self._resolve_blob_path(raw_file) is not None

    def __repr__(self) -> str:
        return f"<BlobReader: {self.extracted_dir}>"
