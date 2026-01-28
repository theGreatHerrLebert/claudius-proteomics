# Precursor Blob Design

Savable object capturing full precursor 4D data + merged fragment spectrum.

## Overview

```
PrecursorWithFragments
├── identifiers (precursor_id, charge, mono_mz, etc.)
├── precursor_signal: Precursor4D
│   ├── rt_profile (XIC)
│   ├── im_profile (mobilogram)
│   └── isotope_envelope
└── fragment_spectrum: TimsFrame (merged MS/MS)
```

## Class Structure

```python
from dataclasses import dataclass
from typing import Optional, List
import numpy as np
from numpy.typing import NDArray
import zstd
import json

from imspy_core.frame import TimsFrame
from imspy_core.dda import PrecursorDDA


@dataclass
class IsotopeEnvelope:
    """Isotope distribution at precursor m/z."""
    mz: NDArray[np.float64]           # Isotope m/z values (M, M+1, M+2, ...)
    intensity: NDArray[np.float64]    # Intensity per isotope

    def to_dict(self) -> dict:
        return {
            "mz": self.mz.tolist(),
            "intensity": self.intensity.tolist()
        }

    @classmethod
    def from_dict(cls, d: dict) -> "IsotopeEnvelope":
        return cls(
            mz=np.array(d["mz"], dtype=np.float64),
            intensity=np.array(d["intensity"], dtype=np.float64)
        )


@dataclass
class SignalProfile:
    """1D signal profile (XIC or mobilogram)."""
    coordinates: NDArray[np.float64]  # RT (seconds) or IM (1/K0)
    intensity: NDArray[np.float64]    # Signal at each coordinate
    apex: float                        # Apex position
    fwhm: float                        # Full width at half maximum

    def to_dict(self) -> dict:
        return {
            "coordinates": self.coordinates.tolist(),
            "intensity": self.intensity.tolist(),
            "apex": self.apex,
            "fwhm": self.fwhm
        }

    @classmethod
    def from_dict(cls, d: dict) -> "SignalProfile":
        return cls(
            coordinates=np.array(d["coordinates"], dtype=np.float64),
            intensity=np.array(d["intensity"], dtype=np.float64),
            apex=d["apex"],
            fwhm=d["fwhm"]
        )


@dataclass
class Precursor4D:
    """
    Full 4D precursor signal.

    Extracted from MS1 frames around the precursor event.
    """
    # Chromatographic dimension
    rt_profile: SignalProfile         # XIC (extracted ion chromatogram)

    # Ion mobility dimension
    im_profile: SignalProfile         # Mobilogram

    # m/z dimension
    isotope_envelope: IsotopeEnvelope # Isotope distribution

    # Apex coordinates (where precursor was selected)
    rt_apex: float                    # seconds
    im_apex: float                    # 1/K0
    mz_mono: float                    # monoisotopic m/z

    # Optional: raw 4D data cube (if we want full flexibility)
    raw_frame: Optional[TimsFrame] = None  # Full MS1 region around precursor

    def to_dict(self) -> dict:
        return {
            "rt_profile": self.rt_profile.to_dict(),
            "im_profile": self.im_profile.to_dict(),
            "isotope_envelope": self.isotope_envelope.to_dict(),
            "rt_apex": self.rt_apex,
            "im_apex": self.im_apex,
            "mz_mono": self.mz_mono,
            # raw_frame serialized separately if present
        }

    @classmethod
    def from_dict(cls, d: dict) -> "Precursor4D":
        return cls(
            rt_profile=SignalProfile.from_dict(d["rt_profile"]),
            im_profile=SignalProfile.from_dict(d["im_profile"]),
            isotope_envelope=IsotopeEnvelope.from_dict(d["isotope_envelope"]),
            rt_apex=d["rt_apex"],
            im_apex=d["im_apex"],
            mz_mono=d["mz_mono"],
            raw_frame=None  # Loaded separately
        )


@dataclass
class PrecursorWithFragments:
    """
    Complete precursor blob: identifiers + 4D signal + fragment spectrum.

    This is the atomic unit stored per precursor in the database.
    """

    # === Identifiers ===
    precursor_id: int                 # Unique ID within dataset
    raw_file: str                     # Source .d file

    # === Precursor properties ===
    charge: int
    mono_mz: float                    # Monoisotopic m/z
    isolation_mz: float               # Quadrupole isolation center
    isolation_width: float            # Quadrupole window width

    # === Precursor 4D signal ===
    precursor_signal: Precursor4D

    # === Fragment spectrum (merged if re-fragmented) ===
    fragment_frame: TimsFrame         # Full 4D MS/MS data
    n_fragments_merged: int           # How many MS/MS events were merged
    collision_energies: List[float]   # CE values used (one per merged spectrum)

    # === Derived values (computed, not stored) ===
    @property
    def ccs(self) -> float:
        """CCS computed from ion mobility apex."""
        # CCS = f(mz, charge, 1/K0) - use imspy calibration
        # TODO: implement using imspy_core calibration
        pass

    # === Serialization ===

    def to_bytes(self) -> bytes:
        """Serialize to compressed bytes for storage."""

        # 1. Metadata as JSON
        metadata = {
            "precursor_id": self.precursor_id,
            "raw_file": self.raw_file,
            "charge": self.charge,
            "mono_mz": self.mono_mz,
            "isolation_mz": self.isolation_mz,
            "isolation_width": self.isolation_width,
            "precursor_signal": self.precursor_signal.to_dict(),
            "n_fragments_merged": self.n_fragments_merged,
            "collision_energies": self.collision_energies,
        }
        metadata_bytes = json.dumps(metadata).encode("utf-8")

        # 2. Fragment frame as numpy arrays
        fragment_arrays = {
            "frame_id": np.array([self.fragment_frame.frame_id], dtype=np.int64),
            "rt": np.array([self.fragment_frame.retention_time], dtype=np.float64),
            "scan": self.fragment_frame.scan,
            "mobility": self.fragment_frame.mobility,
            "tof": self.fragment_frame.tof,
            "mz": self.fragment_frame.mz,
            "intensity": self.fragment_frame.intensity,
        }

        # 3. Pack into single buffer
        # Format: [metadata_len (4 bytes)][metadata][arrays as .npz]
        import io
        npz_buffer = io.BytesIO()
        np.savez_compressed(npz_buffer, **fragment_arrays)
        npz_bytes = npz_buffer.getvalue()

        # 4. Combine and compress
        metadata_len = len(metadata_bytes).to_bytes(4, "little")
        combined = metadata_len + metadata_bytes + npz_bytes

        return zstd.compress(combined)

    @classmethod
    def from_bytes(cls, data: bytes) -> "PrecursorWithFragments":
        """Deserialize from compressed bytes."""

        # 1. Decompress
        combined = zstd.decompress(data)

        # 2. Extract metadata
        metadata_len = int.from_bytes(combined[:4], "little")
        metadata_bytes = combined[4:4+metadata_len]
        metadata = json.loads(metadata_bytes.decode("utf-8"))

        # 3. Extract arrays
        import io
        npz_bytes = combined[4+metadata_len:]
        npz_buffer = io.BytesIO(npz_bytes)
        arrays = np.load(npz_buffer)

        # 4. Reconstruct TimsFrame
        fragment_frame = TimsFrame(
            frame_id=int(arrays["frame_id"][0]),
            ms_type="MS2",
            retention_time=float(arrays["rt"][0]),
            scan=arrays["scan"],
            mobility=arrays["mobility"],
            tof=arrays["tof"],
            mz=arrays["mz"],
            intensity=arrays["intensity"],
        )

        # 5. Reconstruct Precursor4D
        precursor_signal = Precursor4D.from_dict(metadata["precursor_signal"])

        return cls(
            precursor_id=metadata["precursor_id"],
            raw_file=metadata["raw_file"],
            charge=metadata["charge"],
            mono_mz=metadata["mono_mz"],
            isolation_mz=metadata["isolation_mz"],
            isolation_width=metadata["isolation_width"],
            precursor_signal=precursor_signal,
            fragment_frame=fragment_frame,
            n_fragments_merged=metadata["n_fragments_merged"],
            collision_energies=metadata["collision_energies"],
        )

    def summary_dict(self) -> dict:
        """Return summary for index parquet (no heavy data)."""
        return {
            "precursor_id": self.precursor_id,
            "raw_file": self.raw_file,
            "charge": self.charge,
            "mono_mz": self.mono_mz,
            "rt_apex": self.precursor_signal.rt_apex,
            "im_apex": self.precursor_signal.im_apex,
            "n_fragments_merged": self.n_fragments_merged,
            "n_peaks": len(self.fragment_frame.mz),
        }
```

## Storage Layout

```
{accession}/
├── manifest.json
├── fragpipe.parquet          # IDs from FragPipe
├── diann.parquet             # IDs from DIA-NN
├── spectra/
│   ├── index.parquet         # Lightweight index for queries
│   │   └── columns: precursor_id, raw_file, charge, mono_mz,
│   │                rt_apex, im_apex, n_peaks, blob_offset, blob_size
│   └── blobs.bin             # Concatenated PrecursorWithFragments.to_bytes()
```

**index.parquet** enables:
- Fast filtering by m/z, RT, charge without loading blobs
- Random access via (blob_offset, blob_size)

**blobs.bin**:
- Sequential writes during extraction
- Random reads via offset/size from index

## Extraction Workflow

```python
def extract_precursors(dataset: TimsDatasetDDA) -> Iterator[PrecursorWithFragments]:
    """Extract all precursors from a DDA dataset."""

    # Get all fragmented precursors
    precursor_df = dataset.get_fragmented_precursors()
    fragment_df = dataset.get_pasef_fragments()

    # Group fragments by precursor (handles re-fragmentation)
    for precursor_id, group in fragment_df.groupby("precursor_id"):
        precursor_meta = precursor_df[precursor_df.precursor_id == precursor_id].iloc[0]

        # Merge all fragment frames for this precursor
        fragment_frames = group["raw_data"].tolist()  # List[TimsFrame]
        merged_frame = fragment_frames[0]
        for frame in fragment_frames[1:]:
            merged_frame = merged_frame + frame  # TimsFrame.__add__

        # Extract precursor 4D signal (TODO: implement)
        precursor_signal = extract_precursor_4d(
            dataset=dataset,
            precursor_meta=precursor_meta,
        )

        yield PrecursorWithFragments(
            precursor_id=precursor_id,
            raw_file=dataset.path.name,
            charge=precursor_meta.charge or 0,
            mono_mz=precursor_meta.mono_mz,
            isolation_mz=precursor_meta.isolation_mz,
            isolation_width=precursor_meta.isolation_width,
            precursor_signal=precursor_signal,
            fragment_frame=merged_frame,
            n_fragments_merged=len(fragment_frames),
            collision_energies=group["collision_energy"].tolist(),
        )


def extract_precursor_4d(dataset: TimsDatasetDDA, precursor_meta: PrecursorDDA) -> Precursor4D:
    """
    Extract full 4D precursor signal from MS1 frames.

    TODO: Design and implement this.

    Needs to:
    1. Find MS1 frames around precursor RT
    2. Extract XIC (sum intensity at precursor m/z across RT)
    3. Extract mobilogram (sum intensity at precursor m/z across IM)
    4. Extract isotope envelope (m/z pattern at apex)
    """
    raise NotImplementedError("Precursor 4D extraction to be designed")
```

## Size Estimates

Per precursor blob (rough estimates):
- Metadata JSON: ~500 bytes
- Fragment frame (500 peaks): ~20 KB uncompressed → ~5 KB compressed
- **Total: ~5-10 KB per precursor**

Per dataset:
- 100K precursors × 10 KB = **~1 GB per dataset**

### Storage Strategy

Given ~1 GB per dataset, 1000 datasets = 1 TB just for blobs. Consider tiered storage:

| Tier | What | Storage | Use Case |
|------|------|---------|----------|
| **Always** | Index parquet | ~10 MB/dataset | Fast filtering, queries |
| **OVERLAP only** | Blobs for high-confidence | ~300 MB/dataset | Training data |
| **On-demand** | Full blobs | ~1 GB/dataset | Research, debugging |

**Recommendation:** Store blobs only for OVERLAP peptides (found by both FragPipe and DIA-NN). This reduces storage by ~60% while keeping the highest-confidence data.

## Design Decisions

### 1. Precursor 4D extraction window

**Decision:** Adaptive window based on peak shape from search engine output.

```python
# Use FWHM from FragPipe/DIA-NN as extraction window
rt_window = max(precursor.rt_fwhm * 3, 10.0)  # At least 10 seconds
im_window = max(precursor.im_fwhm * 3, 0.05)  # At least 0.05 1/K0
```

**Rationale:** Fixed windows miss narrow peaks or waste space on broad ones. The search engines already computed peak boundaries - reuse them.

### 2. Raw frame storage

**Decision:** Store profiles only, not raw MS1 frame.

**Rationale:**
- Profiles (XIC, mobilogram, isotope envelope) are sufficient for CCS/RT prediction
- Raw frame adds ~10x storage with marginal benefit
- Can always re-extract from .d files if needed for research

### 3. Fragment m/z calibration

**Decision:** Store calibrated m/z, not raw TOF indices.

**Rationale:**
- Self-contained blobs are easier to work with
- ~20% size increase is acceptable
- Avoids needing per-file calibration coefficients at read time
- Recalibration from TOF indices is error-prone across software versions

### 4. Mobilogram per fragment peak

**Decision:** Precursor-level mobilogram only.

**Rationale:**
- Fragment-level IM distributions would 10x the storage
- For CCS prediction, precursor mobilogram is what matters
- Fragment IM is useful for spectral library building (future work, separate pipeline)

## TODO

- [ ] Implement `extract_precursor_4d()` function with adaptive windowing
- [x] ~~Decide on RT/IM window strategy~~ → Adaptive (3× FWHM, minimum floor)
- [x] ~~Decide on raw TOF vs calibrated m/z~~ → Calibrated m/z (self-contained)
- [x] ~~Decide on raw MS1 frame vs extracted profiles~~ → Profiles only
- [ ] Benchmark blob sizes with real PXD019086 data
- [ ] Test serialization round-trip with actual TimsFrame objects
- [ ] Implement tiered storage (OVERLAP-only blobs vs full)
