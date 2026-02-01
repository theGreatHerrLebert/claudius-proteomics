#!/usr/bin/env python3
"""
Zarr-based Precursor Store

Memory-mapped storage for extracted precursor data:
- Fragment spectra (m/z, intensity arrays)
- MS1 signals (XIC, mobilogram, isotope envelope)
- Fast random access without full decompression

Storage layout:
  precursors.zarr/
    ├── index/                    # Metadata (copied from parquet for fast access)
    │   ├── precursor_id
    │   ├── raw_file
    │   ├── mz, charge, rt, mobility
    │   └── ...
    ├── fragments/                # Fragment spectra per precursor
    │   ├── mz/                   # Ragged arrays stored as chunks
    │   ├── intensity/
    │   ├── mobility/
    │   └── offsets              # Index into ragged arrays
    └── ms1/                      # MS1 precursor signals
        ├── xic_rt/              # XIC retention time coords (1D projection)
        ├── xic_intensity/
        ├── mobilogram_im/       # Mobilogram coords (1D projection)
        ├── mobilogram_intensity/
        ├── isotope_mz/          # Isotope envelope (1D projection)
        ├── isotope_intensity/
        ├── raw_rt/              # Raw 4D data: RT per peak
        ├── raw_mz/              # Raw 4D data: m/z per peak
        ├── raw_mobility/        # Raw 4D data: 1/K0 per peak
        ├── raw_intensity/       # Raw 4D data: intensity per peak
        └── *_offsets            # Index into ragged arrays

Usage:
    # Create store from extracted data
    store = PrecursorStore.create_from_extraction(
        raw_data_path="/path/to/data.d",
        output_path="/path/to/precursors.zarr",
    )

    # Open existing store
    store = PrecursorStore("/path/to/precursors.zarr")

    # Get precursor data
    precursor = store.get_precursor(12345)
    print(precursor.fragment_mz)
    print(precursor.xic_rt)
"""

import json
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, List, Dict, Any, Iterator
import numpy as np
import pandas as pd
import zarr

# Handle Zarr v2/v3 codec differences
try:
    from zarr.codecs import BloscCodec
    USE_ZARR_V3 = True
except ImportError:
    from numcodecs import Blosc
    USE_ZARR_V3 = False


@dataclass
class PrecursorData:
    """Complete data for a single precursor."""
    # Identifiers
    precursor_id: int
    raw_file: str

    # Properties
    mz: float
    charge: int
    rt_seconds: float
    mobility: float

    # Search engine IDs (may be None)
    fragpipe_peptide: Optional[str] = None
    diann_peptide: Optional[str] = None
    sage_peptide: Optional[str] = None
    consensus_peptide: Optional[str] = None
    n_engines: int = 0

    # Fragment spectrum
    fragment_mz: Optional[np.ndarray] = None
    fragment_intensity: Optional[np.ndarray] = None
    fragment_mobility: Optional[np.ndarray] = None

    # MS1 signals (1D projections)
    xic_rt: Optional[np.ndarray] = None
    xic_intensity: Optional[np.ndarray] = None
    mobilogram_im: Optional[np.ndarray] = None
    mobilogram_intensity: Optional[np.ndarray] = None
    isotope_mz: Optional[np.ndarray] = None
    isotope_intensity: Optional[np.ndarray] = None

    # Raw 4D MS1 data (all peaks from filtered MS1 frames in RT window, merged)
    # Each peak has (rt, mz, mobility) coordinates + intensity value
    raw_rt: Optional[np.ndarray] = None         # RT per peak (seconds)
    raw_mz: Optional[np.ndarray] = None         # m/z per peak
    raw_mobility: Optional[np.ndarray] = None   # 1/K0 per peak
    raw_intensity: Optional[np.ndarray] = None  # intensity per peak

    # Quality metrics
    fragpipe_probability: Optional[float] = None
    diann_qvalue: Optional[float] = None
    sage_qvalue: Optional[float] = None


class PrecursorStore:
    """Zarr-based precursor data store with memory-mapped access."""

    def __init__(self, store_path: str):
        """Open existing store."""
        self.store_path = Path(store_path)
        self.root = zarr.open(str(self.store_path), mode='r')

        # Load index as DataFrame for fast queries
        self._load_index()

    def _load_index(self):
        """Load index arrays into DataFrame."""
        index_group = self.root['index']

        # Build DataFrame from index arrays
        data = {}
        for key in index_group.array_keys():
            arr = index_group[key][:]
            # Decode bytes to strings if needed
            if arr.dtype.kind == 'S':
                arr = np.char.decode(arr, 'utf-8')
            data[key] = arr

        self.index_df = pd.DataFrame(data)

        # Build precursor_id -> row index lookup
        self._id_to_idx = {
            pid: idx for idx, pid in enumerate(self.index_df['precursor_id'])
        }

    def __len__(self):
        return len(self.index_df)

    def get_precursor(self, precursor_id: int) -> Optional[PrecursorData]:
        """Get complete data for a precursor by ID."""
        if precursor_id not in self._id_to_idx:
            return None

        idx = self._id_to_idx[precursor_id]
        row = self.index_df.iloc[idx]

        # Get fragment data
        frag_mz, frag_int, frag_mob = self._get_fragment_arrays(idx)

        # Get MS1 data with correct offsets
        xic_rt, xic_int = self._get_array_pair('ms1', 'xic_rt', 'xic_intensity', idx, 'xic_offsets')
        mob_im, mob_int = self._get_array_pair('ms1', 'mobilogram_im', 'mobilogram_intensity', idx, 'mobilogram_offsets')
        iso_mz, iso_int = self._get_array_pair('ms1', 'isotope_mz', 'isotope_intensity', idx, 'isotope_offsets')

        # Get raw 2D MS1 data
        raw_rt, raw_mz, raw_mobility, raw_intensity = self._get_raw_ms1_arrays(idx)

        return PrecursorData(
            precursor_id=int(row['precursor_id']),
            raw_file=str(row.get('raw_file', '')),
            mz=float(row.get('mz', row.get('fragpipe_mz', 0))),
            charge=int(row.get('charge', 2) if not pd.isna(row.get('charge', np.nan)) else row.get('fragpipe_charge', 2) if not pd.isna(row.get('fragpipe_charge', np.nan)) else 2),
            rt_seconds=float(row.get('rt_seconds', row.get('fragpipe_rt', 0) * 60)),
            mobility=float(row.get('mobility', row.get('fragpipe_mobility', 0))),
            fragpipe_peptide=row.get('fragpipe_peptide'),
            diann_peptide=row.get('diann_peptide'),
            sage_peptide=row.get('sage_peptide'),
            consensus_peptide=row.get('consensus_peptide'),
            n_engines=int(row.get('n_engines', 0)),
            fragment_mz=frag_mz,
            fragment_intensity=frag_int,
            fragment_mobility=frag_mob,
            xic_rt=xic_rt,
            xic_intensity=xic_int,
            mobilogram_im=mob_im,
            mobilogram_intensity=mob_int,
            isotope_mz=iso_mz,
            isotope_intensity=iso_int,
            raw_rt=raw_rt,
            raw_mz=raw_mz,
            raw_mobility=raw_mobility,
            raw_intensity=raw_intensity,
            fragpipe_probability=row.get('fragpipe_probability'),
            diann_qvalue=row.get('diann_qvalue'),
            sage_qvalue=row.get('sage_qvalue'),
        )

    def _get_fragment_arrays(self, idx: int):
        """Get fragment spectrum arrays for precursor at index."""
        if 'fragments' not in self.root:
            return None, None, None

        frag = self.root['fragments']
        offsets = frag['offsets'][:]

        start = offsets[idx]
        end = offsets[idx + 1] if idx + 1 < len(offsets) else len(frag['mz'])

        if start >= end:
            return np.array([]), np.array([]), np.array([])

        mz = frag['mz'][start:end]
        intensity = frag['intensity'][start:end]
        mobility = frag['mobility'][start:end] if 'mobility' in frag else None

        return mz, intensity, mobility

    def _get_array_pair(self, group_name: str, coord_name: str, value_name: str, idx: int, offset_name: str = None):
        """Get a coordinate-value array pair for precursor at index."""
        if group_name not in self.root:
            return None, None

        group = self.root[group_name]
        if coord_name not in group or value_name not in group:
            return None, None

        # Determine which offset array to use
        if offset_name is None:
            # Try specific offset first, fall back to generic
            if f'{coord_name.split("_")[0]}_offsets' in group:
                offset_name = f'{coord_name.split("_")[0]}_offsets'
            else:
                offset_name = 'offsets'

        if offset_name not in group:
            return None, None

        offsets = group[offset_name][:]
        start = offsets[idx]
        end = offsets[idx + 1] if idx + 1 < len(offsets) else len(group[coord_name])

        if start >= end:
            return np.array([]), np.array([])

        coords = group[coord_name][start:end]
        values = group[value_name][start:end]

        return coords, values

    def _get_raw_ms1_arrays(self, idx: int):
        """Get raw 4D MS1 arrays (rt, mz, mobility, intensity) for precursor at index."""
        if 'ms1' not in self.root:
            return None, None, None, None

        ms1 = self.root['ms1']
        if 'raw_rt' not in ms1:
            return None, None, None, None

        offsets = ms1['raw_offsets'][:]
        start = offsets[idx]
        end = offsets[idx + 1] if idx + 1 < len(offsets) else len(ms1['raw_rt'])

        if start >= end:
            return np.array([]), np.array([]), np.array([]), np.array([])

        raw_rt = ms1['raw_rt'][start:end]
        raw_mz = ms1['raw_mz'][start:end]
        raw_mobility = ms1['raw_mobility'][start:end]
        raw_intensity = ms1['raw_intensity'][start:end]

        return raw_rt, raw_mz, raw_mobility, raw_intensity

    def query(self, **kwargs) -> pd.DataFrame:
        """Query index by column values."""
        mask = pd.Series(True, index=self.index_df.index)

        for col, value in kwargs.items():
            if col in self.index_df.columns:
                if isinstance(value, (list, tuple)):
                    mask &= self.index_df[col].isin(value)
                else:
                    mask &= self.index_df[col] == value

        return self.index_df[mask]

    def get_by_agreement(self, min_engines: int = 1) -> pd.DataFrame:
        """Get precursors sorted by agreement (n_engines descending)."""
        df = self.index_df[self.index_df['n_engines'] >= min_engines].copy()
        return df.sort_values(['n_engines', 'fragpipe_probability'],
                             ascending=[False, False])

    def iter_precursors(self, precursor_ids: Optional[List[int]] = None) -> Iterator[PrecursorData]:
        """Iterate over precursors."""
        if precursor_ids is None:
            precursor_ids = self.index_df['precursor_id'].tolist()

        for pid in precursor_ids:
            data = self.get_precursor(pid)
            if data:
                yield data

    @classmethod
    def create_from_index_and_raw(
        cls,
        index_path: str,
        raw_data_path: str,
        output_path: str,
        num_threads: int = 16,
        rt_window_sec: float = 30.0,
        mz_tol_ppm: float = 20.0,
        im_window: float = 0.1,
        batch_size: int = 2000,
        limit: Optional[int] = None,
    ) -> 'PrecursorStore':
        """
        Create store from unified index + raw data extraction.

        Args:
            index_path: Path to precursor_index.parquet
            raw_data_path: Path to .d folder
            output_path: Where to create .zarr store
            num_threads: Threads for extraction
            rt_window_sec: RT window for MS1 extraction
            mz_tol_ppm: m/z tolerance for MS1
            im_window: Ion mobility window

        Returns:
            PrecursorStore instance
        """
        from imspy_core.timstof import TimsDatasetDDA
        from imspy_connector import py_dda

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Load index
        print(f"Loading index: {index_path}")
        index_df = pd.read_parquet(index_path)

        # Filter to this raw file
        raw_name = Path(raw_data_path).stem
        file_index = index_df[index_df['raw_file'] == raw_name].copy()
        print(f"  {len(file_index)} precursors for {raw_name}")

        if len(file_index) == 0:
            raise ValueError(f"No precursors found for raw file: {raw_name}")

        # Apply limit for testing
        if limit is not None and limit < len(file_index):
            print(f"  Limiting to {limit} precursors")
            file_index = file_index.head(limit)

        # Load dataset
        print(f"Loading dataset: {raw_data_path}")
        dataset = TimsDatasetDDA(str(raw_data_path), in_memory=False, use_bruker_sdk=False)

        # Get fragment data
        print("Extracting fragments...")
        fragments_df = dataset.get_pasef_fragments(num_threads=num_threads)

        # Group by precursor_id
        grouped = fragments_df.groupby('precursor_id').agg({
            'raw_data': 'sum',
            'time': 'first',
        })

        # Build precursor coordinates for MS1 extraction
        print("Preparing MS1 extraction...")
        coords = []
        precursor_data = {}

        for _, row in file_index.iterrows():
            pid = int(row['precursor_id'])

            # Get RT from grouped fragments
            if pid in grouped.index:
                rt_sec = grouped.loc[pid, 'time'] * 60.0
                fragment_frame = grouped.loc[pid, 'raw_data']
            else:
                rt_val = row.get('fragpipe_rt') or 0
                if pd.isna(rt_val):
                    rt_val = 0
                rt_sec = rt_val * 60.0
                fragment_frame = None

            # Get precursor properties - prefer search engine values, fall back to raw metadata
            mz = row.get('fragpipe_mz') or row.get('raw_mz') or 0
            if pd.isna(mz):
                mz = 0
            charge_val = row.get('fragpipe_charge') or row.get('raw_charge') or 2
            if pd.isna(charge_val):
                charge_val = 2
            charge = int(charge_val)
            mobility = row.get('fragpipe_mobility') or row.get('raw_mobility') or 0.9
            if pd.isna(mobility):
                mobility = 0.9

            if mz > 0:
                coord = py_dda.PyPrecursorCoord(
                    precursor_id=pid,
                    mz=float(mz),
                    rt_seconds=float(rt_sec),
                    mobility=float(mobility),
                    charge=charge,
                )
                coords.append(coord)

            precursor_data[pid] = {
                'fragment_frame': fragment_frame,
                'rt_seconds': rt_sec,
            }

        # Create Zarr store first (before extraction)
        print(f"Creating Zarr store: {output_path}")
        if USE_ZARR_V3:
            compressors = BloscCodec(cname='zstd', clevel=3, shuffle='bitshuffle')
        else:
            compressors = Blosc(cname='zstd', clevel=3, shuffle=Blosc.BITSHUFFLE)
        root = zarr.open(str(output_path), mode='w')

        # Write index
        index_group = root.create_group('index')
        for col in file_index.columns:
            arr = file_index[col].values
            # Handle various pandas dtypes
            dtype_name = str(arr.dtype)
            if arr.dtype == object:
                # Convert strings to bytes for zarr
                arr = np.array([str(x) if pd.notna(x) else '' for x in arr], dtype='S256')
            elif 'Int' in dtype_name or 'int' in dtype_name.lower():
                # Convert pandas nullable integers to numpy int64
                arr = pd.array(arr).to_numpy(dtype='float64', na_value=np.nan)
                # Fill NaN with -1 and convert to int64
                arr = np.where(np.isnan(arr), -1, arr).astype(np.int64)
            elif 'Float' in dtype_name or 'float' in dtype_name.lower():
                # Convert pandas nullable floats to numpy float64
                arr = pd.array(arr).to_numpy(dtype='float64', na_value=np.nan)
            index_group.create_array(col, data=arr, compressors=compressors)

        # Write fragments (ragged arrays)
        frag_group = root.create_group('fragments')
        all_mz = []
        all_intensity = []
        all_mobility = []
        offsets = [0]

        for pid in file_index['precursor_id']:
            pdata = precursor_data.get(pid, {})
            frame = pdata.get('fragment_frame')

            if frame is not None and len(frame.mz) > 0:
                all_mz.extend(frame.mz)
                all_intensity.extend(frame.intensity)
                all_mobility.extend(frame.mobility)
                offsets.append(len(all_mz))
            else:
                offsets.append(offsets[-1])

        frag_group.create_array('mz', data=np.array(all_mz, dtype=np.float64), compressors=compressors)
        frag_group.create_array('intensity', data=np.array(all_intensity, dtype=np.float64), compressors=compressors)
        frag_group.create_array('mobility', data=np.array(all_mobility, dtype=np.float64), compressors=compressors)
        frag_group.create_array('offsets', data=np.array(offsets, dtype=np.int64), compressors=compressors)
        n_fragment_peaks = len(all_mz)

        # Clear fragment data from memory
        del precursor_data, all_mz, all_intensity, all_mobility, offsets
        import gc
        gc.collect()

        # Extract MS1 signals in batches (raw 4D data is very large)
        print(f"Extracting MS1 signals for {len(coords)} precursors in batches of {batch_size}...")
        rust_dataset = dataset.get_py_ptr()

        # Accumulators for all MS1 data
        xic_rt = []
        xic_int = []
        mob_im = []
        mob_int = []
        iso_mz = []
        iso_int = []
        offsets_xic = [0]
        offsets_mob = [0]
        offsets_iso = [0]
        raw_rt_all = []
        raw_mz_all = []
        raw_mobility_all = []
        raw_intensity_all = []
        raw_offsets = [0]

        # Build coord lookup for ordering
        coord_lookup = {c.precursor_id: c for c in coords}
        precursor_ids = list(file_index['precursor_id'])

        # Build mapping: precursor_id -> batch index
        precursor_batch = {}
        for i, c in enumerate(coords):
            precursor_batch[c.precursor_id] = i // batch_size

        # Process in batches, storing results temporarily
        n_batches = (len(coords) + batch_size - 1) // batch_size
        all_ms1_lookup = {}

        for batch_idx in range(n_batches):
            batch_start = batch_idx * batch_size
            batch_end = min(batch_start + batch_size, len(coords))
            batch_coords = coords[batch_start:batch_end]

            print(f"  Batch {batch_idx + 1}/{n_batches}: {len(batch_coords)} precursors...")

            # Extract MS1 signals for this batch
            batch_signals = rust_dataset.extract_precursor_ms1_signals(
                batch_coords,
                rt_window_sec=rt_window_sec,
                mz_tol_ppm=mz_tol_ppm,
                im_window=im_window,
                n_isotopes=5,
                num_threads=num_threads,
            )

            # Store only the arrays we need, then clear the signal objects
            for sig in batch_signals:
                all_ms1_lookup[sig.precursor_id] = {
                    'rt_coords': np.asarray(sig.rt_coords),
                    'rt_intensities': np.asarray(sig.rt_intensities),
                    'im_coords': np.asarray(sig.im_coords),
                    'im_intensities': np.asarray(sig.im_intensities),
                    'isotope_mz': np.asarray(sig.isotope_mz),
                    'isotope_intensity': np.asarray(sig.isotope_intensity),
                    'raw_rt': np.asarray(sig.raw_rt) if hasattr(sig, 'raw_rt') else np.array([]),
                    'raw_mz': np.asarray(sig.raw_mz) if hasattr(sig, 'raw_mz') else np.array([]),
                    'raw_mobility': np.asarray(sig.raw_mobility) if hasattr(sig, 'raw_mobility') else np.array([]),
                    'raw_intensity': np.asarray(sig.raw_intensity) if hasattr(sig, 'raw_intensity') else np.array([]),
                }

            # Clear batch data
            del batch_signals
            gc.collect()

        # Now accumulate in order of precursor_ids
        print("Accumulating MS1 data in precursor order...")
        for pid in precursor_ids:
            ms1 = all_ms1_lookup.get(pid)

            if ms1 is not None:
                # 1D projections
                xic_rt.extend(ms1['rt_coords'])
                xic_int.extend(ms1['rt_intensities'])
                offsets_xic.append(len(xic_rt))

                mob_im.extend(ms1['im_coords'])
                mob_int.extend(ms1['im_intensities'])
                offsets_mob.append(len(mob_im))

                iso_mz.extend(ms1['isotope_mz'])
                iso_int.extend(ms1['isotope_intensity'])
                offsets_iso.append(len(iso_mz))

                # Raw 4D data
                raw_rt_all.extend(ms1['raw_rt'])
                raw_mz_all.extend(ms1['raw_mz'])
                raw_mobility_all.extend(ms1['raw_mobility'])
                raw_intensity_all.extend(ms1['raw_intensity'])
                raw_offsets.append(len(raw_rt_all))
            else:
                offsets_xic.append(offsets_xic[-1])
                offsets_mob.append(offsets_mob[-1])
                offsets_iso.append(offsets_iso[-1])
                raw_offsets.append(raw_offsets[-1])

        # Clear ms1 lookup
        del all_ms1_lookup
        gc.collect()

        # Write MS1 signals (ragged arrays)
        print("Writing MS1 data to Zarr...")
        ms1_group = root.create_group('ms1')

        ms1_group.create_array('xic_rt', data=np.array(xic_rt, dtype=np.float64), compressors=compressors)
        ms1_group.create_array('xic_intensity', data=np.array(xic_int, dtype=np.float64), compressors=compressors)
        ms1_group.create_array('xic_offsets', data=np.array(offsets_xic, dtype=np.int64), compressors=compressors)
        ms1_group.create_array('mobilogram_im', data=np.array(mob_im, dtype=np.float64), compressors=compressors)
        ms1_group.create_array('mobilogram_intensity', data=np.array(mob_int, dtype=np.float64), compressors=compressors)
        ms1_group.create_array('mobilogram_offsets', data=np.array(offsets_mob, dtype=np.int64), compressors=compressors)
        ms1_group.create_array('isotope_mz', data=np.array(iso_mz, dtype=np.float64), compressors=compressors)
        ms1_group.create_array('isotope_intensity', data=np.array(iso_int, dtype=np.float64), compressors=compressors)
        ms1_group.create_array('isotope_offsets', data=np.array(offsets_iso, dtype=np.int64), compressors=compressors)
        # Backward compat: also save as 'offsets' (XIC)
        ms1_group.create_array('offsets', data=np.array(offsets_xic, dtype=np.int64), compressors=compressors)

        # Write raw 4D MS1 data
        ms1_group.create_array('raw_rt', data=np.array(raw_rt_all, dtype=np.float64), compressors=compressors)
        ms1_group.create_array('raw_mz', data=np.array(raw_mz_all, dtype=np.float64), compressors=compressors)
        ms1_group.create_array('raw_mobility', data=np.array(raw_mobility_all, dtype=np.float64), compressors=compressors)
        ms1_group.create_array('raw_intensity', data=np.array(raw_intensity_all, dtype=np.float64), compressors=compressors)
        ms1_group.create_array('raw_offsets', data=np.array(raw_offsets, dtype=np.int64), compressors=compressors)

        print(f"Store created: {output_path}")
        print(f"  Precursors: {len(file_index)}")
        print(f"  Fragment peaks: {n_fragment_peaks}")
        print(f"  XIC points: {len(xic_rt)}")
        print(f"  Raw 4D points: {len(raw_rt_all)}")

        return cls(str(output_path))


def main():
    """CLI for creating stores."""
    import argparse

    parser = argparse.ArgumentParser(description="Create Zarr precursor store")
    parser.add_argument("--index", required=True, help="Path to precursor_index.parquet")
    parser.add_argument("--raw", required=True, help="Path to .d folder")
    parser.add_argument("--output", required=True, help="Output .zarr path")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--rt-window", type=float, default=30.0)
    parser.add_argument("--batch-size", type=int, default=2000, help="Batch size for MS1 extraction")
    parser.add_argument("--limit", type=int, default=None, help="Limit number of precursors (for testing)")

    args = parser.parse_args()

    PrecursorStore.create_from_index_and_raw(
        index_path=args.index,
        raw_data_path=args.raw,
        output_path=args.output,
        num_threads=args.threads,
        rt_window_sec=args.rt_window,
        batch_size=args.batch_size,
        limit=args.limit,
    )


if __name__ == "__main__":
    main()
