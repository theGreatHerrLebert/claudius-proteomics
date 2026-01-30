#!/usr/bin/env python3
"""
Build Training Spectra Pipeline

Reads precursor data from existing Parquet stores, applies normalization,
annotates fragments with b/y ions, and writes extended Parquet with
training-ready columns.

This script combines:
1. Raw spectrum data from precursor_store_parquet.py
2. Normalization from spectrum_normalization.py
3. Annotation from annotate_spectrum.py (using fragment_matching.py)

Output adds these columns to the input Parquet:
- fragment_intensity_norm: Normalized intensities (list<float64>)
- normalization_factor: For reversibility (float64)
- fragment_ion_type: "b", "y", "other" per peak (list<string>)
- fragment_ion_number: Ion position, 0 for unmatched (list<int32>)
- fragment_ion_charge: Fragment charge (list<int32>)
- fragment_theoretical_mz: Theoretical m/z (list<float64>)
- fragment_error_ppm: Mass errors (list<float64>)
- is_annotated: Has sequence for annotation (bool)
- n_matched_peaks: Peaks matched to theoretical (int32)
- intensity_explained: Fraction of intensity matched (float64)
- sequence_coverage_b: b ion coverage (float64)
- sequence_coverage_y: y ion coverage (float64)

Usage:
    python build_training_spectra.py \\
        --input data/processed/PXD019086/precursor_index_v3.parquet \\
        --output /tmp/training_spectra.parquet \\
        --normalization base_peak \\
        --mz-tolerance 20.0 \\
        --limit 100

    # Verify output
    python -c "import pyarrow.parquet as pq; print(pq.read_schema('output.parquet'))"
"""

import sys
import gc
from pathlib import Path
from typing import Optional, List
import argparse

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# Add scripts directory to path
SCRIPTS_DIR = Path(__file__).parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from spectrum_normalization import (
    normalize_spectrum,
    NormalizationMethod,
)
from annotate_spectrum import (
    SpectrumAnnotator,
    IMSPY_AVAILABLE,
)


def get_sequence_for_annotation(row: pd.Series) -> Optional[str]:
    """
    Get the best sequence for annotation from a precursor row.

    Priority: sage_modified > sage_peptide > consensus_peptide

    Args:
        row: DataFrame row with peptide columns

    Returns:
        Sequence string or None if no sequence available
    """
    # Try modified sequences first (contain modification info)
    for col in ['sage_modified', 'sage_peptide', 'consensus_peptide',
                'fragpipe_modified', 'diann_modified']:
        val = row.get(col)
        if pd.notna(val) and isinstance(val, str) and len(val) > 0:
            return val

    return None


def build_training_spectra(
    input_path: str,
    output_path: str,
    normalization: str = "base_peak",
    mz_tolerance_ppm: float = 20.0,
    min_intensity: float = 0.0,
    batch_size: int = 1000,
    limit: Optional[int] = None,
    sequence_col: Optional[str] = None,
) -> dict:
    """
    Build training spectra with normalization and annotation.

    Args:
        input_path: Path to input Parquet (precursor store or index)
        output_path: Path for output Parquet
        normalization: Normalization method (base_peak, tic, sqrt, log)
        mz_tolerance_ppm: m/z tolerance for fragment matching
        min_intensity: Minimum intensity threshold
        batch_size: Batch size for processing
        limit: Limit number of precursors (for testing)
        sequence_col: Override sequence column selection

    Returns:
        Dictionary with processing statistics
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Building training spectra")
    print(f"=" * 60)
    print(f"Input:         {input_path}")
    print(f"Output:        {output_path}")
    print(f"Normalization: {normalization}")
    print(f"m/z tolerance: {mz_tolerance_ppm} ppm")

    # Parse normalization method
    try:
        norm_method = NormalizationMethod.from_string(normalization)
    except ValueError as e:
        print(f"ERROR: {e}")
        return {'error': str(e)}

    # Check imspy availability
    if not IMSPY_AVAILABLE:
        print("WARNING: imspy-core not available. Annotation will be skipped.")

    # Initialize annotator if available
    annotator = None
    if IMSPY_AVAILABLE:
        annotator = SpectrumAnnotator(
            mz_tolerance_ppm=mz_tolerance_ppm,
            min_intensity=min_intensity,
        )

    # Read input
    print(f"\nLoading input...")
    pf = pq.ParquetFile(str(input_path))
    schema = pf.schema_arrow

    # Check for required columns
    available_cols = set(schema.names)
    required_cols = {'fragment_mz', 'fragment_intensity'}

    if not required_cols.issubset(available_cols):
        missing = required_cols - available_cols
        print(f"ERROR: Missing required columns: {missing}")
        print(f"Available columns: {sorted(available_cols)}")
        return {'error': f'Missing columns: {missing}'}

    # Get total count
    total_rows = pf.metadata.num_rows
    if limit:
        total_rows = min(total_rows, limit)
    print(f"Processing {total_rows:,} precursors")

    # Define output schema (input schema + new columns)
    new_fields = [
        ('fragment_intensity_norm', pa.list_(pa.float64())),
        ('normalization_factor', pa.float64()),
        ('normalization_method', pa.string()),
        ('fragment_ion_type', pa.list_(pa.string())),
        ('fragment_ion_number', pa.list_(pa.int32())),
        ('fragment_ion_charge', pa.list_(pa.int32())),
        ('fragment_theoretical_mz', pa.list_(pa.float64())),
        ('fragment_error_ppm', pa.list_(pa.float64())),
        ('is_annotated', pa.bool_()),
        ('n_matched_peaks', pa.int32()),
        ('intensity_explained', pa.float64()),
        ('sequence_coverage_b', pa.float64()),
        ('sequence_coverage_y', pa.float64()),
    ]

    # Build output schema: original + new fields
    output_fields = list(schema)
    for name, pa_type in new_fields:
        if name not in available_cols:
            output_fields.append(pa.field(name, pa_type))
    output_schema = pa.schema(output_fields)

    # Stats
    stats = {
        'n_precursors': 0,
        'n_with_spectrum': 0,
        'n_annotated': 0,
        'n_matched_total': 0,
        'intensity_explained_sum': 0.0,
        'coverage_b_sum': 0.0,
        'coverage_y_sum': 0.0,
    }

    # Process in batches
    writer = None
    rows_processed = 0

    print(f"\nProcessing batches...")

    for batch in pf.iter_batches(batch_size=batch_size):
        if limit and rows_processed >= limit:
            break

        # Convert to pandas for processing
        df = batch.to_pandas()

        if limit:
            remaining = limit - rows_processed
            if len(df) > remaining:
                df = df.head(remaining)

        # Process each row
        new_data = {
            'fragment_intensity_norm': [],
            'normalization_factor': [],
            'normalization_method': [],
            'fragment_ion_type': [],
            'fragment_ion_number': [],
            'fragment_ion_charge': [],
            'fragment_theoretical_mz': [],
            'fragment_error_ppm': [],
            'is_annotated': [],
            'n_matched_peaks': [],
            'intensity_explained': [],
            'sequence_coverage_b': [],
            'sequence_coverage_y': [],
        }

        for _, row in df.iterrows():
            stats['n_precursors'] += 1

            # Get spectrum data
            frag_mz = row.get('fragment_mz')
            frag_int = row.get('fragment_intensity')

            # Handle missing/empty spectrum
            if frag_mz is None or not isinstance(frag_mz, (list, np.ndarray)) or len(frag_mz) == 0:
                # Empty spectrum - add placeholder values
                new_data['fragment_intensity_norm'].append([])
                new_data['normalization_factor'].append(0.0)
                new_data['normalization_method'].append(norm_method.value)
                new_data['fragment_ion_type'].append([])
                new_data['fragment_ion_number'].append([])
                new_data['fragment_ion_charge'].append([])
                new_data['fragment_theoretical_mz'].append([])
                new_data['fragment_error_ppm'].append([])
                new_data['is_annotated'].append(False)
                new_data['n_matched_peaks'].append(0)
                new_data['intensity_explained'].append(0.0)
                new_data['sequence_coverage_b'].append(0.0)
                new_data['sequence_coverage_y'].append(0.0)
                continue

            stats['n_with_spectrum'] += 1

            # Convert to numpy arrays
            frag_mz = np.array(frag_mz, dtype=np.float64)
            frag_int = np.array(frag_int, dtype=np.float64)

            # Normalize
            norm_result = normalize_spectrum(frag_mz, frag_int, norm_method, min_intensity)
            new_data['fragment_intensity_norm'].append(norm_result.normalized_intensity.tolist())
            new_data['normalization_factor'].append(norm_result.normalization_factor)
            new_data['normalization_method'].append(norm_method.value)

            # Annotate
            if annotator:
                sequence = sequence_col and row.get(sequence_col) or get_sequence_for_annotation(row)
                charge = int(row.get('charge', 2))

                annotation = annotator.annotate(
                    sequence=sequence,
                    charge=charge,
                    fragment_mz=frag_mz,
                    fragment_intensity=frag_int,
                )

                new_data['fragment_ion_type'].append(annotation.ion_type)
                new_data['fragment_ion_number'].append(annotation.ion_number.tolist())
                new_data['fragment_ion_charge'].append(annotation.ion_charge.tolist())
                new_data['fragment_theoretical_mz'].append(annotation.theoretical_mz.tolist())
                new_data['fragment_error_ppm'].append(annotation.error_ppm.tolist())
                new_data['is_annotated'].append(annotation.is_annotated)
                new_data['n_matched_peaks'].append(annotation.n_matched_peaks)
                new_data['intensity_explained'].append(annotation.intensity_explained)
                new_data['sequence_coverage_b'].append(annotation.sequence_coverage_b)
                new_data['sequence_coverage_y'].append(annotation.sequence_coverage_y)

                if annotation.is_annotated:
                    stats['n_annotated'] += 1
                    stats['n_matched_total'] += annotation.n_matched_peaks
                    stats['intensity_explained_sum'] += annotation.intensity_explained
                    stats['coverage_b_sum'] += annotation.sequence_coverage_b
                    stats['coverage_y_sum'] += annotation.sequence_coverage_y
            else:
                # No annotator - add empty annotations
                n_peaks = len(frag_mz)
                new_data['fragment_ion_type'].append(['other'] * n_peaks)
                new_data['fragment_ion_number'].append([0] * n_peaks)
                new_data['fragment_ion_charge'].append([0] * n_peaks)
                new_data['fragment_theoretical_mz'].append([0.0] * n_peaks)
                new_data['fragment_error_ppm'].append([0.0] * n_peaks)
                new_data['is_annotated'].append(False)
                new_data['n_matched_peaks'].append(0)
                new_data['intensity_explained'].append(0.0)
                new_data['sequence_coverage_b'].append(0.0)
                new_data['sequence_coverage_y'].append(0.0)

        # Add new columns to dataframe
        for col, values in new_data.items():
            df[col] = values

        # Convert to PyArrow table
        output_table = pa.Table.from_pandas(df, schema=output_schema, preserve_index=False)

        # Write batch
        if writer is None:
            writer = pq.ParquetWriter(
                str(output_path),
                output_schema,
                compression='zstd',
                compression_level=3,
            )
        writer.write_table(output_table)

        rows_processed += len(df)
        print(f"  Processed {rows_processed:,} / {total_rows:,} "
              f"({100 * rows_processed / total_rows:.1f}%)")

        # Clean up
        del df, output_table
        gc.collect()

    # Close writer
    if writer:
        writer.close()

    # Compute summary stats
    if stats['n_annotated'] > 0:
        stats['mean_intensity_explained'] = stats['intensity_explained_sum'] / stats['n_annotated']
        stats['mean_coverage_b'] = stats['coverage_b_sum'] / stats['n_annotated']
        stats['mean_coverage_y'] = stats['coverage_y_sum'] / stats['n_annotated']
    else:
        stats['mean_intensity_explained'] = 0.0
        stats['mean_coverage_b'] = 0.0
        stats['mean_coverage_y'] = 0.0

    # Print summary
    print(f"\n" + "=" * 60)
    print(f"Processing complete")
    print(f"=" * 60)
    print(f"Total precursors:      {stats['n_precursors']:,}")
    print(f"With spectrum:         {stats['n_with_spectrum']:,}")
    print(f"Annotated:             {stats['n_annotated']:,}")
    print(f"Total matched peaks:   {stats['n_matched_total']:,}")
    print(f"Mean intensity explained: {stats['mean_intensity_explained']:.1%}")
    print(f"Mean b coverage:       {stats['mean_coverage_b']:.1%}")
    print(f"Mean y coverage:       {stats['mean_coverage_y']:.1%}")
    print(f"\nOutput: {output_path}")
    print(f"Size:   {output_path.stat().st_size / 1e6:.1f} MB")

    return stats


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Build training spectra with normalization and annotation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Process with default settings
    python build_training_spectra.py \\
        --input data/processed/PXD019086/precursors.parquet \\
        --output training_spectra.parquet

    # Test with limit
    python build_training_spectra.py \\
        --input precursor_index_v3.parquet \\
        --output /tmp/test.parquet \\
        --limit 100

    # Different normalization
    python build_training_spectra.py \\
        --input precursors.parquet \\
        --output training.parquet \\
        --normalization sqrt \\
        --mz-tolerance 10.0
        """,
    )

    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input Parquet file (precursor store or index)",
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output Parquet file",
    )
    parser.add_argument(
        "--normalization", "-n",
        default="base_peak",
        choices=["base_peak", "tic", "sqrt", "log"],
        help="Normalization method (default: base_peak)",
    )
    parser.add_argument(
        "--mz-tolerance",
        type=float,
        default=20.0,
        help="m/z tolerance in ppm for fragment matching (default: 20.0)",
    )
    parser.add_argument(
        "--min-intensity",
        type=float,
        default=0.0,
        help="Minimum intensity threshold (default: 0.0)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=1000,
        help="Batch size for processing (default: 1000)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit number of precursors (for testing)",
    )
    parser.add_argument(
        "--sequence-col",
        default=None,
        help="Override sequence column for annotation",
    )

    args = parser.parse_args()

    stats = build_training_spectra(
        input_path=args.input,
        output_path=args.output,
        normalization=args.normalization,
        mz_tolerance_ppm=args.mz_tolerance,
        min_intensity=args.min_intensity,
        batch_size=args.batch_size,
        limit=args.limit,
        sequence_col=args.sequence_col,
    )

    if 'error' in stats:
        sys.exit(1)


if __name__ == "__main__":
    main()
