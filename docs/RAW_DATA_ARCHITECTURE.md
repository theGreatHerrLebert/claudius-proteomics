# Raw Data Architecture

This document describes how raw timsTOF 4D mass spectrometry data is extracted, stored, and served for the Precursor Browser dashboard.

## Overview

The pipeline extracts full 4D signal data (RT × IM × m/z × intensity) from timsTOF raw files and stores it in a training-optimized Parquet format with variable-length list columns.

```
Raw .d files → Rust extraction → Parquet store → FastAPI backend → React frontend
```

## Data Extraction

### Source Data

Each precursor in a timsTOF DDA experiment has:
1. **Fragment spectrum**: MS2 peaks from PASEF fragmentation (merged across re-captured frames)
2. **MS1 precursor signal**: Extracted from surrounding MS1 frames in a configurable RT window

### Extraction Process

The extraction is performed by `scripts/precursor_store_parquet.py`:

```python
PrecursorStoreParquet.create_from_index_and_raw(
    index_path="precursor_index.parquet",  # Unified index from search engines
    raw_data_path="data.d",                 # timsTOF raw folder
    output_path="precursors.parquet",
    rt_window_sec=30.0,                     # ±15 sec around precursor RT
    mz_tol_ppm=20.0,                        # m/z tolerance for MS1 extraction
    im_window=0.1,                          # ±0.05 1/K0 around precursor
    batch_size=5000,                        # Process in batches to manage memory
)
```

### Rust Backend (imspy-connector)

The heavy lifting is done in Rust via `rustdf`:

1. **Fragment extraction**: `get_pasef_fragments()` returns merged fragment spectra per precursor
2. **MS1 signal extraction**: `extract_precursor_ms1_signals()` extracts:
   - 1D projections: XIC (RT profile), mobilogram (IM profile), isotope envelope
   - **Raw 4D data**: All individual peaks from MS1 frames in the RT window

The raw 4D data preserves every detected peak with its full coordinates:
- `raw_rt`: Retention time in seconds
- `raw_mz`: m/z value
- `raw_mobility`: Ion mobility (1/K0)
- `raw_intensity`: Signal intensity

## Storage Format

### Parquet with List Columns

We use Apache Parquet with list columns for variable-length data. This provides:
- **Columnar storage**: Only read columns you need
- **Predicate pushdown**: Filter by scalar columns before loading arrays
- **Row group batching**: Natural batching for training (e.g., 10k precursors per group)
- **Compression**: ZSTD compression, ~5-6 KB per precursor

### Schema

```
precursor_id:         int64
raw_file:             string
mz:                   float64
charge:               int32
rt_seconds:           float64
mobility:             float64

# Search engine results
fragpipe_peptide:     string (nullable)
sage_peptide:         string (nullable)
diann_peptide:        string (nullable)
consensus_peptide:    string (nullable)
n_engines:            int32

# Quality metrics
fragpipe_probability: float64 (nullable)
diann_qvalue:         float64 (nullable)
sage_qvalue:          float64 (nullable)
raw_intensity_meta:   float64 (nullable)  # timsTOF metadata intensity

# Fragment spectrum (variable length)
fragment_mz:          list<float64>
fragment_intensity:   list<float64>
fragment_mobility:    list<float64>

# MS1 1D projections (variable length)
xic_rt:               list<float64>
xic_intensity:        list<float64>
mobilogram_im:        list<float64>
mobilogram_intensity: list<float64>
isotope_mz:           list<float64>
isotope_intensity:    list<float64>

# Raw 4D MS1 data (variable length, ~1-10k points per precursor)
raw_rt:               list<float64>
raw_mz:               list<float64>
raw_mobility:         list<float64>
raw_intensity:        list<float64>
```

### Training Spectra Schema (Extended)

The `training_spectra.parquet` file extends the base schema with normalized intensities
and fragment ion annotations. These columns are added by `build_training_spectra.py`.

```
# Normalized intensities (same length as fragment_mz)
fragment_intensity_norm:   list<float64>   # Normalized to 0-100 scale
normalization_factor:      float64         # Factor for reversibility
normalization_method:      string          # "base_peak", "tic", "sqrt", "log"

# Fragment annotations (parallel arrays, same length as fragment_mz)
fragment_ion_type:         list<string>    # "b", "y", or "other" (unmatched)
fragment_ion_number:       list<int32>     # Ion position (1-indexed), 0 for unmatched
fragment_ion_charge:       list<int32>     # Fragment charge state
fragment_theoretical_mz:   list<float64>   # Theoretical m/z (0 for unmatched)
fragment_error_ppm:        list<float64>   # Mass error in ppm (0 for unmatched)

# Annotation summary metrics
is_annotated:              bool            # Has sequence for annotation
n_matched_peaks:           int32           # Number of peaks matched to theoretical
intensity_explained:       float64         # Fraction of total intensity matched (0-1)
sequence_coverage_b:       float64         # Fraction of b ions found (0-1)
sequence_coverage_y:       float64         # Fraction of y ions found (0-1)
```

#### Normalization Methods

| Method | Formula | Use Case |
|--------|---------|----------|
| `base_peak` | I / max(I) × 100 | Standard, preserves relative intensities |
| `tic` | I / sum(I) × 100 | Accounts for spectrum complexity |
| `sqrt` | sqrt(I / max(I)) × 100 | Compresses dynamic range |
| `log` | log10(I+1) / log10(max+1) × 100 | Very wide dynamic range |

#### Building Training Spectra

```bash
# Via Snakemake
snakemake build_training_spectra --config accession=PXD019086

# Directly
python scripts/build_training_spectra.py \
    --input data/processed/PXD019086/precursors.parquet \
    --output data/processed/PXD019086/training_spectra.parquet \
    --normalization base_peak \
    --mz-tolerance 20.0
```

#### Example: Loading Annotated Spectra

```python
import pyarrow.parquet as pq
import numpy as np

# Read training spectra
table = pq.read_table(
    "training_spectra.parquet",
    columns=[
        'precursor_id', 'sage_peptide', 'charge',
        'fragment_mz', 'fragment_intensity_norm',
        'fragment_ion_type', 'fragment_ion_number',
        'intensity_explained', 'sequence_coverage_b'
    ],
    filters=[('is_annotated', '=', True), ('n_engines', '>=', 2)]
)

df = table.to_pandas()
print(f"Annotated spectra: {len(df)}")
print(f"Mean intensity explained: {df['intensity_explained'].mean():.1%}")
print(f"Mean b-ion coverage: {df['sequence_coverage_b'].mean():.1%}")

# Access parallel arrays for a spectrum
row = df.iloc[0]
for i, (mz, ion_type, ion_num) in enumerate(zip(
    row['fragment_mz'],
    row['fragment_ion_type'],
    row['fragment_ion_number']
)):
    if ion_type != 'other':
        print(f"  m/z {mz:.4f}: {ion_type}{ion_num}")
```

### Typical Sizes

| Precursors | Fragment Peaks | Raw 4D Points | File Size |
|------------|----------------|---------------|-----------|
| 1,000      | 197k           | 400k          | 2.8 MB    |
| 10,000     | 4.1M           | 8.5M          | 53 MB     |
| 214,824    | 218M           | 1.4M*         | 1.3 GB    |

*Only precursors with valid coordinates get raw 4D extraction

## Training Data Access

### Batch Iteration

```python
import pyarrow.parquet as pq

pf = pq.ParquetFile("precursors.parquet")

# Iterate in batches (memory efficient)
for batch in pf.iter_batches(
    batch_size=1000,
    columns=['precursor_id', 'mz', 'charge', 'raw_rt', 'raw_mz', 'raw_mobility', 'raw_intensity']
):
    # batch is a PyArrow RecordBatch
    for i in range(len(batch)):
        raw_rt = batch.column('raw_rt')[i].as_py()  # List[float]
        # Process...
```

### Filtered Reads

```python
# Only load charge 2+ precursors with 2+ search engine hits
table = pq.read_table(
    "precursors.parquet",
    columns=['precursor_id', 'raw_rt', 'raw_mz', 'raw_mobility', 'raw_intensity'],
    filters=[('charge', '=', 2), ('n_engines', '>=', 2)]
)
```

### PyTorch DataLoader Integration

```python
class PrecursorDataset(torch.utils.data.Dataset):
    def __init__(self, parquet_path, columns):
        self.table = pq.read_table(parquet_path, columns=columns)

    def __len__(self):
        return len(self.table)

    def __getitem__(self, idx):
        row = self.table.slice(idx, 1).to_pandas().iloc[0]
        return {
            'raw_rt': torch.tensor(row['raw_rt']),
            'raw_mz': torch.tensor(row['raw_mz']),
            'raw_mobility': torch.tensor(row['raw_mobility']),
            'raw_intensity': torch.tensor(row['raw_intensity']),
        }
```

## Dashboard Architecture

### Backend (FastAPI)

Location: `dashboard/backend/main.py`

The backend serves data from the Parquet store via REST API:

```
GET  /                      # Health check
POST /load?path=...         # Load a Parquet store
GET  /info                  # Store metadata
GET  /stats                 # Summary statistics
GET  /precursors            # List with pagination, filtering, sorting
GET  /precursor/{id}        # Full detail for single precursor
```

#### Key Endpoints

**List precursors** (scalar columns only, fast):
```
GET /precursors?offset=0&limit=100&min_engines=2&charge=2&sort_by=raw_intensity_meta&sort_desc=true
```

**Get precursor detail** (includes all array data):
```
GET /precursor/12345
```

Returns full fragment spectrum, MS1 projections, and raw 4D data.

### Frontend (React + TypeScript)

Location: `dashboard/frontend/`

Tech stack:
- **Vite**: Build tool
- **React**: UI framework
- **TypeScript**: Type safety
- **Tailwind CSS**: Styling
- **TanStack Query**: Data fetching with caching
- **TanStack Table**: Virtualized table
- **Deck.gl**: WebGL scatter plots for large point clouds

#### Components

```
App.tsx
├── Header (filters, stats)
├── PrecursorTable (paginated, sortable, filterable)
└── PrecursorViz
    ├── Precursor info header
    ├── Fragment spectrum (SVG, log scale)
    ├── IM vs m/z scatter (Deck.gl WebGL)
    ├── Raw 4D: IM vs RT scatter (Deck.gl WebGL)
    ├── XIC profile (SVG)
    ├── Mobilogram profile (SVG)
    └── Isotope envelope (SVG)
```

#### Data Flow

```
User clicks row → React Query fetches /precursor/{id} →
  Response cached → PrecursorViz renders WebGL plots
```

### Running the Dashboard

```bash
# Terminal 1: Backend
cd dashboard/backend
python main.py --store ../data/processed/PXD019086/frac01_full.parquet --port 8000

# Terminal 2: Frontend
cd dashboard/frontend
npm run dev  # http://localhost:5173
```

Or use the run script:
```bash
./dashboard/run.sh data/processed/PXD019086/frac01_full.parquet
```

## Design Decisions

### Why Parquet over Zarr?

| Aspect | Parquet | Zarr |
|--------|---------|------|
| Batch reads | Excellent (row groups) | Poor (random access) |
| Column selection | Native | Manual |
| Filtering | Predicate pushdown | Post-load |
| Variable-length | List columns | Ragged arrays + offsets |
| Ecosystem | Arrow, DuckDB, Polars | NumPy-centric |

For training workloads with batch iteration and column selection, Parquet is more efficient.

### Why preserve raw 4D data?

1. **No information loss**: Binning/projection loses signal structure
2. **Model flexibility**: Let the model learn representations
3. **QC/debugging**: Visual inspection of extraction quality
4. **Future-proof**: Can derive any projection from raw data

### Memory Management

The extraction processes precursors in batches to avoid OOM:
1. Fragment extraction loads all at once (unavoidable), but stores minimal data
2. MS1 extraction processes `batch_size` precursors at a time
3. ParquetWriter streams to disk incrementally
4. Each batch is garbage collected before the next

## Future Improvements

1. **Streaming extraction**: Process raw files without loading all fragments
2. **Delta encoding**: Compress RT/IM coordinates (monotonic within precursor)
3. **Quantization**: Store intensities as float16 or log-scaled int16
4. **Sharding**: Split large datasets across multiple Parquet files
5. **Arrow IPC**: Direct memory-mapped access for training
