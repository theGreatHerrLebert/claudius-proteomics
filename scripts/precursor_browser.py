#!/usr/bin/env python3
"""
Precursor Browser GUI

Interactive browser for visualizing precursors with search engine identifications.

Features:
- Table view of all targeted precursors
- Sorted by agreement (n_engines) first
- Filter by engine, peptide, m/z range
- Click to visualize: fragment spectrum + MS1 signal

Usage:
    # With Zarr store (fast, recommended)
    python precursor_browser.py --store /path/to/precursors.zarr

    # With index + raw data (extracts on-the-fly)
    python precursor_browser.py --index /path/to/index.parquet --raw /path/to/data.d

    # Serve as web app
    panel serve precursor_browser.py --args --store /path/to/precursors.zarr
"""

import sys
import argparse
from pathlib import Path
from typing import Optional
import numpy as np
import pandas as pd

import panel as pn
import param
import holoviews as hv
from holoviews import opts

# Enable Panel extensions
pn.extension('tabulator')
hv.extension('bokeh')

# Try to import precursor stores (Zarr and Parquet)
try:
    from precursor_store import PrecursorStore, PrecursorData
    from precursor_store_parquet import PrecursorStoreParquet
except ImportError:
    sys.path.insert(0, str(Path(__file__).parent))
    from precursor_store import PrecursorStore, PrecursorData
    from precursor_store_parquet import PrecursorStoreParquet


def load_store(store_path: str):
    """Load either Zarr or Parquet store based on file extension."""
    path = Path(store_path)
    if path.suffix == '.parquet':
        return PrecursorStoreParquet(str(path))
    elif path.suffix == '.zarr' or path.is_dir():
        return PrecursorStore(str(path))
    else:
        raise ValueError(f"Unknown store format: {path.suffix}")


class PrecursorBrowser(param.Parameterized):
    """Interactive precursor browser with Panel GUI."""

    # Parameters
    selected_precursor_id = param.Integer(default=None, allow_None=True)
    min_engines = param.Integer(default=1, bounds=(0, 3))
    search_peptide = param.String(default="")
    search_protein = param.String(default="")

    def __init__(self, store: PrecursorStore, **kwargs):
        super().__init__(**kwargs)
        self.store = store
        self._setup_data()
        self._setup_widgets()
        self._setup_layout()

    def _setup_data(self):
        """Prepare data for display."""
        # Get sorted index
        self.full_df = self.store.get_by_agreement(min_engines=0)

        # Add signal quality column for sorting
        # Prefer raw_intensity from timsTOF metadata, fall back to fragpipe_probability
        if 'raw_intensity' in self.full_df.columns:
            self.full_df['signal_quality'] = self.full_df['raw_intensity'].fillna(0)
        elif 'fragpipe_probability' in self.full_df.columns:
            self.full_df['signal_quality'] = self.full_df['fragpipe_probability'].fillna(0) * 1e6
        else:
            self.full_df['signal_quality'] = 1.0

        # Add display columns
        self.full_df['agreement'] = self.full_df['n_engines'].map({
            3: '3',
            2: '2',
            1: '1',
            0: '0'
        })

        # Determine match status
        def get_match_status(row):
            fp = str(row.get('fragpipe_peptide', '')).replace('I', 'L').upper()
            sage = str(row.get('sage_peptide', '')).replace('I', 'L').upper()
            diann = str(row.get('diann_peptide', '')).replace('I', 'L').upper()

            peptides = [p for p in [fp, sage, diann] if p and p != 'NAN' and p != 'NONE']

            if len(peptides) <= 1:
                return '-'
            if len(set(peptides)) == 1:
                return 'Y'
            return 'N'

        self.full_df['match'] = self.full_df.apply(get_match_status, axis=1)

        # Select display columns
        self.display_cols = [
            'precursor_id', 'agreement', 'match', 'raw_intensity',
            'consensus_peptide', 'fragpipe_peptide', 'sage_peptide', 'diann_peptide',
            'fragpipe_probability', 'sage_qvalue', 'diann_qvalue'
        ]
        self.display_cols = [c for c in self.display_cols if c in self.full_df.columns]

    def _setup_widgets(self):
        """Create Panel widgets."""
        # Filters - horizontal layout
        self.identified_filter = pn.widgets.Checkbox(
            name='Identified only',
            value=True,  # Default to showing only identified
            width=120
        )

        self.engine_filter = pn.widgets.Select(
            name='Min Engines',
            options={'All': 0, '1+': 1, '2+': 2, '3': 3},
            value=0,
            width=100
        )

        self.agreement_filter = pn.widgets.Select(
            name='Agreement',
            options={'All': 'all', 'Agree': 'agree', 'Disagree': 'disagree'},
            value='all',
            width=100
        )

        self.peptide_search = pn.widgets.TextInput(
            name='Search Peptide',
            placeholder='Peptide sequence...',
            width=200
        )

        # Sort options
        self.sort_by = pn.widgets.Select(
            name='Sort',
            options={'Intensity (high)': 'intensity', 'Engines (desc)': 'engines', 'Precursor ID': 'pid'},
            value='intensity',
            width=130
        )

        # Table - full width, more rows, NOT editable
        self.table = pn.widgets.Tabulator(
            self._get_filtered_df(),
            selectable='toggle',
            selection=[],
            pagination='remote',
            page_size=30,
            sizing_mode='stretch_both',
            show_index=False,
            disabled=True,  # Not editable
            formatters={
                'fragpipe_probability': {'type': 'progress', 'max': 1.0, 'color': '#3498db'},
            },
            configuration={
                'columnDefaults': {
                    'headerFilter': True,
                },
                'selectable': 1,
            }
        )

        # Link callbacks
        self.identified_filter.param.watch(self._on_filter_change, 'value')
        self.engine_filter.param.watch(self._on_filter_change, 'value')
        self.agreement_filter.param.watch(self._on_filter_change, 'value')
        self.peptide_search.param.watch(self._on_filter_change, 'value')
        self.sort_by.param.watch(self._on_filter_change, 'value')
        self.table.param.watch(self._on_selection, 'selection')

        # Plots (initially empty)
        self.fragment_plot = pn.pane.HoloViews(hv.Curve([]), sizing_mode='stretch_width', height=250)
        self.im_mz_plot = pn.pane.HoloViews(hv.Points([]), sizing_mode='stretch_width', height=200)
        self.raw_im_rt_plot = pn.pane.HoloViews(hv.Points([]), sizing_mode='stretch_width', height=200)
        self.xic_plot = pn.pane.HoloViews(hv.Curve([]), sizing_mode='stretch_width', height=150)
        self.mobilogram_plot = pn.pane.HoloViews(hv.Curve([]), sizing_mode='stretch_width', height=150)
        self.isotope_plot = pn.pane.HoloViews(hv.Bars([]), sizing_mode='stretch_width', height=150)

        # Info panel
        self.info_panel = pn.pane.Markdown("*Select a precursor*", width=300)

    def _setup_layout(self):
        """Create panel layout - filters on top, table 2/3, viz 1/3."""
        # Filter bar on top
        self.filter_bar = pn.Row(
            pn.pane.Markdown(f"**{len(self.full_df):,} precursors**", width=150),
            self.identified_filter,
            self.engine_filter,
            self.agreement_filter,
            self.sort_by,
            self.peptide_search,
            sizing_mode='stretch_width',
            height=60,
        )

        # Visualization panel (right side, 1/3)
        self.viz_panel = pn.Column(
            self.info_panel,
            pn.layout.Divider(),
            pn.pane.Markdown("**Fragment Spectrum (log)**"),
            self.fragment_plot,
            pn.pane.Markdown("**IM vs m/z (fragments)**"),
            self.im_mz_plot,
            pn.pane.Markdown("**Raw 4D: IM vs RT (precursor)**"),
            self.raw_im_rt_plot,
            pn.Row(
                pn.Column(pn.pane.Markdown("**XIC**"), self.xic_plot),
                pn.Column(pn.pane.Markdown("**Mobilogram**"), self.mobilogram_plot),
            ),
            pn.pane.Markdown("**Isotopes**"),
            self.isotope_plot,
            sizing_mode='stretch_both',
            width=400,
        )

    def _get_filtered_df(self) -> pd.DataFrame:
        """Get filtered DataFrame based on current filter settings."""
        df = self.full_df.copy()

        # Filter identified only
        if hasattr(self, 'identified_filter') and self.identified_filter.value:
            df = df[df['n_engines'] > 0]

        # Filter by min engines
        min_eng = self.engine_filter.value if hasattr(self, 'engine_filter') else 0
        df = df[df['n_engines'] >= min_eng]

        # Filter by agreement
        if hasattr(self, 'agreement_filter'):
            if self.agreement_filter.value == 'agree':
                df = df[df['match'] == 'Y']
            elif self.agreement_filter.value == 'disagree':
                df = df[df['match'] == 'N']

        # Filter by peptide search
        if hasattr(self, 'peptide_search') and self.peptide_search.value:
            search = self.peptide_search.value.upper()
            mask = (
                df['fragpipe_peptide'].str.upper().str.contains(search, na=False) |
                df['sage_peptide'].str.upper().str.contains(search, na=False) |
                df['diann_peptide'].str.upper().str.contains(search, na=False) |
                df['consensus_peptide'].str.upper().str.contains(search, na=False)
            )
            df = df[mask]

        # Sort by selected option
        if hasattr(self, 'sort_by'):
            if self.sort_by.value == 'intensity':
                df = df.sort_values('signal_quality', ascending=False)
            elif self.sort_by.value == 'engines':
                df = df.sort_values(['n_engines', 'signal_quality'], ascending=[False, False])
            elif self.sort_by.value == 'pid':
                df = df.sort_values('precursor_id', ascending=True)

        return df[self.display_cols]

    def _on_filter_change(self, event):
        """Handle filter changes."""
        self.table.value = self._get_filtered_df()

    def _on_selection(self, event):
        """Handle table row selection."""
        if not event.new:
            return

        row_idx = event.new[0]
        filtered_df = self._get_filtered_df()

        if row_idx >= len(filtered_df):
            return

        pid = int(filtered_df.iloc[row_idx]['precursor_id'])
        self._update_visualization(pid)

    def _update_visualization(self, precursor_id: int):
        """Update visualization for selected precursor."""
        data = self.store.get_precursor(precursor_id)
        if data is None:
            self.info_panel.object = f"*Precursor {precursor_id} not found*"
            return

        # Update info panel - compact format
        fp_prob = f"{data.fragpipe_probability:.3f}" if data.fragpipe_probability else "N/A"
        sage_q = f"{data.sage_qvalue:.2e}" if data.sage_qvalue else "N/A"
        diann_q = f"{data.diann_qvalue:.2e}" if data.diann_qvalue else "N/A"

        info_md = f"""**Precursor {data.precursor_id}**

m/z: {data.mz:.4f} | z: {data.charge}+ | RT: {data.rt_seconds/60:.2f}min | IM: {data.mobility:.3f}

**IDs ({data.n_engines} engines):**
- FP: `{data.fragpipe_peptide or '-'}` ({fp_prob})
- Sage: `{data.sage_peptide or '-'}` ({sage_q})
- DIA-NN: `{data.diann_peptide or '-'}` ({diann_q})
"""
        self.info_panel.object = info_md

        # Update fragment spectrum with LOG SCALE
        if data.fragment_mz is not None and len(data.fragment_mz) > 0:
            # Log transform: log10(intensity + 1)
            log_intensity = np.log10(data.fragment_intensity + 1)
            frag_df = pd.DataFrame({
                'mz': data.fragment_mz,
                'intensity': log_intensity
            })
            frag_plot = hv.Spikes(frag_df, 'mz', 'intensity').opts(
                opts.Spikes(color='#3498db', line_width=1, tools=['hover'],
                           xlabel='m/z', ylabel='log10(intensity+1)',
                           title=f'{len(data.fragment_mz)} peaks')
            )
            self.fragment_plot.object = frag_plot
        else:
            self.fragment_plot.object = hv.Text(0, 0, 'No fragment data')

        # Update 2D IM vs m/z plot
        if (data.fragment_mz is not None and len(data.fragment_mz) > 0 and
            data.fragment_mobility is not None and len(data.fragment_mobility) > 0):
            # Log intensity for point size/color
            log_int = np.log10(data.fragment_intensity + 1)
            im_mz_df = pd.DataFrame({
                'mz': data.fragment_mz,
                'im': data.fragment_mobility,
                'log_int': log_int
            })
            im_mz_plot = hv.Points(im_mz_df, ['mz', 'im'], ['log_int']).opts(
                opts.Points(color='log_int', cmap='viridis', size=3,
                           tools=['hover'], colorbar=True,
                           xlabel='m/z', ylabel='Ion Mobility (1/K0)',
                           title='IM vs m/z')
            )
            self.im_mz_plot.object = im_mz_plot
        else:
            self.im_mz_plot.object = hv.Text(0, 0, 'No IM data')

        # Update raw 4D IM vs RT plot (precursor signal from MS1 frames)
        if data.raw_rt is not None and len(data.raw_rt) > 0:
            # Log intensity for color
            log_int = np.log10(np.array(data.raw_intensity) + 1)
            raw_df = pd.DataFrame({
                'rt_min': np.array(data.raw_rt) / 60.0,  # Convert to minutes
                'im': np.array(data.raw_mobility),
                'log_int': log_int
            })
            # Create scatter plot with heatmap coloring
            raw_im_rt = hv.Points(raw_df, ['rt_min', 'im'], ['log_int']).opts(
                opts.Points(color='log_int', cmap='plasma', size=3,
                           tools=['hover'], colorbar=True,
                           xlabel='RT (min)', ylabel='Ion Mobility (1/K0)',
                           title=f'{len(data.raw_rt):,} raw 4D points')
            )
            self.raw_im_rt_plot.object = raw_im_rt
        else:
            self.raw_im_rt_plot.object = hv.Text(0, 0, 'No raw 4D data')

        # Update XIC
        if data.xic_rt is not None and len(data.xic_rt) > 0:
            xic_df = pd.DataFrame({
                'rt_min': np.array(data.xic_rt) / 60.0,
                'intensity': data.xic_intensity
            })
            xic_plot = hv.Area(xic_df, 'rt_min', 'intensity').opts(
                opts.Area(color='#2ecc71', alpha=0.5, line_width=2,
                         xlabel='RT (min)', ylabel='Int')
            )
            self.xic_plot.object = xic_plot
        else:
            self.xic_plot.object = hv.Text(0, 0, 'No XIC')

        # Update mobilogram
        if data.mobilogram_im is not None and len(data.mobilogram_im) > 0:
            mob_df = pd.DataFrame({
                'mobility': data.mobilogram_im,
                'intensity': data.mobilogram_intensity
            })
            mob_plot = hv.Area(mob_df, 'mobility', 'intensity').opts(
                opts.Area(color='#9b59b6', alpha=0.5, line_width=2,
                         xlabel='1/K0', ylabel='Int')
            )
            self.mobilogram_plot.object = mob_plot
        else:
            self.mobilogram_plot.object = hv.Text(0, 0, 'No IM')

        # Update isotope envelope
        if data.isotope_mz is not None and len(data.isotope_mz) > 0:
            # Filter to relevant m/z range: mono_mz - 1 to mono_mz + 4
            mono_mz = data.mz
            mask = (data.isotope_mz >= mono_mz - 1.5) & (data.isotope_mz <= mono_mz + 5)
            if mask.any():
                filtered_mz = data.isotope_mz[mask]
                filtered_int = data.isotope_intensity[mask]
                # Further filter non-zero
                nonzero = filtered_int > 0
                if nonzero.any():
                    iso_df = pd.DataFrame({
                        'mz': filtered_mz[nonzero],
                        'intensity': filtered_int[nonzero]
                    })
                    iso_plot = hv.Spikes(iso_df, 'mz', 'intensity').opts(
                        opts.Spikes(color='#e74c3c', line_width=3,
                                   xlabel='m/z', ylabel='Int',
                                   title=f'Isotopes (mono={mono_mz:.2f})')
                    )
                    self.isotope_plot.object = iso_plot
                else:
                    self.isotope_plot.object = hv.Text(0, 0, 'No isotope signal')
            else:
                self.isotope_plot.object = hv.Text(0, 0, 'No isotopes in range')
        else:
            self.isotope_plot.object = hv.Text(0, 0, 'No isotope data')

    def view(self):
        """Return the full Panel layout."""
        return pn.Column(
            self.filter_bar,
            pn.Row(
                self.table,
                self.viz_panel,
                sizing_mode='stretch_both',
            ),
            sizing_mode='stretch_both',
        )


def create_app(store_path: Optional[str] = None,
               index_path: Optional[str] = None,
               raw_path: Optional[str] = None):
    """Create the Panel app."""

    # Load or create store
    if store_path and Path(store_path).exists():
        print(f"Loading store: {store_path}")
        store = load_store(store_path)
    elif index_path and raw_path:
        # Create store on-the-fly
        store_path = Path(raw_path).with_suffix('.zarr')
        print(f"Creating store: {store_path}")
        store = PrecursorStore.create_from_index_and_raw(
            index_path=index_path,
            raw_data_path=raw_path,
            output_path=str(store_path),
        )
    else:
        raise ValueError("Provide either --store or both --index and --raw")

    browser = PrecursorBrowser(store)
    return browser.view()


def main():
    parser = argparse.ArgumentParser(description="Precursor Browser GUI")
    parser.add_argument("--store", help="Path to .zarr store")
    parser.add_argument("--index", help="Path to precursor_index.parquet")
    parser.add_argument("--raw", help="Path to .d folder")
    parser.add_argument("--port", type=int, default=5006, help="Server port")
    parser.add_argument("--show", action="store_true", help="Open browser automatically")

    args = parser.parse_args()

    app = create_app(
        store_path=args.store,
        index_path=args.index,
        raw_path=args.raw,
    )

    # Serve the app
    pn.serve(app, port=args.port, show=args.show, title="Precursor Browser")


if __name__ == "__main__":
    main()


# For panel serve: create app from command line args
if __name__.startswith("bokeh"):
    import sys
    # Parse args passed via --args
    store_path = None
    index_path = None
    raw_path = None

    args = sys.argv[1:]
    for i, arg in enumerate(args):
        if arg == '--store' and i + 1 < len(args):
            store_path = args[i + 1]
        elif arg == '--index' and i + 1 < len(args):
            index_path = args[i + 1]
        elif arg == '--raw' and i + 1 < len(args):
            raw_path = args[i + 1]

    if store_path or (index_path and raw_path):
        app = create_app(store_path, index_path, raw_path)
        app.servable()
