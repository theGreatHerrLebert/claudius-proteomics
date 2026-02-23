import { useState, useEffect } from 'react';
import { QueryClient, QueryClientProvider, useQuery, useQueryClient } from '@tanstack/react-query';
import PrecursorTable from './components/PrecursorTable';
import PrecursorViz from './components/PrecursorViz';
import CollectionBrowser from './components/CollectionBrowser';
import {
  listPrecursors,
  getPrecursor,
  getStats,
  getRawFiles,
  getAppStatus,
} from './api';

const queryClient = new QueryClient();

type ViewMode = 'collection' | 'dataset';

interface BreadcrumbProps {
  mode: ViewMode;
  activeDataset: string | null;
  onNavigateToCollection: () => void;
}

function Breadcrumb({ mode, activeDataset, onNavigateToCollection }: BreadcrumbProps) {
  if (mode === 'collection') {
    return (
      <div className="flex items-center gap-2 text-sm reveal-up">
        <span className="metric-pill">
          <span className="subtle-label text-[0.62rem] text-[#9fb4d4]">View</span>
          <span className="metric-value">Collection</span>
        </span>
      </div>
    );
  }

  return (
    <div className="flex items-center gap-2 text-sm reveal-up">
      <button
        onClick={onNavigateToCollection}
        className="text-cyan-300 hover:text-cyan-100 transition-colors"
      >
        Collection
      </button>
      <span className="text-slate-500">/</span>
      <span className="metric-pill">
        <span className="subtle-label text-[0.62rem] text-[#9fb4d4]">Dataset</span>
        <span className="mono text-[0.77rem] tracking-wide">{activeDataset}</span>
      </span>
    </div>
  );
}

function DatasetView({
  activeDataset,
  onNavigateToCollection,
}: {
  activeDataset: string | null;
  onNavigateToCollection: () => void;
}) {
  const [selectedId, setSelectedId] = useState<number | null>(null);
  const [selectedRawFile, setSelectedRawFile] = useState<string | undefined>(undefined);
  const [filters, setFilters] = useState({
    minEngines: 0,
    maxEngines: undefined as number | undefined,
    charge: undefined as number | undefined,
    rawFile: undefined as string | undefined,
    hasMs1: false,
    hasRawData: false,
    sortBy: 'n_engines',
    sortDesc: true,
  });
  const [page, setPage] = useState(0);
  const pageSize = 100;

  // Resizable panel state
  const [tableWidth, setTableWidth] = useState(60); // percentage
  const [isResizing, setIsResizing] = useState(false);

  // Handle resize drag
  const handleMouseDown = (e: React.MouseEvent) => {
    e.preventDefault();
    setIsResizing(true);
  };

  useEffect(() => {
    const handleMouseMove = (e: MouseEvent) => {
      if (!isResizing) return;
      const container = document.getElementById('main-content');
      if (!container) return;
      const rect = container.getBoundingClientRect();
      const newWidth = ((e.clientX - rect.left) / rect.width) * 100;
      setTableWidth(Math.min(Math.max(newWidth, 20), 80)); // Clamp between 20-80%
    };

    const handleMouseUp = () => {
      setIsResizing(false);
    };

    if (isResizing) {
      document.addEventListener('mousemove', handleMouseMove);
      document.addEventListener('mouseup', handleMouseUp);
    }

    return () => {
      document.removeEventListener('mousemove', handleMouseMove);
      document.removeEventListener('mouseup', handleMouseUp);
    };
  }, [isResizing]);

  // Handle sorting - clicking same column toggles direction, new column sorts desc
  const handleSort = (column: string) => {
    setFilters((prev) => ({
      ...prev,
      sortBy: column,
      sortDesc: prev.sortBy === column ? !prev.sortDesc : true,
    }));
    setPage(0); // Reset to first page on sort change
  };

  // Fetch precursor list
  const {
    data: precursors,
    isLoading: isLoadingList,
  } = useQuery({
    queryKey: ['precursors', filters, page, activeDataset],
    queryFn: () =>
      listPrecursors({
        offset: page * pageSize,
        limit: pageSize,
        min_engines: filters.minEngines,
        max_engines: filters.maxEngines,
        charge: filters.charge,
        raw_file: filters.rawFile,
        has_ms1: filters.hasMs1 || undefined,
        has_raw_data: filters.hasRawData || undefined,
        sort_by: filters.sortBy,
        sort_desc: filters.sortDesc,
      }),
  });

  // Fetch selected precursor detail
  const { data: selectedPrecursor, isLoading: isLoadingDetail } = useQuery({
    queryKey: ['precursor', selectedId, selectedRawFile],
    queryFn: () => (selectedId ? getPrecursor(selectedId, selectedRawFile) : null),
    enabled: selectedId !== null,
  });

  // Fetch stats
  const { data: stats } = useQuery({
    queryKey: ['stats', activeDataset],
    queryFn: getStats,
  });

  // Fetch available raw files
  const { data: rawFiles } = useQuery({
    queryKey: ['raw_files', activeDataset],
    queryFn: getRawFiles,
  });

  return (
    <div className="app-shell h-screen flex flex-col">
      {/* Header */}
      <header className="chrome-header flex-none px-3 py-2 md:px-4 md:py-3 flex flex-wrap items-center gap-2 md:gap-4">
        <Breadcrumb
          mode="dataset"
          activeDataset={activeDataset}
          onNavigateToCollection={onNavigateToCollection}
        />

        {stats && (
          <div className="flex items-center gap-2 text-sm reveal-up">
            <span className="metric-pill">
              Precursors <span className="metric-value">{stats.total_precursors.toLocaleString()}</span>
            </span>
            <span className="metric-pill">
              m/z <span className="metric-value">{stats.mz_range[0].toFixed(0)}-{stats.mz_range[1].toFixed(0)}</span>
            </span>
          </div>
        )}

        {/* Filters */}
        <div className="flex flex-wrap items-center gap-2 ml-auto">
          <label className="flex items-center gap-2 text-sm">
            <span className="control-label">Engines</span>
            <select
              value={filters.maxEngines === 0 ? 'unid' : (filters.minEngines > 0 ? `min${filters.minEngines}` : 'all')}
              onChange={(e) => {
                const val = e.target.value;
                if (val === 'unid') {
                  setFilters({ ...filters, minEngines: 0, maxEngines: 0 });
                } else if (val === 'all') {
                  setFilters({ ...filters, minEngines: 0, maxEngines: undefined });
                } else {
                  setFilters({ ...filters, minEngines: Number(val.slice(3)), maxEngines: undefined });
                }
                setPage(0);
              }}
              className="control-select"
            >
              <option value="all">All</option>
              <option value="unid">Unidentified (0)</option>
              <option value="min1">≥1</option>
              <option value="min2">≥2</option>
              <option value="min3">3 only</option>
            </select>
          </label>

          <label className="flex items-center gap-2 text-sm">
            <span className="control-label">Charge</span>
            <select
              value={filters.charge ?? ''}
              onChange={(e) => {
                setFilters({
                  ...filters,
                  charge: e.target.value ? Number(e.target.value) : undefined,
                });
                setPage(0);
              }}
              className="control-select"
            >
              <option value="">All</option>
              <option value={1}>1+</option>
              <option value={2}>2+</option>
              <option value={3}>3+</option>
              <option value={4}>4+</option>
              <option value={5}>5+</option>
            </select>
          </label>

          {rawFiles && rawFiles.length > 1 && (
            <label className="flex items-center gap-2 text-sm">
              <span className="control-label">File</span>
              <select
                value={filters.rawFile ?? ''}
                onChange={(e) => {
                  setFilters({
                    ...filters,
                    rawFile: e.target.value || undefined,
                  });
                  setPage(0);
                }}
                className="control-select max-w-[220px]"
              >
                <option value="">All files</option>
                {rawFiles.map((f) => (
                  <option key={f.name} value={f.name}>
                    {f.name.slice(-30)} ({f.count.toLocaleString()})
                  </option>
                ))}
              </select>
            </label>
          )}

          <label className="flex items-center gap-2 text-sm">
            <span className="control-label">Sort</span>
            <select
              value={filters.sortBy}
              onChange={(e) => setFilters({ ...filters, sortBy: e.target.value })}
              className="control-select"
            >
              <option value="n_engines">Engines</option>
              <option value="raw_intensity_meta">Intensity</option>
              <option value="mz">m/z</option>
              <option value="rt_seconds">RT</option>
              <option value="precursor_id">ID</option>
            </select>
          </label>

          <label className="control-check flex items-center gap-2 text-sm cursor-pointer rounded px-2 py-1">
            <input
              type="checkbox"
              checked={filters.hasMs1}
              onChange={(e) => {
                setFilters({ ...filters, hasMs1: e.target.checked });
                setPage(0);
              }}
              className="rounded border-slate-500/70 bg-slate-900/80"
            />
            <span>MS1 data</span>
          </label>

          <label className="control-check flex items-center gap-2 text-sm cursor-pointer rounded px-2 py-1">
            <input
              type="checkbox"
              checked={filters.hasRawData}
              onChange={(e) => {
                setFilters({ ...filters, hasRawData: e.target.checked });
                setPage(0);
              }}
              className="rounded border-slate-500/70 bg-slate-900/80"
            />
            <span>Raw data</span>
          </label>
        </div>
      </header>

      {/* Main content */}
      <div id="main-content" className={`flex-1 flex min-h-0 p-2 gap-2 ${isResizing ? 'select-none' : ''}`}>
        {/* Table panel - resizable */}
        <div
          className="chrome-panel flex-none flex flex-col overflow-x-auto min-w-[280px]"
          style={{ width: `${tableWidth}%` }}
        >
          <div className="flex-1 min-h-0">
            <PrecursorTable
              data={precursors || []}
              selectedId={selectedId}
              onSelect={(id, rawFile) => { setSelectedId(id); setSelectedRawFile(rawFile); }}
              isLoading={isLoadingList}
              sortBy={filters.sortBy}
              sortDesc={filters.sortDesc}
              onSort={handleSort}
            />
          </div>

          {/* Pagination */}
          <div className="flex-none h-11 bg-slate-950/25 border-t border-[#25415f] flex items-center justify-between px-3">
            <button
              onClick={() => setPage(Math.max(0, page - 1))}
              disabled={page === 0}
              className="btn-secondary"
            >
              Previous
            </button>
            <span className="text-xs md:text-sm text-slate-300">
              Page {page + 1} ({page * pageSize + 1} - {page * pageSize + (precursors?.length || 0)})
            </span>
            <button
              onClick={() => setPage(page + 1)}
              disabled={(precursors?.length || 0) < pageSize}
              className="btn-secondary"
            >
              Next
            </button>
          </div>
        </div>

        {/* Resize handle */}
        <div
          className="split-handle w-2 rounded cursor-col-resize flex-none transition-colors active:brightness-110"
          onMouseDown={handleMouseDown}
          title="Drag to resize"
        >
          <div className="h-full w-full flex items-center justify-center">
            <div className="h-9 w-0.5 bg-cyan-200/50 rounded" />
          </div>
        </div>

        {/* Visualization panel */}
        <div className="chrome-panel flex-1 min-h-0 overflow-hidden">
          <PrecursorViz precursor={selectedPrecursor || null} isLoading={isLoadingDetail} />
        </div>
      </div>
    </div>
  );
}

function Dashboard() {
  const [viewMode, setViewMode] = useState<ViewMode>('dataset');
  const [activeDataset, setActiveDataset] = useState<string | null>(null);
  const [isCollectionMode, setIsCollectionMode] = useState(false);
  const queryClientHook = useQueryClient();

  // Check app status on mount
  const { data: appStatus, isLoading: isLoadingStatus } = useQuery({
    queryKey: ['appStatus'],
    queryFn: getAppStatus,
    retry: false,
  });

  // Update state when app status loads
  useEffect(() => {
    if (appStatus) {
      setIsCollectionMode(appStatus.mode === 'collection');
      setActiveDataset(appStatus.active_dataset);
      // If collection mode but no dataset loaded, show collection browser
      if (appStatus.mode === 'collection' && !appStatus.store_loaded) {
        setViewMode('collection');
      } else {
        setViewMode('dataset');
      }
    }
  }, [appStatus]);

  const handleDatasetLoaded = (accession: string) => {
    setActiveDataset(accession);
    setViewMode('dataset');
    // Invalidate queries to refresh with new data
    queryClientHook.invalidateQueries({ queryKey: ['precursors'] });
    queryClientHook.invalidateQueries({ queryKey: ['stats'] });
    queryClientHook.invalidateQueries({ queryKey: ['raw_files'] });
  };

  const handleNavigateToCollection = () => {
    setViewMode('collection');
  };

  if (isLoadingStatus) {
    return (
      <div className="app-shell h-screen flex items-center justify-center">
        <div className="metric-pill">Loading dashboard status...</div>
      </div>
    );
  }

  // Collection mode views
  if (isCollectionMode) {
    if (viewMode === 'collection') {
      return (
        <div className="app-shell h-screen flex flex-col surface-grid">
          <CollectionBrowser
            onDatasetLoaded={handleDatasetLoaded}
            activeDataset={activeDataset}
            onReturnToDataset={activeDataset ? () => setViewMode('dataset') : undefined}
          />
        </div>
      );
    }
    return (
      <DatasetView
        activeDataset={activeDataset}
        onNavigateToCollection={handleNavigateToCollection}
      />
    );
  }

  // Single dataset mode (legacy)
  return (
    <DatasetView
      activeDataset={null}
      onNavigateToCollection={() => {}}
    />
  );
}

function App() {
  return (
    <QueryClientProvider client={queryClient}>
      <Dashboard />
    </QueryClientProvider>
  );
}

export default App;
