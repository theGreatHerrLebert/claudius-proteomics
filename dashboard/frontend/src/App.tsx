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
      <div className="flex items-center gap-2 text-sm">
        <span className="text-white font-medium">Collection</span>
      </div>
    );
  }

  return (
    <div className="flex items-center gap-2 text-sm">
      <button
        onClick={onNavigateToCollection}
        className="text-blue-400 hover:text-blue-300"
      >
        Collection
      </button>
      <span className="text-gray-500">/</span>
      <span className="text-white font-medium">{activeDataset}</span>
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
  const [filters, setFilters] = useState({
    minEngines: 0,
    charge: undefined as number | undefined,
    rawFile: undefined as string | undefined,
    hasMs1: false,
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
        charge: filters.charge,
        raw_file: filters.rawFile,
        has_ms1: filters.hasMs1 || undefined,
        sort_by: filters.sortBy,
        sort_desc: filters.sortDesc,
      }),
  });

  // Fetch selected precursor detail
  const { data: selectedPrecursor, isLoading: isLoadingDetail } = useQuery({
    queryKey: ['precursor', selectedId],
    queryFn: () => (selectedId ? getPrecursor(selectedId) : null),
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
    <div className="h-screen flex flex-col bg-gray-900 text-gray-100">
      {/* Header */}
      <header className="flex-none h-14 bg-gray-800 border-b border-gray-700 flex items-center px-4 gap-6">
        <Breadcrumb
          mode="dataset"
          activeDataset={activeDataset}
          onNavigateToCollection={onNavigateToCollection}
        />

        {stats && (
          <div className="flex items-center gap-4 text-sm text-gray-400">
            <span>{stats.total_precursors.toLocaleString()} precursors</span>
            <span className="text-gray-600">|</span>
            <span>
              m/z: {stats.mz_range[0].toFixed(0)} - {stats.mz_range[1].toFixed(0)}
            </span>
          </div>
        )}

        {/* Filters */}
        <div className="flex items-center gap-3 ml-auto">
          <label className="flex items-center gap-2 text-sm">
            <span className="text-gray-400">Min engines:</span>
            <select
              value={filters.minEngines}
              onChange={(e) => {
                setFilters({ ...filters, minEngines: Number(e.target.value) });
                setPage(0);
              }}
              className="bg-gray-700 border border-gray-600 rounded px-2 py-1 text-sm"
            >
              <option value={0}>All</option>
              <option value={1}>≥1</option>
              <option value={2}>≥2</option>
              <option value={3}>3</option>
            </select>
          </label>

          <label className="flex items-center gap-2 text-sm">
            <span className="text-gray-400">Charge:</span>
            <select
              value={filters.charge ?? ''}
              onChange={(e) => {
                setFilters({
                  ...filters,
                  charge: e.target.value ? Number(e.target.value) : undefined,
                });
                setPage(0);
              }}
              className="bg-gray-700 border border-gray-600 rounded px-2 py-1 text-sm"
            >
              <option value="">All</option>
              <option value={2}>2+</option>
              <option value={3}>3+</option>
              <option value={4}>4+</option>
            </select>
          </label>

          {rawFiles && rawFiles.length > 1 && (
            <label className="flex items-center gap-2 text-sm">
              <span className="text-gray-400">File:</span>
              <select
                value={filters.rawFile ?? ''}
                onChange={(e) => {
                  setFilters({
                    ...filters,
                    rawFile: e.target.value || undefined,
                  });
                  setPage(0);
                }}
                className="bg-gray-700 border border-gray-600 rounded px-2 py-1 text-sm max-w-[200px]"
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
            <span className="text-gray-400">Sort:</span>
            <select
              value={filters.sortBy}
              onChange={(e) => setFilters({ ...filters, sortBy: e.target.value })}
              className="bg-gray-700 border border-gray-600 rounded px-2 py-1 text-sm"
            >
              <option value="n_engines">Engines</option>
              <option value="raw_intensity_meta">Intensity</option>
              <option value="mz">m/z</option>
              <option value="rt_seconds">RT</option>
              <option value="precursor_id">ID</option>
            </select>
          </label>

          <label className="flex items-center gap-2 text-sm cursor-pointer">
            <input
              type="checkbox"
              checked={filters.hasMs1}
              onChange={(e) => {
                setFilters({ ...filters, hasMs1: e.target.checked });
                setPage(0);
              }}
              className="rounded bg-gray-700 border-gray-600"
            />
            <span className="text-gray-400">Has MS1 data</span>
          </label>
        </div>
      </header>

      {/* Main content */}
      <div id="main-content" className={`flex-1 flex min-h-0 ${isResizing ? 'select-none' : ''}`}>
        {/* Table panel - resizable */}
        <div
          className="flex-none border-r border-gray-700 flex flex-col overflow-x-auto"
          style={{ width: `${tableWidth}%` }}
        >
          <div className="flex-1 min-h-0">
            <PrecursorTable
              data={precursors || []}
              selectedId={selectedId}
              onSelect={setSelectedId}
              isLoading={isLoadingList}
              sortBy={filters.sortBy}
              sortDesc={filters.sortDesc}
              onSort={handleSort}
            />
          </div>

          {/* Pagination */}
          <div className="flex-none h-10 bg-gray-800 border-t border-gray-700 flex items-center justify-between px-4">
            <button
              onClick={() => setPage(Math.max(0, page - 1))}
              disabled={page === 0}
              className="px-3 py-1 bg-gray-700 rounded text-sm disabled:opacity-50 disabled:cursor-not-allowed hover:bg-gray-600"
            >
              Previous
            </button>
            <span className="text-sm text-gray-400">
              Page {page + 1} ({page * pageSize + 1} - {page * pageSize + (precursors?.length || 0)})
            </span>
            <button
              onClick={() => setPage(page + 1)}
              disabled={(precursors?.length || 0) < pageSize}
              className="px-3 py-1 bg-gray-700 rounded text-sm disabled:opacity-50 disabled:cursor-not-allowed hover:bg-gray-600"
            >
              Next
            </button>
          </div>
        </div>

        {/* Resize handle */}
        <div
          className="w-2 bg-gray-600 hover:bg-blue-500 cursor-col-resize flex-none transition-colors active:bg-blue-600"
          onMouseDown={handleMouseDown}
          title="Drag to resize"
        >
          <div className="h-full w-full flex items-center justify-center">
            <div className="h-8 w-0.5 bg-gray-400 rounded" />
          </div>
        </div>

        {/* Visualization panel */}
        <div className="flex-1 min-h-0">
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
      <div className="h-screen flex items-center justify-center bg-gray-900 text-gray-100">
        <div className="text-gray-400">Loading...</div>
      </div>
    );
  }

  // Collection mode views
  if (isCollectionMode) {
    if (viewMode === 'collection') {
      return (
        <div className="h-screen flex flex-col bg-gray-900 text-gray-100">
          <CollectionBrowser onDatasetLoaded={handleDatasetLoaded} />
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
