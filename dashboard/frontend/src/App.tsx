import { useState, useEffect, useMemo } from 'react';
import { QueryClient, QueryClientProvider, useQuery, useQueryClient } from '@tanstack/react-query';
import PrecursorTable from './components/PrecursorTable';
import PrecursorViz from './components/PrecursorViz';
import CollectionBrowser from './components/CollectionBrowser';
import DatasetSummaryPage from './components/DatasetSummary';
import LandingPage from './components/LandingPage';
import VisitSummaryPage from './components/VisitSummaryPage';
import BlueprintPage from './components/BlueprintPage';
import {
  listPrecursors,
  getPrecursor,
  getStats,
  getRawFiles,
  getAppStatus,
  getStudies,
  getDatasetInfo,
  getDatasetOverlap,
} from './api';
import type { StudySummary } from './api';

const queryClient = new QueryClient();

type ViewMode = 'landing' | 'visit' | 'blueprint' | 'collection' | 'dataset';

// External URL for the data explorer (hosted on a separate machine).
// Set to '' to fall back to in-app collection browser navigation.
const EXPLORER_URL = '';
type SubView = 'summary' | 'browse';
type LinkCopyState = 'idle' | 'copied' | 'failed';
type OverlapFocusMode = 'all' | 'shared' | 'unique_a';

interface DatasetFilters {
  minEngines: number;
  maxEngines: number | undefined;
  charge: number | undefined;
  rawFile: string | undefined;
  hasMs1: boolean;
  hasRawData: boolean;
  sortBy: string;
  sortDesc: boolean;
}

interface SavedView {
  id: string;
  name: string;
  filters: DatasetFilters;
  page: number;
  tableWidth: number;
  compareEnabled: boolean;
  compareTarget: CompareTarget | null;
  datasetCompareEnabled: boolean;
  datasetCompareTarget: string | null;
  overlapFocusMode: OverlapFocusMode;
  createdAt: string;
}

interface DatasetUrlState {
  subView: SubView;
  filters: DatasetFilters;
  page: number;
  tableWidth: number;
  compareEnabled: boolean;
  compareTarget: CompareTarget | null;
  datasetCompareEnabled: boolean;
  datasetCompareTarget: string | null;
  overlapFocusMode: OverlapFocusMode;
}

interface CompareTarget {
  id: number;
  rawFile: string | undefined;
}

interface OverlapVisualProps {
  title: string;
  totalA: number;
  totalB: number;
  shared: number;
  uniqueA: number;
  uniqueB: number;
  jaccard: number;
  labelA: string;
  labelB: string;
  activeFocus?: OverlapFocusMode;
  onFocusShared?: () => void;
  onFocusUniqueA?: () => void;
  onClearFocus?: () => void;
}

const SAVED_VIEWS_STORAGE_KEY = 'sanjose.dashboard.savedViews';
const DEFAULT_TABLE_WIDTH = 60;
const DEFAULT_FILTERS: DatasetFilters = {
  minEngines: 0,
  maxEngines: undefined,
  charge: undefined,
  rawFile: undefined,
  hasMs1: false,
  hasRawData: false,
  sortBy: 'n_engines',
  sortDesc: true,
};

function parseNumberParam(value: string | null): number | undefined {
  if (!value) return undefined;
  const parsed = Number(value);
  if (Number.isNaN(parsed)) return undefined;
  return parsed;
}

function clampTableWidth(width: number): number {
  return Math.min(Math.max(width, 20), 80);
}

function buildCompareKey(target: CompareTarget | null): string {
  if (!target) return '';
  return `${target.id}::${target.rawFile ?? ''}`;
}

function parseCompareKey(value: string): CompareTarget | null {
  if (!value) return null;
  const [idPart, rawFilePart] = value.split('::');
  if (!idPart) return null;
  const parsedId = Number(idPart);
  if (!Number.isFinite(parsedId) || parsedId <= 0) return null;
  return {
    id: parsedId,
    rawFile: rawFilePart || undefined,
  };
}

function formatSigned(value: number, digits: number): string {
  const abs = Math.abs(value).toFixed(digits);
  return `${value >= 0 ? '+' : '-'}${abs}`;
}

function formatNullableNumber(value: number | null | undefined, digits = 2): string {
  if (value === null || value === undefined) return '-';
  return value.toFixed(digits);
}

function formatNullableInt(value: number | null | undefined): string {
  if (value === null || value === undefined) return '-';
  return value.toLocaleString();
}

function formatPercent(value: number): string {
  return `${(value * 100).toFixed(1)}%`;
}

function OverlapVisualCard({
  title,
  totalA,
  totalB,
  shared,
  uniqueA,
  uniqueB,
  jaccard,
  labelA,
  labelB,
  activeFocus,
  onFocusShared,
  onFocusUniqueA,
  onClearFocus,
}: OverlapVisualProps) {
  const union = Math.max(1, uniqueA + shared + uniqueB);
  const aCoverage = totalA > 0 ? shared / totalA : 0;
  const bCoverage = totalB > 0 ? shared / totalB : 0;

  const sharedWidth = `${Math.max((shared / union) * 100, shared > 0 ? 2 : 0)}%`;
  const uniqueAWidth = `${Math.max((uniqueA / union) * 100, uniqueA > 0 ? 2 : 0)}%`;
  const uniqueBWidth = `${Math.max((uniqueB / union) * 100, uniqueB > 0 ? 2 : 0)}%`;

  return (
    <div className="panel-inset p-3">
      <div className="flex items-center justify-between mb-2">
        <div className="control-label">{title}</div>
        <div className="metric-pill">Jaccard <span className="metric-value mono">{formatPercent(jaccard)}</span></div>
      </div>

      {(onFocusShared || onFocusUniqueA || onClearFocus) && (
        <div className="flex flex-wrap items-center gap-2 mb-2">
          {onFocusShared && (
            <button
              type="button"
              onClick={onFocusShared}
              className={`btn-secondary ${activeFocus === 'shared' ? 'border-cyan-300 text-cyan-100' : ''}`}
            >
              Focus Shared
            </button>
          )}
          {onFocusUniqueA && (
            <button
              type="button"
              onClick={onFocusUniqueA}
              className={`btn-secondary ${activeFocus === 'unique_a' ? 'border-cyan-300 text-cyan-100' : ''}`}
            >
              Focus Unique {labelA}
            </button>
          )}
          {onClearFocus && (
            <button
              type="button"
              onClick={onClearFocus}
              className={`btn-secondary ${activeFocus === 'all' ? 'border-cyan-300 text-cyan-100' : ''}`}
            >
              Clear Focus
            </button>
          )}
        </div>
      )}

      <div className="h-5 rounded-md overflow-hidden border border-[#2a4368] flex bg-slate-900/70">
        <div className="bg-cyan-700/80 h-full" style={{ width: uniqueAWidth }} title={`Unique ${labelA}: ${uniqueA.toLocaleString()}`} />
        <div className="bg-emerald-500/85 h-full" style={{ width: sharedWidth }} title={`Shared: ${shared.toLocaleString()}`} />
        <div className="bg-orange-500/80 h-full" style={{ width: uniqueBWidth }} title={`Unique ${labelB}: ${uniqueB.toLocaleString()}`} />
      </div>

      <div className="grid grid-cols-3 gap-2 text-xs mt-2">
        <div className="metric-pill justify-between">Unique {labelA} <span className="metric-value mono">{uniqueA.toLocaleString()}</span></div>
        <div className="metric-pill justify-between">Shared <span className="metric-value mono">{shared.toLocaleString()}</span></div>
        <div className="metric-pill justify-between">Unique {labelB} <span className="metric-value mono">{uniqueB.toLocaleString()}</span></div>
      </div>

      <div className="mt-3 grid grid-cols-1 gap-2">
        <div>
          <div className="flex items-center justify-between text-[11px] text-slate-300 mb-1">
            <span>{labelA} coverage by shared</span>
            <span className="mono">{formatPercent(aCoverage)}</span>
          </div>
          <div className="h-2 rounded bg-slate-800 overflow-hidden">
            <div className="h-2 bg-cyan-400/85" style={{ width: `${aCoverage * 100}%` }} />
          </div>
          <div className="text-[10px] text-slate-400 mt-0.5">{shared.toLocaleString()} / {totalA.toLocaleString()}</div>
        </div>
        <div>
          <div className="flex items-center justify-between text-[11px] text-slate-300 mb-1">
            <span>{labelB} coverage by shared</span>
            <span className="mono">{formatPercent(bCoverage)}</span>
          </div>
          <div className="h-2 rounded bg-slate-800 overflow-hidden">
            <div className="h-2 bg-orange-400/85" style={{ width: `${bCoverage * 100}%` }} />
          </div>
          <div className="text-[10px] text-slate-400 mt-0.5">{shared.toLocaleString()} / {totalB.toLocaleString()}</div>
        </div>
      </div>
    </div>
  );
}

function areCompareTargetsEqual(a: CompareTarget | null, b: CompareTarget | null): boolean {
  if (a === null && b === null) return true;
  if (a === null || b === null) return false;
  return a.id === b.id && a.rawFile === b.rawFile;
}

function areDatasetTargetsEqual(a: string | null, b: string | null): boolean {
  return a === b;
}

function areFiltersEqual(a: DatasetFilters, b: DatasetFilters): boolean {
  return (
    a.minEngines === b.minEngines &&
    a.maxEngines === b.maxEngines &&
    a.charge === b.charge &&
    a.rawFile === b.rawFile &&
    a.hasMs1 === b.hasMs1 &&
    a.hasRawData === b.hasRawData &&
    a.sortBy === b.sortBy &&
    a.sortDesc === b.sortDesc
  );
}

function readDatasetStateFromUrl(): DatasetUrlState {
  const defaultState: DatasetUrlState = {
    subView: 'summary',
    filters: { ...DEFAULT_FILTERS },
    page: 0,
    tableWidth: DEFAULT_TABLE_WIDTH,
    compareEnabled: false,
    compareTarget: null,
    datasetCompareEnabled: false,
    datasetCompareTarget: null,
    overlapFocusMode: 'all' as OverlapFocusMode,
  };
  if (typeof window === 'undefined') return defaultState;

  const params = new URLSearchParams(window.location.search);
  const viewParam = params.get('view');
  const subView: SubView = viewParam === 'browse' ? 'browse' : 'summary';
  const filters: DatasetFilters = { ...DEFAULT_FILTERS };

  const minEngines = parseNumberParam(params.get('minEng'));
  if (minEngines !== undefined && minEngines >= 0) {
    filters.minEngines = minEngines;
  }

  const maxEngines = parseNumberParam(params.get('maxEng'));
  if (maxEngines !== undefined && maxEngines >= 0) {
    filters.maxEngines = maxEngines;
  }

  const charge = parseNumberParam(params.get('charge'));
  if (charge !== undefined && charge > 0) {
    filters.charge = charge;
  }

  const rawFile = params.get('rawFile');
  if (rawFile) {
    filters.rawFile = rawFile;
  }

  filters.hasMs1 = params.get('ms1') === '1';
  filters.hasRawData = params.get('rawData') === '1';

  const sortBy = params.get('sort');
  if (sortBy) {
    filters.sortBy = sortBy;
  }

  const sortDesc = params.get('desc');
  if (sortDesc === '0') {
    filters.sortDesc = false;
  }

  const pageParam = parseNumberParam(params.get('page'));
  const page = pageParam !== undefined && pageParam >= 0 ? Math.floor(pageParam) : 0;

  const widthParam = parseNumberParam(params.get('table'));
  const tableWidth = widthParam !== undefined ? clampTableWidth(widthParam) : DEFAULT_TABLE_WIDTH;

  const compareEnabled = params.get('cmp') === '1';
  const compareId = parseNumberParam(params.get('cmpId'));
  const compareRaw = params.get('cmpRaw');
  const compareTarget = compareId !== undefined
    ? {
      id: Math.floor(compareId),
      rawFile: compareRaw || undefined,
    }
    : null;

  const datasetCompareEnabled = params.get('dcmp') === '1';
  const datasetCompareTarget = params.get('dcmpDs');
  const overlapFocusRaw = params.get('ofm');
  const overlapFocusMode: OverlapFocusMode =
    overlapFocusRaw === 'shared' || overlapFocusRaw === 'unique_a' ? overlapFocusRaw : 'all';

  return {
    subView,
    filters,
    page,
    tableWidth,
    compareEnabled,
    compareTarget,
    datasetCompareEnabled,
    datasetCompareTarget: datasetCompareTarget || null,
    overlapFocusMode,
  };
}

function writeDatasetStateToUrl({
  subView,
  filters,
  page,
  tableWidth,
  activeDataset,
  compareEnabled,
  compareTarget,
  datasetCompareEnabled,
  datasetCompareTarget,
  overlapFocusMode,
}: {
  subView: SubView;
  filters: DatasetFilters;
  page: number;
  tableWidth: number;
  activeDataset: string | null;
  compareEnabled: boolean;
  compareTarget: CompareTarget | null;
  datasetCompareEnabled: boolean;
  datasetCompareTarget: string | null;
  overlapFocusMode: OverlapFocusMode;
}) {
  if (typeof window === 'undefined') return;

  const params = new URLSearchParams(window.location.search);
  const setParam = (key: string, value: string | undefined) => {
    if (value === undefined || value === '') {
      params.delete(key);
      return;
    }
    params.set(key, value);
  };

  setParam('view', subView !== 'summary' ? subView : undefined);
  setParam('ds', activeDataset ?? undefined);
  setParam('minEng', filters.minEngines > 0 ? String(filters.minEngines) : undefined);
  setParam('maxEng', filters.maxEngines !== undefined ? String(filters.maxEngines) : undefined);
  setParam('charge', filters.charge !== undefined ? String(filters.charge) : undefined);
  setParam('rawFile', filters.rawFile);
  setParam('ms1', filters.hasMs1 ? '1' : undefined);
  setParam('rawData', filters.hasRawData ? '1' : undefined);
  setParam('sort', filters.sortBy !== DEFAULT_FILTERS.sortBy ? filters.sortBy : undefined);
  setParam('desc', filters.sortDesc ? undefined : '0');
  setParam('page', page > 0 ? String(page) : undefined);
  setParam('table', tableWidth !== DEFAULT_TABLE_WIDTH ? String(Math.round(tableWidth)) : undefined);
  setParam('cmp', compareEnabled ? '1' : undefined);
  setParam('cmpId', compareEnabled && compareTarget ? String(compareTarget.id) : undefined);
  setParam('cmpRaw', compareEnabled && compareTarget ? compareTarget.rawFile : undefined);
  setParam('dcmp', datasetCompareEnabled ? '1' : undefined);
  setParam('dcmpDs', datasetCompareEnabled ? datasetCompareTarget ?? undefined : undefined);
  setParam('ofm', overlapFocusMode !== 'all' ? overlapFocusMode : undefined);

  const nextQuery = params.toString();
  const currentQuery = window.location.search.startsWith('?')
    ? window.location.search.slice(1)
    : window.location.search;
  if (nextQuery === currentQuery) return;

  const nextUrl = `${window.location.pathname}${nextQuery ? `?${nextQuery}` : ''}${window.location.hash}`;
  window.history.replaceState(null, '', nextUrl);
}

function loadSavedViewsFromStorage(): SavedView[] {
  if (typeof window === 'undefined') return [];
  try {
    const raw = window.localStorage.getItem(SAVED_VIEWS_STORAGE_KEY);
    if (!raw) return [];
    const parsed: unknown = JSON.parse(raw);
    if (!Array.isArray(parsed)) return [];
    return parsed
      .filter((item): item is SavedView =>
        typeof item === 'object' &&
        item !== null &&
        typeof (item as SavedView).id === 'string' &&
        typeof (item as SavedView).name === 'string' &&
        typeof (item as SavedView).createdAt === 'string' &&
        typeof (item as SavedView).page === 'number' &&
        typeof (item as SavedView).tableWidth === 'number' &&
        typeof (item as SavedView).filters === 'object' &&
        (item as SavedView).filters !== null
      )
      .map((item) => ({
        ...item,
        compareEnabled: item.compareEnabled ?? false,
        compareTarget: item.compareTarget ?? null,
        datasetCompareEnabled: item.datasetCompareEnabled ?? false,
        datasetCompareTarget: item.datasetCompareTarget ?? null,
        overlapFocusMode: item.overlapFocusMode ?? 'all',
      }))
      .slice(0, 25);
  } catch {
    return [];
  }
}

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
  const [initialUrlState] = useState<DatasetUrlState>(readDatasetStateFromUrl);
  const [subView, setSubView] = useState<SubView>(initialUrlState.subView);
  const [selectedId, setSelectedId] = useState<number | null>(null);
  const [selectedRawFile, setSelectedRawFile] = useState<string | undefined>(undefined);
  const [filters, setFilters] = useState<DatasetFilters>(initialUrlState.filters);
  const [page, setPage] = useState(initialUrlState.page);
  const pageSize = 100;

  // Resizable panel state
  const [tableWidth, setTableWidth] = useState(initialUrlState.tableWidth); // percentage
  const [isResizing, setIsResizing] = useState(false);
  const [savedViews, setSavedViews] = useState<SavedView[]>(loadSavedViewsFromStorage);
  const [selectedViewId, setSelectedViewId] = useState<string>('');
  const [linkCopyState, setLinkCopyState] = useState<LinkCopyState>('idle');
  const [compareEnabled, setCompareEnabled] = useState(initialUrlState.compareEnabled);
  const [compareTarget, setCompareTarget] = useState<CompareTarget | null>(initialUrlState.compareTarget);
  const [datasetCompareEnabled, setDatasetCompareEnabled] = useState(initialUrlState.datasetCompareEnabled);
  const [datasetCompareTarget, setDatasetCompareTarget] = useState<string | null>(initialUrlState.datasetCompareTarget);
  const [overlapFocusMode, setOverlapFocusMode] = useState<OverlapFocusMode>(initialUrlState.overlapFocusMode);
  const [showAdvancedControls, setShowAdvancedControls] = useState(false);

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

  useEffect(() => {
    if (isResizing) return;
    writeDatasetStateToUrl({
      subView,
      filters,
      page,
      tableWidth,
      activeDataset,
      compareEnabled,
      compareTarget,
      datasetCompareEnabled,
      datasetCompareTarget,
      overlapFocusMode,
    });
  }, [subView, filters, page, tableWidth, activeDataset, compareEnabled, compareTarget, datasetCompareEnabled, datasetCompareTarget, overlapFocusMode, isResizing]);

  const persistSavedViews = (nextViews: SavedView[]) => {
    setSavedViews(nextViews);
    if (typeof window !== 'undefined') {
      window.localStorage.setItem(SAVED_VIEWS_STORAGE_KEY, JSON.stringify(nextViews));
    }
  };

  const handleResetToDefaultView = () => {
    setFilters({ ...DEFAULT_FILTERS });
    setPage(0);
    setTableWidth(DEFAULT_TABLE_WIDTH);
    setCompareEnabled(false);
    setCompareTarget(null);
    setDatasetCompareEnabled(false);
    setDatasetCompareTarget(null);
    setOverlapFocusMode('all');
    setSelectedViewId('');
  };

  const handleApplySavedView = (id: string) => {
    const target = savedViews.find((view) => view.id === id);
    if (!target) return;
    setFilters({ ...target.filters });
    setPage(Math.max(0, Math.floor(target.page)));
    setTableWidth(clampTableWidth(target.tableWidth));
    setCompareEnabled(target.compareEnabled);
    setCompareTarget(target.compareTarget);
    setDatasetCompareEnabled(target.datasetCompareEnabled);
    setDatasetCompareTarget(target.datasetCompareTarget);
    setOverlapFocusMode(target.overlapFocusMode);
    setSelectedViewId(id);
  };

  const handleSaveView = () => {
    if (typeof window === 'undefined') return;
    const suggestedName = activeDataset ? `${activeDataset} view` : `Saved view ${savedViews.length + 1}`;
    const name = window.prompt('Name this saved view:', suggestedName);
    if (!name || !name.trim()) return;

    const entry: SavedView = {
      id: `view-${Date.now()}-${Math.floor(Math.random() * 100000)}`,
      name: name.trim(),
      filters: { ...filters },
      page,
      tableWidth,
      compareEnabled,
      compareTarget,
      datasetCompareEnabled,
      datasetCompareTarget,
      overlapFocusMode,
      createdAt: new Date().toISOString(),
    };
    const nextViews = [entry, ...savedViews].slice(0, 25);
    persistSavedViews(nextViews);
    setSelectedViewId(entry.id);
  };

  const handleDeleteView = () => {
    if (typeof window === 'undefined' || !selectedViewId) return;
    const target = savedViews.find((view) => view.id === selectedViewId);
    if (!target) return;
    const confirmed = window.confirm(`Delete saved view "${target.name}"?`);
    if (!confirmed) return;

    const nextViews = savedViews.filter((view) => view.id !== selectedViewId);
    persistSavedViews(nextViews);
    setSelectedViewId('');
  };

  const handleCopyLink = async () => {
    if (typeof window === 'undefined') return;
    try {
      if (navigator.clipboard?.writeText) {
        await navigator.clipboard.writeText(window.location.href);
      } else {
        const helper = document.createElement('textarea');
        helper.value = window.location.href;
        helper.style.position = 'fixed';
        helper.style.opacity = '0';
        document.body.appendChild(helper);
        helper.select();
        document.execCommand('copy');
        document.body.removeChild(helper);
      }
      setLinkCopyState('copied');
    } catch {
      setLinkCopyState('failed');
    } finally {
      window.setTimeout(() => {
        setLinkCopyState('idle');
      }, 1800);
    }
  };

  const handleSetCurrentAsCompare = () => {
    if (selectedId === null) return;
    setCompareEnabled(true);
    setCompareTarget({
      id: selectedId,
      rawFile: selectedRawFile,
    });
  };

  const handleOverlapFocusChange = (mode: OverlapFocusMode) => {
    setOverlapFocusMode(mode);
    setPage(0);
  };

  useEffect(() => {
    if (!activeDataset) {
      if (datasetCompareEnabled) setDatasetCompareEnabled(false);
      if (datasetCompareTarget !== null) setDatasetCompareTarget(null);
      if (overlapFocusMode !== 'all') setOverlapFocusMode('all');
      return;
    }
    if (datasetCompareTarget === activeDataset) {
      setDatasetCompareTarget(null);
    }
  }, [activeDataset, datasetCompareEnabled, datasetCompareTarget, overlapFocusMode]);

  useEffect(() => {
    if (!datasetCompareEnabled || !datasetCompareTarget) {
      if (overlapFocusMode !== 'all') setOverlapFocusMode('all');
    }
  }, [datasetCompareEnabled, datasetCompareTarget, overlapFocusMode]);

  useEffect(() => {
    if (subView !== 'browse' && showAdvancedControls) {
      setShowAdvancedControls(false);
    }
  }, [subView, showAdvancedControls]);

  useEffect(() => {
    const matching = savedViews.find(
      (view) =>
        areFiltersEqual(view.filters, filters) &&
        view.page === page &&
        clampTableWidth(view.tableWidth) === clampTableWidth(tableWidth) &&
        view.compareEnabled === compareEnabled &&
        areCompareTargetsEqual(view.compareTarget, compareTarget) &&
        view.datasetCompareEnabled === datasetCompareEnabled &&
        areDatasetTargetsEqual(view.datasetCompareTarget, datasetCompareTarget) &&
        view.overlapFocusMode === overlapFocusMode
    );
    if (!matching) {
      if (selectedViewId) setSelectedViewId('');
      return;
    }
    if (matching.id !== selectedViewId) {
      setSelectedViewId(matching.id);
    }
  }, [filters, page, tableWidth, compareEnabled, compareTarget, datasetCompareEnabled, datasetCompareTarget, overlapFocusMode, savedViews, selectedViewId]);

  // Reset to summary view when a new dataset is loaded
  const [prevDataset, setPrevDataset] = useState(activeDataset);
  if (activeDataset !== prevDataset) {
    setPrevDataset(activeDataset);
    setSubView('summary');
  }

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
    queryKey: ['precursors', filters, page, activeDataset, overlapFocusMode, datasetCompareEnabled, datasetCompareTarget],
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
        overlap_mode: datasetCompareEnabled && datasetCompareTarget && overlapFocusMode !== 'all'
          ? overlapFocusMode
          : undefined,
        overlap_dataset: datasetCompareEnabled && datasetCompareTarget && overlapFocusMode !== 'all'
          ? datasetCompareTarget
          : undefined,
      }),
  });

  // Fetch selected precursor detail
  const { data: selectedPrecursor, isLoading: isLoadingDetail } = useQuery({
    queryKey: ['precursor', selectedId, selectedRawFile],
    queryFn: () => (selectedId ? getPrecursor(selectedId, selectedRawFile) : null),
    enabled: selectedId !== null,
  });

  const {
    data: comparisonPrecursor,
    isLoading: isLoadingComparison,
  } = useQuery({
    queryKey: ['precursor-compare', compareTarget?.id, compareTarget?.rawFile],
    queryFn: () => (compareTarget ? getPrecursor(compareTarget.id, compareTarget.rawFile) : null),
    enabled: compareEnabled && compareTarget !== null,
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

  const compareOptions = useMemo(() => {
    const optionMap = new Map<string, CompareTarget>();
    (precursors || []).forEach((p) => {
      const target = { id: p.precursor_id, rawFile: p.raw_file };
      optionMap.set(buildCompareKey(target), target);
    });
    if (selectedId !== null) {
      const selectedTarget = { id: selectedId, rawFile: selectedRawFile };
      optionMap.set(buildCompareKey(selectedTarget), selectedTarget);
    }
    if (compareTarget) {
      optionMap.set(buildCompareKey(compareTarget), compareTarget);
    }
    return Array.from(optionMap.values());
  }, [precursors, selectedId, selectedRawFile, compareTarget]);

  const compareDelta = useMemo(() => {
    if (!selectedPrecursor || !comparisonPrecursor) return null;
    const mzDelta = comparisonPrecursor.mz - selectedPrecursor.mz;
    const rtDeltaMin = (comparisonPrecursor.rt_seconds - selectedPrecursor.rt_seconds) / 60;
    const mobilityDelta = comparisonPrecursor.mobility - selectedPrecursor.mobility;
    const ppmDelta = selectedPrecursor.mz !== 0 ? (mzDelta / selectedPrecursor.mz) * 1_000_000 : 0;
    return {
      mzDelta,
      rtDeltaMin,
      mobilityDelta,
      ppmDelta,
      sameCharge: selectedPrecursor.charge === comparisonPrecursor.charge,
    };
  }, [selectedPrecursor, comparisonPrecursor]);

  const { data: compareStudies } = useQuery({
    queryKey: ['compare-studies'],
    queryFn: getStudies,
    enabled: activeDataset !== null,
  });

  const { data: activeDatasetInfo } = useQuery({
    queryKey: ['dataset-info', activeDataset],
    queryFn: () => (activeDataset ? getDatasetInfo(activeDataset) : null),
    enabled: activeDataset !== null,
  });

  const {
    data: datasetComparisonInfo,
    isLoading: isLoadingDatasetComparison,
  } = useQuery({
    queryKey: ['dataset-info-compare', datasetCompareTarget],
    queryFn: () => (datasetCompareTarget ? getDatasetInfo(datasetCompareTarget) : null),
    enabled: datasetCompareEnabled && datasetCompareTarget !== null,
  });

  const {
    data: datasetOverlap,
    isLoading: isLoadingDatasetOverlap,
  } = useQuery({
    queryKey: ['dataset-overlap', activeDataset, datasetCompareTarget],
    queryFn: () => (activeDataset && datasetCompareTarget ? getDatasetOverlap(activeDataset, datasetCompareTarget) : null),
    enabled: datasetCompareEnabled && activeDataset !== null && datasetCompareTarget !== null,
  });

  const datasetCompareOptions = useMemo(() => {
    const options = new Set<string>();
    (compareStudies || []).forEach((study: StudySummary) => {
      study.datasets.forEach((accession) => {
        if (accession && accession !== activeDataset) {
          options.add(accession);
        }
      });
    });
    return Array.from(options).sort();
  }, [compareStudies, activeDataset]);

  const datasetDelta = useMemo(() => {
    if (!activeDatasetInfo || !datasetComparisonInfo) return null;
    const nPrecursors = datasetComparisonInfo.n_precursors - activeDatasetInfo.n_precursors;
    const nAllThree = (datasetComparisonInfo.n_all_three ?? 0) - (activeDatasetInfo.n_all_three ?? 0);
    const nAtLeastTwo = (datasetComparisonInfo.n_at_least_two ?? 0) - (activeDatasetInfo.n_at_least_two ?? 0);
    const qualityA = activeDatasetInfo.quality?.pct_high_quality ?? null;
    const qualityB = datasetComparisonInfo.quality?.pct_high_quality ?? null;
    const rtA = activeDatasetInfo.quality?.rt_r2_median ?? null;
    const rtB = datasetComparisonInfo.quality?.rt_r2_median ?? null;
    const imA = activeDatasetInfo.quality?.im_r2_median ?? null;
    const imB = datasetComparisonInfo.quality?.im_r2_median ?? null;
    return {
      nPrecursors,
      nAllThree,
      nAtLeastTwo,
      qualityPct: qualityA !== null && qualityB !== null ? qualityB - qualityA : null,
      rtR2: rtA !== null && rtB !== null ? rtB - rtA : null,
      imR2: imA !== null && imB !== null ? imB - imA : null,
      sameStudy: activeDatasetInfo.study_id === datasetComparisonInfo.study_id,
    };
  }, [activeDatasetInfo, datasetComparisonInfo]);

  const isDefaultSingleView = useMemo(
    () =>
      areFiltersEqual(filters, DEFAULT_FILTERS) &&
      page === 0 &&
      clampTableWidth(tableWidth) === DEFAULT_TABLE_WIDTH &&
      !compareEnabled &&
      compareTarget === null &&
      !datasetCompareEnabled &&
      datasetCompareTarget === null &&
      overlapFocusMode === 'all',
    [
      filters,
      page,
      tableWidth,
      compareEnabled,
      compareTarget,
      datasetCompareEnabled,
      datasetCompareTarget,
      overlapFocusMode,
    ]
  );
  const viewSelectorValue = selectedViewId || (isDefaultSingleView ? '__default__' : '');

  return (
    <div className="app-shell h-screen flex flex-col">
      {/* Header */}
      <header className="chrome-header flex-none px-3 py-2 md:px-4 md:py-3">
        <div className="toolbar-lane">
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

          <div className="segmented-toggle ml-auto">
            <button
              type="button"
              onClick={() => setSubView('summary')}
              className={subView === 'summary' ? 'active' : ''}
            >
              Summary
            </button>
            <button
              type="button"
              onClick={() => setSubView('browse')}
              className={subView === 'browse' ? 'active' : ''}
            >
              Browse
            </button>
          </div>
        </div>

        {subView === 'browse' && (
          <div className="toolbar-stack mt-2 md:mt-3 reveal-up">
            <div className="toolbar-lane">
              <div className="toolbar-cluster toolbar-cluster--accent md:flex-1">
                <div className="toolbar-field">
                  <span className="toolbar-label">View</span>
                  <select
                    value={viewSelectorValue}
                    onChange={(e) => {
                      const id = e.target.value;
                      if (id === '__default__') {
                        handleResetToDefaultView();
                        return;
                      }
                      if (!id) {
                        setSelectedViewId('');
                        return;
                      }
                      handleApplySavedView(id);
                    }}
                    className="control-select max-w-[180px]"
                  >
                    <option value="__default__">Default (Single)</option>
                    <option value="">Custom</option>
                    {savedViews.map((view) => (
                      <option key={view.id} value={view.id}>
                        {view.name}
                      </option>
                    ))}
                  </select>
                </div>
                <button onClick={handleResetToDefaultView} className="btn-secondary" type="button">
                  Reset
                </button>
                <button onClick={handleSaveView} className="btn-secondary" type="button">
                  Save
                </button>
                <button
                  onClick={handleDeleteView}
                  disabled={!selectedViewId}
                  className="btn-secondary"
                  type="button"
                >
                  Delete
                </button>
                <button onClick={handleCopyLink} className="btn-secondary" type="button">
                  {linkCopyState === 'copied' ? 'Copied' : linkCopyState === 'failed' ? 'Retry' : 'Copy Link'}
                </button>
              </div>

              <div className="hidden md:flex toolbar-cluster">
                <label className="control-check flex items-center gap-2 text-sm cursor-pointer rounded px-1">
                  <input
                    type="checkbox"
                    checked={compareEnabled}
                    onChange={(e) => setCompareEnabled(e.target.checked)}
                    className="rounded border-slate-500/70 bg-slate-900/80"
                  />
                  <span>Precursor Compare</span>
                </label>
                <button
                  onClick={handleSetCurrentAsCompare}
                  disabled={selectedId === null}
                  className="btn-secondary"
                  type="button"
                >
                  Set Current as B
                </button>
                <select
                  value={buildCompareKey(compareTarget)}
                  onChange={(e) => setCompareTarget(parseCompareKey(e.target.value))}
                  disabled={!compareEnabled}
                  className="control-select max-w-[210px]"
                >
                  <option value="">B precursor</option>
                  {compareOptions.map((target) => (
                    <option key={buildCompareKey(target)} value={buildCompareKey(target)}>
                      #{target.id} {target.rawFile ? `(${target.rawFile.slice(-16)})` : ''}
                    </option>
                  ))}
                </select>
                <button
                  onClick={() => setCompareTarget(null)}
                  disabled={!compareEnabled || compareTarget === null}
                  className="btn-secondary"
                  type="button"
                >
                  Clear
                </button>
              </div>

              <div className="hidden md:flex toolbar-cluster">
                <label className="control-check flex items-center gap-2 text-sm cursor-pointer rounded px-1">
                  <input
                    type="checkbox"
                    checked={datasetCompareEnabled}
                    onChange={(e) => setDatasetCompareEnabled(e.target.checked)}
                    className="rounded border-slate-500/70 bg-slate-900/80"
                    disabled={!activeDataset}
                  />
                  <span>Dataset Compare</span>
                </label>
                <select
                  value={datasetCompareTarget ?? ''}
                  onChange={(e) => setDatasetCompareTarget(e.target.value || null)}
                  disabled={!datasetCompareEnabled || !activeDataset}
                  className="control-select max-w-[220px]"
                >
                  <option value="">B dataset</option>
                  {datasetCompareOptions.map((accession) => (
                    <option key={accession} value={accession}>
                      {accession}
                    </option>
                  ))}
                </select>
                <button
                  onClick={() => setDatasetCompareTarget(null)}
                  disabled={!datasetCompareEnabled || datasetCompareTarget === null}
                  className="btn-secondary"
                  type="button"
                >
                  Clear
                </button>
              </div>

              <button
                type="button"
                className="btn-secondary md:hidden ml-auto"
                onClick={() => setShowAdvancedControls((prev) => !prev)}
                aria-expanded={showAdvancedControls}
              >
                {showAdvancedControls ? 'Hide filters' : 'More filters'}
              </button>
            </div>

            <div className={`${showAdvancedControls ? 'block' : 'hidden'} md:block space-y-[0.55rem]`}>
              <div className="toolbar-lane md:hidden">
                <div className="toolbar-cluster">
                  <label className="control-check flex items-center gap-2 text-sm cursor-pointer rounded px-1">
                    <input
                      type="checkbox"
                      checked={compareEnabled}
                      onChange={(e) => setCompareEnabled(e.target.checked)}
                      className="rounded border-slate-500/70 bg-slate-900/80"
                    />
                    <span>Precursor Compare</span>
                  </label>
                  <button
                    onClick={handleSetCurrentAsCompare}
                    disabled={selectedId === null}
                    className="btn-secondary"
                    type="button"
                  >
                    Set Current as B
                  </button>
                  <select
                    value={buildCompareKey(compareTarget)}
                    onChange={(e) => setCompareTarget(parseCompareKey(e.target.value))}
                    disabled={!compareEnabled}
                    className="control-select max-w-[210px]"
                  >
                    <option value="">B precursor</option>
                    {compareOptions.map((target) => (
                      <option key={buildCompareKey(target)} value={buildCompareKey(target)}>
                        #{target.id} {target.rawFile ? `(${target.rawFile.slice(-16)})` : ''}
                      </option>
                    ))}
                  </select>
                  <button
                    onClick={() => setCompareTarget(null)}
                    disabled={!compareEnabled || compareTarget === null}
                    className="btn-secondary"
                    type="button"
                  >
                    Clear
                  </button>
                </div>

                <div className="toolbar-cluster">
                  <label className="control-check flex items-center gap-2 text-sm cursor-pointer rounded px-1">
                    <input
                      type="checkbox"
                      checked={datasetCompareEnabled}
                      onChange={(e) => setDatasetCompareEnabled(e.target.checked)}
                      className="rounded border-slate-500/70 bg-slate-900/80"
                      disabled={!activeDataset}
                    />
                    <span>Dataset Compare</span>
                  </label>
                  <select
                    value={datasetCompareTarget ?? ''}
                    onChange={(e) => setDatasetCompareTarget(e.target.value || null)}
                    disabled={!datasetCompareEnabled || !activeDataset}
                    className="control-select max-w-[220px]"
                  >
                    <option value="">B dataset</option>
                    {datasetCompareOptions.map((accession) => (
                      <option key={accession} value={accession}>
                        {accession}
                      </option>
                    ))}
                  </select>
                  <button
                    onClick={() => setDatasetCompareTarget(null)}
                    disabled={!datasetCompareEnabled || datasetCompareTarget === null}
                    className="btn-secondary"
                    type="button"
                  >
                    Clear
                  </button>
                </div>
              </div>

              <div className="toolbar-lane">
                <div className="toolbar-cluster md:flex-1">
                  <label className="toolbar-field text-sm">
                    <span className="toolbar-label">Engines</span>
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

                  <label className="toolbar-field text-sm">
                    <span className="toolbar-label">Charge</span>
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
                    <label className="toolbar-field text-sm">
                      <span className="toolbar-label">File</span>
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

                  <label className="toolbar-field text-sm">
                    <span className="toolbar-label">Sort</span>
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
                </div>

                <div className="toolbar-cluster">
                  <label className="control-check flex items-center gap-2 text-sm cursor-pointer rounded px-1">
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

                  <label className="control-check flex items-center gap-2 text-sm cursor-pointer rounded px-1">
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
              </div>
            </div>
          </div>
        )}
      </header>

      <div className="flex-1 min-h-0 flex flex-col">
        {subView === 'summary' && (
          <div className="flex-1 min-h-0">
            <DatasetSummaryPage
              activeDataset={activeDataset}
              onBrowse={() => setSubView('browse')}
            />
          </div>
        )}

        {subView === 'browse' && datasetCompareEnabled && (
          <div className="px-2 pb-0">
            <div className="chrome-panel p-3">
              <div className="flex flex-wrap items-center gap-2 mb-2">
                <span className="metric-pill">Dataset A/B Compare</span>
                {datasetDelta && (
                  <>
                    <span className="metric-pill">delta precursors <span className="metric-value mono">{formatSigned(datasetDelta.nPrecursors, 0)}</span></span>
                    <span className="metric-pill">delta 3-engine <span className="metric-value mono">{formatSigned(datasetDelta.nAllThree, 0)}</span></span>
                    <span className="metric-pill">delta 2+ engine <span className="metric-value mono">{formatSigned(datasetDelta.nAtLeastTwo, 0)}</span></span>
                    {datasetDelta.qualityPct !== null && (
                      <span className="metric-pill">delta HQ% <span className="metric-value mono">{formatSigned(datasetDelta.qualityPct, 2)}</span></span>
                    )}
                    <span className="metric-pill">
                      Study <span className={datasetDelta.sameStudy ? 'status-good' : 'status-warn'}>{datasetDelta.sameStudy ? 'match' : 'different'}</span>
                    </span>
                  </>
                )}
                {datasetOverlap && (
                  <>
                    <span className="metric-pill">shared precursors <span className="metric-value mono">{datasetOverlap.shared_precursors.toLocaleString()}</span></span>
                    <span className="metric-pill">shared peptides <span className="metric-value mono">{datasetOverlap.shared_peptides.toLocaleString()}</span></span>
                    <span className="metric-pill">precursor jaccard <span className="metric-value mono">{(datasetOverlap.precursor_jaccard * 100).toFixed(1)}%</span></span>
                    <span className="metric-pill">peptide jaccard <span className="metric-value mono">{(datasetOverlap.peptide_jaccard * 100).toFixed(1)}%</span></span>
                  </>
                )}
                {!datasetDelta && datasetCompareTarget && (
                  <span className="text-xs text-slate-300">
                    {isLoadingDatasetComparison ? 'Loading dataset comparison...' : 'Unable to compute dataset delta.'}
                  </span>
                )}
                {!datasetOverlap && datasetCompareTarget && (
                  <span className="text-xs text-slate-300">
                    {isLoadingDatasetOverlap ? 'Computing overlap...' : 'Overlap not available yet.'}
                  </span>
                )}
                {overlapFocusMode !== 'all' && (
                  <span className="metric-pill">
                    Table filter
                    <span className="metric-value">
                      {overlapFocusMode === 'shared' ? 'shared (A∩B)' : 'unique A (A-B)'}
                    </span>
                  </span>
                )}
                {!datasetCompareTarget && (
                  <span className="text-xs text-slate-300">Select B dataset to compare summary metrics.</span>
                )}
              </div>

              <div className="grid grid-cols-1 xl:grid-cols-2 gap-2">
                <div className="panel-inset">
                  <div className="panel-inset-head">
                    <span className="control-label">Dataset A</span>
                    <span className="mono text-xs text-slate-300">{activeDataset ?? '-'}</span>
                  </div>
                  <table className="data-table text-sm">
                    <tbody>
                      <tr><td className="px-3 py-2">Precursors</td><td className="px-3 py-2 text-right">{formatNullableInt(activeDatasetInfo?.n_precursors)}</td></tr>
                      <tr><td className="px-3 py-2">3-engine</td><td className="px-3 py-2 text-right">{formatNullableInt(activeDatasetInfo?.n_all_three)}</td></tr>
                      <tr><td className="px-3 py-2">2+ engine</td><td className="px-3 py-2 text-right">{formatNullableInt(activeDatasetInfo?.n_at_least_two)}</td></tr>
                      <tr><td className="px-3 py-2">High quality %</td><td className="px-3 py-2 text-right">{formatNullableNumber(activeDatasetInfo?.quality?.pct_high_quality, 2)}</td></tr>
                      <tr><td className="px-3 py-2">RT R² median</td><td className="px-3 py-2 text-right">{formatNullableNumber(activeDatasetInfo?.quality?.rt_r2_median, 3)}</td></tr>
                      <tr><td className="px-3 py-2">IM R² median</td><td className="px-3 py-2 text-right">{formatNullableNumber(activeDatasetInfo?.quality?.im_r2_median, 3)}</td></tr>
                    </tbody>
                  </table>
                </div>

                <div className="panel-inset">
                  <div className="panel-inset-head">
                    <span className="control-label">Dataset B</span>
                    <span className="mono text-xs text-slate-300">{datasetCompareTarget ?? '-'}</span>
                  </div>
                  <table className="data-table text-sm">
                    <tbody>
                      <tr><td className="px-3 py-2">Precursors</td><td className="px-3 py-2 text-right">{formatNullableInt(datasetComparisonInfo?.n_precursors)}</td></tr>
                      <tr><td className="px-3 py-2">3-engine</td><td className="px-3 py-2 text-right">{formatNullableInt(datasetComparisonInfo?.n_all_three)}</td></tr>
                      <tr><td className="px-3 py-2">2+ engine</td><td className="px-3 py-2 text-right">{formatNullableInt(datasetComparisonInfo?.n_at_least_two)}</td></tr>
                      <tr><td className="px-3 py-2">High quality %</td><td className="px-3 py-2 text-right">{formatNullableNumber(datasetComparisonInfo?.quality?.pct_high_quality, 2)}</td></tr>
                      <tr><td className="px-3 py-2">RT R² median</td><td className="px-3 py-2 text-right">{formatNullableNumber(datasetComparisonInfo?.quality?.rt_r2_median, 3)}</td></tr>
                      <tr><td className="px-3 py-2">IM R² median</td><td className="px-3 py-2 text-right">{formatNullableNumber(datasetComparisonInfo?.quality?.im_r2_median, 3)}</td></tr>
                    </tbody>
                  </table>
                </div>
              </div>

              {datasetOverlap && (
                <div className="grid grid-cols-1 xl:grid-cols-2 gap-2 mt-2">
                  <OverlapVisualCard
                    title="Precursor Overlap"
                    labelA="A"
                    labelB="B"
                    totalA={datasetOverlap.precursors_a}
                    totalB={datasetOverlap.precursors_b}
                    shared={datasetOverlap.shared_precursors}
                    uniqueA={datasetOverlap.unique_precursors_a}
                    uniqueB={datasetOverlap.unique_precursors_b}
                    jaccard={datasetOverlap.precursor_jaccard}
                    activeFocus={overlapFocusMode}
                    onFocusShared={() => handleOverlapFocusChange('shared')}
                    onFocusUniqueA={() => handleOverlapFocusChange('unique_a')}
                    onClearFocus={() => handleOverlapFocusChange('all')}
                  />

                  <OverlapVisualCard
                    title="Peptide Overlap"
                    labelA="A"
                    labelB="B"
                    totalA={datasetOverlap.peptides_a}
                    totalB={datasetOverlap.peptides_b}
                    shared={datasetOverlap.shared_peptides}
                    uniqueA={datasetOverlap.unique_peptides_a}
                    uniqueB={datasetOverlap.unique_peptides_b}
                    jaccard={datasetOverlap.peptide_jaccard}
                  />
                </div>
              )}
            </div>
          </div>
        )}

        {/* Main content */}
        {subView === 'browse' && <div id="main-content" className={`flex-1 flex min-h-0 p-2 gap-2 ${isResizing ? 'select-none' : ''}`}>
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
            {!compareEnabled && (
              <PrecursorViz precursor={selectedPrecursor || null} isLoading={isLoadingDetail} />
            )}

            {compareEnabled && (
              <div className="h-full flex flex-col min-h-0">
                <div className="flex-none px-3 py-2 border-b border-[#25415f] bg-slate-950/25 flex flex-wrap items-center gap-2">
                  <span className="metric-pill">A/B Compare Mode</span>
                  {compareDelta && (
                    <>
                      <span className="metric-pill">delta m/z <span className="metric-value mono">{formatSigned(compareDelta.mzDelta, 4)}</span></span>
                      <span className="metric-pill">delta ppm <span className="metric-value mono">{formatSigned(compareDelta.ppmDelta, 1)}</span></span>
                      <span className="metric-pill">delta RT (min) <span className="metric-value mono">{formatSigned(compareDelta.rtDeltaMin, 2)}</span></span>
                      <span className="metric-pill">delta 1/K0 <span className="metric-value mono">{formatSigned(compareDelta.mobilityDelta, 3)}</span></span>
                      <span className="metric-pill">
                        Charge <span className={compareDelta.sameCharge ? 'status-good' : 'status-warn'}>{compareDelta.sameCharge ? 'match' : 'mismatch'}</span>
                      </span>
                    </>
                  )}
                  {!compareDelta && (
                    <span className="text-xs text-slate-300">Select precursor B to calculate deltas.</span>
                  )}
                </div>

                <div className="flex-1 min-h-0 grid grid-cols-1 xl:grid-cols-2 gap-2 p-2 overflow-auto">
                  <div className="panel-inset flex flex-col min-h-[34rem] xl:min-h-0">
                    <div className="panel-inset-head">
                      <span className="control-label">A / Current Selection</span>
                      {selectedId !== null && (
                        <span className="mono text-xs text-slate-300">
                          #{selectedId}{selectedRawFile ? ` (${selectedRawFile.slice(-16)})` : ''}
                        </span>
                      )}
                    </div>
                    <div className="flex-1 min-h-0">
                      <PrecursorViz precursor={selectedPrecursor || null} isLoading={isLoadingDetail} />
                    </div>
                  </div>

                  <div className="panel-inset flex flex-col min-h-[34rem] xl:min-h-0">
                    <div className="panel-inset-head">
                      <span className="control-label">B / Comparison</span>
                      {compareTarget && (
                        <span className="mono text-xs text-slate-300">
                          #{compareTarget.id}{compareTarget.rawFile ? ` (${compareTarget.rawFile.slice(-16)})` : ''}
                        </span>
                      )}
                    </div>
                    <div className="flex-1 min-h-0">
                      {compareTarget ? (
                        <PrecursorViz precursor={comparisonPrecursor || null} isLoading={isLoadingComparison} />
                      ) : (
                        <div className="h-full flex items-center justify-center p-6">
                          <div className="metric-pill">Choose B precursor from the compare controls.</div>
                        </div>
                      )}
                    </div>
                  </div>
                </div>
              </div>
            )}
          </div>
        </div>}
      </div>
    </div>
  );
}

function Dashboard() {
  const [viewMode, setViewMode] = useState<ViewMode>('landing');
  const [activeDataset, setActiveDataset] = useState<string | null>(null);
  const [isCollectionMode, setIsCollectionMode] = useState(false);
  const queryClientHook = useQueryClient();

  // Check app status on mount
  const { data: appStatus, isLoading: isLoadingStatus } = useQuery({
    queryKey: ['appStatus'],
    queryFn: getAppStatus,
    retry: false,
  });

  // Update state when app status loads — check URL for direct page navigation
  useEffect(() => {
    if (appStatus) {
      setIsCollectionMode(appStatus.mode === 'collection');
      setActiveDataset(appStatus.active_dataset);

      const params = new URLSearchParams(window.location.search);
      const pageParam = params.get('page');

      if (pageParam === 'visit') {
        setViewMode('visit');
      } else if (pageParam === 'blueprint') {
        setViewMode('blueprint');
      } else if (pageParam === 'collection') {
        if (appStatus.mode === 'collection' && !appStatus.store_loaded) {
          setViewMode('collection');
        } else {
          setViewMode('dataset');
        }
      } else if (pageParam === 'dataset') {
        setViewMode('dataset');
      } else {
        // Default: show landing page
        setViewMode('landing');
      }
    }
  }, [appStatus]);

  // Sync viewMode to URL
  useEffect(() => {
    const params = new URLSearchParams(window.location.search);
    if (viewMode === 'landing') {
      params.delete('page');
    } else if (viewMode === 'visit') {
      params.set('page', 'visit');
    } else if (viewMode === 'blueprint') {
      params.set('page', 'blueprint');
    } else if (viewMode === 'collection') {
      params.set('page', 'collection');
    } else {
      params.set('page', 'dataset');
    }
    const nextQuery = params.toString();
    const nextUrl = `${window.location.pathname}${nextQuery ? `?${nextQuery}` : ''}`;
    window.history.replaceState(null, '', nextUrl);
  }, [viewMode]);

  const handleDatasetLoaded = (accession: string) => {
    setActiveDataset(accession);
    setViewMode('dataset');
    // Invalidate queries to refresh with new data
    queryClientHook.invalidateQueries({ queryKey: ['precursors'] });
    queryClientHook.invalidateQueries({ queryKey: ['stats'] });
    queryClientHook.invalidateQueries({ queryKey: ['raw_files'] });
    queryClientHook.invalidateQueries({ queryKey: ['dataset-summary'] });
  };

  const handleNavigateToCollection = () => {
    setViewMode('collection');
  };

  const handleNavigateToLanding = () => setViewMode('landing');
  const handleNavigateToVisit = () => setViewMode('visit');
  const handleNavigateToBlueprint = () => setViewMode('blueprint');

  const handleExploreFromPublic = () => {
    if (isCollectionMode) setViewMode('collection');
    else setViewMode('dataset');
  };

  if (isLoadingStatus) {
    return (
      <div className="app-shell h-screen flex items-center justify-center">
        <div className="metric-pill">Loading dashboard status...</div>
      </div>
    );
  }

  // Landing page
  if (viewMode === 'landing') {
    return (
      <LandingPage
        onViewSummary={handleNavigateToVisit}
        onNavigateBlueprint={handleNavigateToBlueprint}
        explorerUrl={EXPLORER_URL || undefined}
        onExploreData={!EXPLORER_URL ? handleExploreFromPublic : undefined}
      />
    );
  }

  // Visit summary page
  if (viewMode === 'visit') {
    return (
      <VisitSummaryPage
        onBack={handleNavigateToLanding}
        onNavigateBlueprint={handleNavigateToBlueprint}
        explorerUrl={EXPLORER_URL || undefined}
        onExploreData={!EXPLORER_URL ? handleExploreFromPublic : undefined}
      />
    );
  }

  // Blueprint documentation page
  if (viewMode === 'blueprint') {
    return (
      <BlueprintPage
        onBack={handleNavigateToLanding}
        onNavigateVisit={handleNavigateToVisit}
        explorerUrl={EXPLORER_URL || undefined}
        onExploreData={!EXPLORER_URL ? handleExploreFromPublic : undefined}
      />
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
