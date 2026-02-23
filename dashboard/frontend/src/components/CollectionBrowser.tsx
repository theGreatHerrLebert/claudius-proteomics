import { useState, useEffect } from 'react';
import {
  getStudies,
  getStudyDatasets,
  loadDataset,
  getCollectionInfo,
} from '../api';
import type { StudySummary, DatasetSummary, CollectionInfo } from '../api';

interface CollectionBrowserProps {
  onDatasetLoaded: (accession: string) => void;
  activeDataset?: string | null;
  onReturnToDataset?: () => void;
}

function formatNumber(n: number): string {
  if (n >= 1000000) return `${(n / 1000000).toFixed(1)}M`;
  if (n >= 1000) return `${(n / 1000).toFixed(0)}K`;
  return n.toString();
}

function QualityBar({ pct }: { pct: number | null }) {
  const p = pct || 0;
  const bars = Math.round(p / 10);
  return (
    <div className="flex items-center gap-1">
      <div className="flex gap-0.5">
        {Array.from({ length: 10 }, (_, i) => (
          <div
            key={i}
            className={`w-1.5 h-3 rounded-sm ${
              i < bars ? 'bg-emerald-300' : 'bg-slate-700/80'
            }`}
          />
        ))}
      </div>
      <span className="text-xs text-slate-300">{p.toFixed(1)}%</span>
    </div>
  );
}

function StudyCard({
  study,
  onClick,
}: {
  study: StudySummary;
  onClick: () => void;
}) {
  return (
    <button
      onClick={onClick}
      className="chrome-panel reveal-up p-4 text-left w-full transition-all duration-150 hover:border-cyan-400/55 hover:-translate-y-0.5 hover:shadow-[0_16px_30px_#0410218a]"
    >
      <div className="flex items-start justify-between mb-2">
        <h3 className="font-semibold text-slate-100 truncate flex-1">{study.title}</h3>
      </div>
      {study.organism && (
        <p className="text-sm text-cyan-200/90 italic mb-2">{study.organism}</p>
      )}
      {study.description && (
        <p className="text-sm text-slate-400 mb-3 line-clamp-2">{study.description}</p>
      )}
      <div className="flex items-center gap-2 text-sm text-slate-300">
        <span className="metric-pill">
          {study.n_datasets} dataset{study.n_datasets !== 1 ? 's' : ''}
        </span>
        <span className="metric-pill">
          {formatNumber(study.n_total_precursors)} precursors
        </span>
      </div>
    </button>
  );
}

function DatasetRow({
  dataset,
  onLoad,
  isLoading,
}: {
  dataset: DatasetSummary;
  onLoad: () => void;
  isLoading: boolean;
}) {
  return (
    <tr className="border-b border-slate-700/45 hover:bg-slate-700/15 transition-colors">
      <td className="px-4 py-3">
        <span className="mono text-cyan-200">{dataset.accession}</span>
        <span className="text-slate-400 text-sm ml-2">v{dataset.version}</span>
      </td>
      <td className="px-4 py-3 text-right">
        {formatNumber(dataset.n_precursors)}
      </td>
      <td className="px-4 py-3 text-right">
        <span className="text-emerald-300">
          {dataset.n_all_three !== null ? formatNumber(dataset.n_all_three) : '-'}
        </span>
      </td>
      <td className="px-4 py-3 text-right">
        <span className="text-cyan-300">
          {dataset.n_at_least_two !== null ? formatNumber(dataset.n_at_least_two) : '-'}
        </span>
      </td>
      <td className="px-4 py-3">
        <QualityBar pct={dataset.quality?.pct_high_quality ?? null} />
      </td>
      <td className="px-4 py-3 text-right">
        <button
          onClick={onLoad}
          disabled={isLoading}
          className="btn-primary"
        >
          {isLoading ? 'Loading...' : 'Open'}
        </button>
      </td>
    </tr>
  );
}

export default function CollectionBrowser({ onDatasetLoaded, activeDataset, onReturnToDataset }: CollectionBrowserProps) {
  const [collectionInfo, setCollectionInfo] = useState<CollectionInfo | null>(null);
  const [studies, setStudies] = useState<StudySummary[]>([]);
  const [selectedStudy, setSelectedStudy] = useState<StudySummary | null>(null);
  const [datasets, setDatasets] = useState<DatasetSummary[]>([]);
  const [loading, setLoading] = useState(true);
  const [loadingDataset, setLoadingDataset] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  // Load initial data
  useEffect(() => {
    async function loadData() {
      try {
        setLoading(true);
        setError(null);
        const [info, studyList] = await Promise.all([
          getCollectionInfo(),
          getStudies(),
        ]);
        setCollectionInfo(info);
        setStudies(studyList);
      } catch (e) {
        setError(e instanceof Error ? e.message : 'Failed to load collection');
      } finally {
        setLoading(false);
      }
    }
    loadData();
  }, []);

  // Load datasets when study is selected
  useEffect(() => {
    async function loadStudyDatasets() {
      if (!selectedStudy) {
        setDatasets([]);
        return;
      }
      try {
        setLoading(true);
        const datasetList = await getStudyDatasets(selectedStudy.id);
        setDatasets(datasetList);
      } catch (e) {
        setError(e instanceof Error ? e.message : 'Failed to load datasets');
      } finally {
        setLoading(false);
      }
    }
    loadStudyDatasets();
  }, [selectedStudy]);

  const handleLoadDataset = async (accession: string) => {
    try {
      setLoadingDataset(accession);
      await loadDataset(accession);
      onDatasetLoaded(accession);
    } catch (e) {
      setError(e instanceof Error ? e.message : 'Failed to load dataset');
    } finally {
      setLoadingDataset(null);
    }
  };

  if (loading && !collectionInfo) {
    return (
      <div className="flex items-center justify-center h-64">
        <div className="metric-pill">Loading collection...</div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex items-center justify-center h-64">
        <div className="metric-pill border-rose-400/35 text-rose-200">Error: {error}</div>
      </div>
    );
  }

  return (
    <div className="p-4 md:p-6 max-w-[1240px] mx-auto reveal-up">
      {/* Active dataset banner */}
      {activeDataset && onReturnToDataset && (
        <button
          onClick={onReturnToDataset}
          className="chrome-panel mb-4 w-full px-4 py-3 flex items-center justify-between hover:border-cyan-400/55 transition-colors"
        >
          <span className="text-sm text-slate-200">
            Currently viewing: <span className="mono text-cyan-200 font-medium">{activeDataset}</span>
          </span>
          <span className="text-sm text-cyan-300 flex items-center gap-1">
            Return to dataset <span>&rarr;</span>
          </span>
        </button>
      )}

      {/* Header */}
      <div className="mb-6">
        <div className="flex items-center gap-4 mb-2">
          {selectedStudy && (
            <button
              onClick={() => setSelectedStudy(null)}
              className="btn-secondary flex items-center gap-1"
            >
              <span>&larr;</span>
              <span>Back</span>
            </button>
          )}
          <h1 className="text-2xl md:text-3xl font-semibold text-slate-100 tracking-tight">
            {selectedStudy ? selectedStudy.title : 'San Jose Collection'}
          </h1>
        </div>
        {!selectedStudy && collectionInfo && (
          <div className="flex flex-wrap items-center gap-2">
            <span className="metric-pill">{collectionInfo.n_studies} studies</span>
            <span className="metric-pill">{collectionInfo.n_datasets} datasets</span>
            <span className="metric-pill">{formatNumber(collectionInfo.n_total_precursors)} precursors</span>
          </div>
        )}
        {selectedStudy && (
          <>
            {selectedStudy.organism && (
              <p className="text-slate-300 italic">{selectedStudy.organism}</p>
            )}
            {selectedStudy.publication && (
              <a
                href={selectedStudy.publication}
                target="_blank"
                rel="noopener noreferrer"
                className="text-cyan-300 hover:text-cyan-100 underline decoration-cyan-400/40 text-sm"
              >
                View Publication
              </a>
            )}
          </>
        )}
      </div>

      {/* Content */}
      {!selectedStudy ? (
        // Studies grid
        <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-3 gap-4">
          {studies.map((study) => (
            <StudyCard
              key={study.id}
              study={study}
              onClick={() => setSelectedStudy(study)}
            />
          ))}
          {studies.length === 0 && (
            <div className="col-span-full text-center text-slate-400 py-12 chrome-panel">
              No studies found in this collection.
            </div>
          )}
        </div>
      ) : (
        // Datasets table
        <div className="chrome-panel overflow-hidden">
          <table className="data-table">
            <thead className="text-sm">
              <tr>
                <th className="px-4 py-3 text-left">Accession</th>
                <th className="px-4 py-3 text-right">Precursors</th>
                <th className="px-4 py-3 text-right">3-Engine</th>
                <th className="px-4 py-3 text-right">2+ Engine</th>
                <th className="px-4 py-3 text-left">Quality</th>
                <th className="px-4 py-3 text-right">Action</th>
              </tr>
            </thead>
            <tbody className="text-slate-200">
              {datasets.map((dataset) => (
                <DatasetRow
                  key={dataset.accession}
                  dataset={dataset}
                  onLoad={() => handleLoadDataset(dataset.accession)}
                  isLoading={loadingDataset === dataset.accession}
                />
              ))}
              {datasets.length === 0 && !loading && (
                <tr>
                  <td colSpan={6} className="px-4 py-12 text-center text-slate-400">
                    No datasets found in this study.
                  </td>
                </tr>
              )}
            </tbody>
          </table>
        </div>
      )}
    </div>
  );
}
