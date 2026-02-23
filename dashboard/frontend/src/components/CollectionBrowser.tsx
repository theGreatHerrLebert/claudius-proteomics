import { useState, useEffect } from 'react';
import {
  getStudies,
  getStudyDatasets,
  loadDataset,
  getCollectionInfo,
  type StudySummary,
  type DatasetSummary,
  type CollectionInfo,
} from '../api';

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
              i < bars ? 'bg-green-500' : 'bg-gray-700'
            }`}
          />
        ))}
      </div>
      <span className="text-xs text-gray-400">{p.toFixed(1)}%</span>
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
      className="bg-gray-800 border border-gray-700 rounded-lg p-4 text-left hover:border-blue-500 hover:bg-gray-750 transition-colors w-full"
    >
      <div className="flex items-start justify-between mb-2">
        <h3 className="font-semibold text-white truncate flex-1">{study.title}</h3>
      </div>
      {study.organism && (
        <p className="text-sm text-gray-400 italic mb-2">{study.organism}</p>
      )}
      {study.description && (
        <p className="text-sm text-gray-500 mb-3 line-clamp-2">{study.description}</p>
      )}
      <div className="flex items-center gap-4 text-sm text-gray-400">
        <span>{study.n_datasets} dataset{study.n_datasets !== 1 ? 's' : ''}</span>
        <span>{formatNumber(study.n_total_precursors)} precursors</span>
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
    <tr className="border-b border-gray-700 hover:bg-gray-800">
      <td className="px-4 py-3">
        <span className="font-mono text-blue-400">{dataset.accession}</span>
        <span className="text-gray-500 text-sm ml-2">v{dataset.version}</span>
      </td>
      <td className="px-4 py-3 text-right">
        {formatNumber(dataset.n_precursors)}
      </td>
      <td className="px-4 py-3 text-right">
        <span className="text-green-400">
          {dataset.n_all_three !== null ? formatNumber(dataset.n_all_three) : '-'}
        </span>
      </td>
      <td className="px-4 py-3 text-right">
        <span className="text-blue-400">
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
          className="px-3 py-1 bg-blue-600 hover:bg-blue-500 disabled:bg-gray-600 text-white text-sm rounded transition-colors"
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
        <div className="text-gray-400">Loading collection...</div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex items-center justify-center h-64">
        <div className="text-red-400">Error: {error}</div>
      </div>
    );
  }

  return (
    <div className="p-6 max-w-6xl mx-auto">
      {/* Active dataset banner */}
      {activeDataset && onReturnToDataset && (
        <button
          onClick={onReturnToDataset}
          className="mb-4 w-full bg-blue-900/40 border border-blue-700 rounded-lg px-4 py-3 flex items-center justify-between hover:bg-blue-900/60 transition-colors"
        >
          <span className="text-sm text-gray-300">
            Currently viewing: <span className="font-mono text-blue-400 font-medium">{activeDataset}</span>
          </span>
          <span className="text-sm text-blue-400 flex items-center gap-1">
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
              className="text-blue-400 hover:text-blue-300 flex items-center gap-1"
            >
              <span>&larr;</span>
              <span>Back</span>
            </button>
          )}
          <h1 className="text-2xl font-bold text-white">
            {selectedStudy ? selectedStudy.title : 'San Jose Collection'}
          </h1>
        </div>
        {!selectedStudy && collectionInfo && (
          <p className="text-gray-400">
            {collectionInfo.n_studies} studies, {collectionInfo.n_datasets} datasets,{' '}
            {formatNumber(collectionInfo.n_total_precursors)} total precursors
          </p>
        )}
        {selectedStudy && (
          <>
            {selectedStudy.organism && (
              <p className="text-gray-400 italic">{selectedStudy.organism}</p>
            )}
            {selectedStudy.publication && (
              <a
                href={selectedStudy.publication}
                target="_blank"
                rel="noopener noreferrer"
                className="text-blue-400 hover:underline text-sm"
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
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
          {studies.map((study) => (
            <StudyCard
              key={study.id}
              study={study}
              onClick={() => setSelectedStudy(study)}
            />
          ))}
          {studies.length === 0 && (
            <div className="col-span-full text-center text-gray-500 py-12">
              No studies found in this collection.
            </div>
          )}
        </div>
      ) : (
        // Datasets table
        <div className="bg-gray-800 rounded-lg overflow-hidden border border-gray-700">
          <table className="w-full">
            <thead className="bg-gray-900 text-gray-400 text-sm">
              <tr>
                <th className="px-4 py-3 text-left">Accession</th>
                <th className="px-4 py-3 text-right">Precursors</th>
                <th className="px-4 py-3 text-right">3-Engine</th>
                <th className="px-4 py-3 text-right">2+ Engine</th>
                <th className="px-4 py-3 text-left">Quality</th>
                <th className="px-4 py-3 text-right">Action</th>
              </tr>
            </thead>
            <tbody className="text-gray-200">
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
                  <td colSpan={6} className="px-4 py-12 text-center text-gray-500">
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
