import { useQuery } from '@tanstack/react-query';
import { getCollectionInfo, getStudies, getStats } from '../api';
import type { StudySummary } from '../api';
import SiteNav from './SiteNav';
import PipelineOverview from './PipelineOverview';

interface VisitSummaryPageProps {
  onBack: () => void;
  explorerUrl?: string;
  onExploreData?: () => void;
}

const techStack = [
  { name: 'Snakemake', role: 'Orchestration' },
  { name: 'rustims / imspy', role: 'Raw processing' },
  { name: 'FragPipe', role: 'Spectrum-centric search' },
  { name: 'DIA-NN', role: 'Peptide-centric search' },
  { name: 'Sage', role: 'Fast Rust search engine' },
  { name: 'DuckDB / Parquet', role: 'Columnar storage' },
  { name: 'TensorFlow', role: 'Model training' },
  { name: 'Singularity / Docker', role: 'Containers' },
];

const roadmap = [
  { label: '10 HeLa Datasets', detail: 'Initial QC and pipeline validation', status: 'completed' as const },
  { label: 'Raw Feature Extraction', detail: '4D signal extraction via rustims', status: 'completed' as const },
  { label: 'Containerized Workers', detail: 'HPC deployment on Mogon2 SLURM', status: 'current' as const },
  { label: '100 Datasets / v1.0 Snapshot', detail: 'First citable versioned release', status: 'future' as const },
  { label: 'All timsTOF on PRIDE', detail: 'Comprehensive reference layer', status: 'future' as const },
];

function formatNumber(n: number): string {
  if (n >= 1_000_000) return `${(n / 1_000_000).toFixed(1)}M`;
  if (n >= 1_000) return `${(n / 1_000).toFixed(0)}K`;
  return n.toLocaleString();
}

export default function VisitSummaryPage({
  onBack,
  explorerUrl,
  onExploreData,
}: VisitSummaryPageProps) {
  const { data: collectionInfo } = useQuery({
    queryKey: ['collection-info-visit'],
    queryFn: getCollectionInfo,
    retry: false,
    staleTime: 60_000,
  });

  const { data: studies } = useQuery({
    queryKey: ['studies-visit'],
    queryFn: getStudies,
    retry: false,
    staleTime: 60_000,
  });

  const { data: stats } = useQuery({
    queryKey: ['stats-visit'],
    queryFn: getStats,
    retry: false,
    enabled: !collectionInfo,
    staleTime: 60_000,
  });

  const totalPrecursors = collectionInfo?.n_total_precursors ?? stats?.total_precursors;
  const totalDatasets = collectionInfo?.n_datasets;
  const totalStudies = collectionInfo?.n_studies;
  const uniqueOrganisms = studies
    ? new Set(studies.map((s: StudySummary) => s.organism).filter(Boolean)).size
    : undefined;

  return (
    <div className="app-shell min-h-screen flex flex-col">
      <SiteNav
        currentPage="visit"
        onNavigateLanding={onBack}
        onNavigateVisit={() => {}}
        explorerUrl={explorerUrl}
        onExploreData={onExploreData}
      />

      <main className="flex-1 overflow-y-auto">
        {/* Header */}
        <section className="px-6 pt-12 pb-6 text-center max-w-4xl mx-auto">
          <h1
            className="text-3xl md:text-4xl font-bold tracking-tight mb-2"
            style={{ color: 'var(--accent-1)' }}
          >
            Project Summary
          </h1>
          <p className="text-sm md:text-base" style={{ color: 'var(--text-3)' }}>
            San Jose systematically reprocesses public timsTOF data from PRIDE through triple
            orthogonal validation, extracts raw 4D signal features, and builds versioned
            snapshots for model training.
          </p>
        </section>

        {/* Live Stats Cards */}
        <section className="px-6 py-8 max-w-4xl mx-auto">
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
            <StatCard label="Precursors" value={totalPrecursors} />
            <StatCard label="Datasets" value={totalDatasets} />
            <StatCard label="Studies" value={totalStudies} />
            <StatCard label="Organisms" value={uniqueOrganisms} />
          </div>
        </section>

        {/* Pipeline Overview */}
        <section className="px-6 py-8 max-w-5xl mx-auto">
          <h2 className="subtle-label mb-2 text-center">Processing Pipeline</h2>
          <p className="text-sm text-center mb-6" style={{ color: 'var(--text-3)' }}>
            Each PRIDE accession is processed through a 6-step runner with checkpointing.
          </p>
          <PipelineOverview compact />
        </section>

        {/* Studies Overview */}
        {studies && studies.length > 0 && (
          <section className="px-6 py-8 max-w-5xl mx-auto">
            <h2 className="subtle-label mb-4 text-center">Processed Studies</h2>
            <div className="panel-inset">
              <table className="data-table">
                <thead>
                  <tr>
                    <th className="text-left px-4 py-2">Study</th>
                    <th className="text-left px-4 py-2">Organism</th>
                    <th className="text-right px-4 py-2">Datasets</th>
                    <th className="text-right px-4 py-2">Precursors</th>
                  </tr>
                </thead>
                <tbody>
                  {studies.map((study: StudySummary) => (
                    <tr key={study.id}>
                      <td className="px-4 py-2 text-sm">
                        <span style={{ color: 'var(--text-1)' }}>{study.title || study.id}</span>
                        {study.description && (
                          <span className="block text-xs" style={{ color: 'var(--text-3)' }}>
                            {study.description}
                          </span>
                        )}
                      </td>
                      <td className="px-4 py-2 text-sm" style={{ color: 'var(--text-2)' }}>
                        {study.organism || '\u2014'}
                      </td>
                      <td className="px-4 py-2 text-sm text-right mono" style={{ color: 'var(--accent-2)' }}>
                        {study.n_datasets}
                      </td>
                      <td className="px-4 py-2 text-sm text-right mono" style={{ color: 'var(--accent-2)' }}>
                        {formatNumber(study.n_total_precursors)}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </section>
        )}

        {/* Technology Stack */}
        <section className="px-6 py-8 max-w-4xl mx-auto">
          <h2 className="subtle-label mb-4 text-center">Technology Stack</h2>
          <div className="flex flex-wrap justify-center gap-3">
            {techStack.map((tech) => (
              <div key={tech.name} className="tech-badge">
                <span className="font-semibold" style={{ color: 'var(--text-1)' }}>
                  {tech.name}
                </span>
                <span className="text-xs" style={{ color: 'var(--text-3)' }}>
                  {tech.role}
                </span>
              </div>
            ))}
          </div>
        </section>

        {/* Roadmap */}
        <section className="px-6 py-8 max-w-3xl mx-auto">
          <h2 className="subtle-label mb-6 text-center">Roadmap</h2>
          <div className="flex flex-col">
            {roadmap.map((item, i) => (
              <div key={item.label} className="timeline-item">
                <div className={`timeline-dot ${item.status}`} />
                <div className={`pb-${i < roadmap.length - 1 ? '6' : '0'}`}>
                  <h3
                    className="text-sm font-semibold"
                    style={{
                      color:
                        item.status === 'completed'
                          ? 'var(--success-1)'
                          : item.status === 'current'
                          ? 'var(--accent-1)'
                          : 'var(--text-3)',
                    }}
                  >
                    {item.label}
                  </h3>
                  <p className="text-xs" style={{ color: 'var(--text-3)' }}>
                    {item.detail}
                  </p>
                </div>
              </div>
            ))}
          </div>
        </section>

        {/* Institutional Context */}
        <section className="px-6 py-8 max-w-3xl mx-auto text-center">
          <h2 className="subtle-label mb-4">About the Name</h2>
          <p className="text-sm leading-relaxed mb-3" style={{ color: 'var(--text-2)' }}>
            Named after the <em>San Jose</em>, a sunken ship whose rediscovery revealed immense,
            carefully preserved value beneath the surface. PRIDE is a vast repository where the
            real value is hidden in raw experimental data, waiting to be systematically recovered,
            catalogued, and understood.
          </p>
          <p className="text-xs" style={{ color: 'var(--text-3)' }}>
            Johannes Gutenberg University Mainz &middot; HPC: Mogon2 SLURM Cluster
          </p>
        </section>

        {/* CTA Footer */}
        <footer
          className="px-6 py-10 text-center"
          style={{ borderTop: '1px solid var(--border-2)' }}
        >
          {(explorerUrl || onExploreData) && (
            <button
              className="btn-primary"
              style={{ padding: '0.6rem 1.5rem', fontSize: '0.9rem' }}
              onClick={() => {
                if (explorerUrl) window.open(explorerUrl, '_blank', 'noopener');
                else onExploreData?.();
              }}
            >
              Explore the Data
            </button>
          )}
        </footer>
      </main>
    </div>
  );
}

function StatCard({ label, value }: { label: string; value: number | undefined }) {
  return (
    <div className="chrome-panel stat-card-large">
      <div className="stat-number mono">
        {value != null ? formatNumber(value) : '\u2014'}
      </div>
      <div className="subtle-label mt-1">{label}</div>
    </div>
  );
}
