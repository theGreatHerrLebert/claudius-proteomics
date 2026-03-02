import { useQuery } from '@tanstack/react-query';
import { getCollectionInfo, getStudies, getStats } from '../api';
import type { StudySummary } from '../api';
import SiteNav from './SiteNav';
import PipelineOverview from './PipelineOverview';

interface VisitSummaryPageProps {
  onBack: () => void;
  onNavigateBlueprint?: () => void;
  explorerUrl?: string;
  onExploreData?: () => void;
}

const techStack = [
  { name: 'Snakemake', role: 'Orchestration', url: 'https://github.com/snakemake/snakemake' },
  { name: 'rustims / imspy', role: 'Raw processing', url: 'https://github.com/theGreatHerrLebert/rustims' },
  { name: 'FragPipe', role: 'Spectrum-centric search', url: 'https://github.com/Nesvilab/FragPipe' },
  { name: 'DIA-NN', role: 'Peptide-centric search', url: 'https://github.com/vdemichev/DiaNN' },
  { name: 'Sage', role: 'Fast Rust search engine', url: 'https://github.com/lazear/sage' },
  { name: 'DuckDB / Parquet', role: 'Columnar storage', url: 'https://github.com/duckdb/duckdb' },
  { name: 'PyTorch', role: 'Model training', url: 'https://github.com/pytorch/pytorch' },
  { name: 'Singularity / Docker', role: 'Containers', url: 'https://github.com/apptainer/apptainer' },
];

const roadmap = [
  { label: '10 HeLa Datasets', detail: 'Initial QC and pipeline validation', status: 'completed' as const },
  { label: 'Raw Feature Extraction', detail: '4D signal extraction via rustims', status: 'completed' as const },
  { label: 'Containerized Workers', detail: 'HPC deployment on Mogon2 SLURM', status: 'current' as const },
  { label: '100 Datasets / v1.0 Snapshot', detail: 'First citable versioned release', status: 'future' as const },
  { label: 'All timsTOF on PRIDE', detail: 'Comprehensive reference layer', status: 'future' as const },
];

const summaryHighlights = [
  'Triple-engine orthogonal validation',
  'Bias-aware dataset stratification',
  'Versioned, citable release snapshots',
];

function formatNumber(n: number): string {
  if (n >= 1_000_000) return `${(n / 1_000_000).toFixed(1)}M`;
  if (n >= 1_000) return `${(n / 1_000).toFixed(0)}K`;
  return n.toLocaleString();
}

function roadmapTextColor(status: 'completed' | 'current' | 'future'): string {
  if (status === 'completed') return 'var(--success-1)';
  if (status === 'current') return 'var(--accent-1)';
  return 'var(--text-3)';
}

export default function VisitSummaryPage({
  onBack,
  onNavigateBlueprint,
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
  const hasExploreAction = Boolean(explorerUrl || onExploreData);

  const handleExplore = () => {
    if (explorerUrl) window.open(explorerUrl, '_blank', 'noopener');
    else onExploreData?.();
  };

  const metricCards = [
    { label: 'Precursors', value: totalPrecursors },
    { label: 'Datasets', value: totalDatasets },
    { label: 'Studies', value: totalStudies },
    { label: 'Organisms', value: uniqueOrganisms },
  ];

  return (
    <div className="app-shell min-h-screen flex flex-col visit-page">
      <SiteNav
        currentPage="visit"
        onNavigateLanding={onBack}
        onNavigateVisit={() => {}}
        onNavigateBlueprint={onNavigateBlueprint}
        explorerUrl={explorerUrl}
        onExploreData={onExploreData}
      />

      <main className="flex-1 overflow-y-auto">
        <section className="visit-hero px-6 pt-14 pb-12 md:pt-20 md:pb-16">
          <div className="max-w-6xl mx-auto visit-hero-grid reveal-up">
            <div className="visit-hero-content">
              <p className="subtle-label mb-3">Project Summary</p>
              <h1 className="visit-hero-title mb-4">San José in Context</h1>
              <p className="text-sm md:text-base mb-7" style={{ color: 'var(--text-2)' }}>
                San José systematically reprocesses public timsTOF data from PRIDE through triple
                orthogonal validation, extracts raw 4D signal features, and builds versioned
                snapshots designed for reliable, reproducible model training.
              </p>

              <div className="visit-hero-actions">
                {hasExploreAction && (
                  <button className="btn-primary landing-cta-primary" onClick={handleExplore}>
                    Explore the Data
                  </button>
                )}
                <button className="btn-secondary landing-cta-secondary" onClick={onBack}>
                  Back to Home
                </button>
              </div>

              <div className="flex gap-2 mt-5 flex-wrap">
                {summaryHighlights.map((highlight) => (
                  <span key={highlight} className="metric-pill">
                    {highlight}
                  </span>
                ))}
              </div>
            </div>

            <aside className="chrome-panel visit-overview-panel p-5 md:p-6">
              <p className="subtle-label mb-3">Live Program Metrics</p>
              <div className="visit-kpi-grid mb-4">
                {metricCards.map((metric) => (
                  <KpiCard key={metric.label} label={metric.label} value={metric.value} />
                ))}
              </div>
              <div className="landing-divider" />
              <p className="text-sm mb-2" style={{ color: 'var(--text-2)' }}>
                End-to-end pipeline coverage
              </p>
              <ul className="landing-checklist">
                <li>Automated checkpointing across the full runner</li>
                <li>Release workflows reviewed through quality gates</li>
                <li>Traceable provenance from PRIDE to final snapshot</li>
              </ul>
            </aside>
          </div>
        </section>

        <section className="px-6 py-10 max-w-6xl mx-auto">
          <div className="grid grid-cols-1 xl:grid-cols-[1.45fr,1fr] gap-4">
            <div className="chrome-panel p-5 md:p-6">
              <h2 className="subtle-label mb-2 text-center">Processing Pipeline</h2>
              <p className="text-sm text-center mb-6" style={{ color: 'var(--text-3)' }}>
                Each PRIDE accession is processed through a 6-step runner with checkpointing.
              </p>
              <PipelineOverview compact />
            </div>

            <div className="chrome-panel p-5 md:p-6 visit-roadmap-panel">
              <h2 className="subtle-label mb-4">Roadmap</h2>
              <div className="flex flex-col">
                {roadmap.map((item) => (
                  <div key={item.label} className="timeline-item">
                    <div className={`timeline-dot ${item.status}`} />
                    <div>
                      <h3 className="text-sm font-semibold" style={{ color: roadmapTextColor(item.status) }}>
                        {item.label}
                      </h3>
                      <p className="text-xs" style={{ color: 'var(--text-3)' }}>
                        {item.detail}
                      </p>
                    </div>
                  </div>
                ))}
              </div>
            </div>
          </div>
        </section>

        {studies && studies.length > 0 && (
          <section className="px-6 py-10 max-w-6xl mx-auto">
            <div className="flex items-center justify-between gap-3 mb-4 flex-wrap">
              <h2 className="subtle-label">Processed Studies</h2>
              <span className="metric-pill">
                <span className="metric-value mono">{studies.length}</span> studies in this snapshot
              </span>
            </div>
            <div className="panel-inset visit-studies-wrap">
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
                          <span className="block text-xs mt-0.5" style={{ color: 'var(--text-3)' }}>
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

        <section className="px-6 py-10 max-w-6xl mx-auto">
          <div className="grid grid-cols-1 xl:grid-cols-[1.2fr,1fr] gap-4">
            <div className="chrome-panel p-5 md:p-6">
              <h2 className="subtle-label mb-4">Technology Stack</h2>
              <div className="visit-tech-grid">
                {techStack.map((tech) => (
                  <a
                    key={tech.name}
                    href={tech.url}
                    target="_blank"
                    rel="noopener noreferrer"
                    className="visit-tech-card"
                    style={{ textDecoration: 'none' }}
                  >
                    <p className="text-sm font-semibold" style={{ color: 'var(--text-1)' }}>
                      {tech.name}
                    </p>
                    <p className="text-xs mt-0.5" style={{ color: 'var(--text-3)' }}>
                      {tech.role}
                    </p>
                  </a>
                ))}
              </div>
            </div>

            <div className="chrome-panel p-5 md:p-6 visit-about-panel">
              <h2 className="subtle-label mb-4">About the Name</h2>
              <p className="text-sm leading-relaxed mb-3" style={{ color: 'var(--text-2)' }}>
                Named after the <em>San José</em>, a sunken ship whose rediscovery revealed immense,
                carefully preserved value beneath the surface. PRIDE is a vast repository where the
                real value is hidden in raw experimental data, waiting to be systematically recovered,
                catalogued, and understood.
              </p>
              <p className="text-xs mb-4" style={{ color: 'var(--text-3)' }}>
                Johannes Gutenberg University Mainz | HPC: Mogon2 SLURM Cluster
              </p>
              <div className="flex gap-2 flex-wrap">
                <span className="metric-pill">Open science infrastructure</span>
                <span className="metric-pill">Containerized reproducibility</span>
                <span className="metric-pill">Model-ready exports</span>
              </div>
            </div>
          </div>
        </section>

        <footer className="px-6 py-10 text-center landing-footer">
          {hasExploreAction && (
            <button className="btn-primary landing-cta-primary" onClick={handleExplore}>
              Explore the Data
            </button>
          )}
          <div className="mt-3">
            <button className="btn-secondary landing-cta-secondary" onClick={onBack}>
              Back to Home
            </button>
          </div>
        </footer>
      </main>
    </div>
  );
}

function KpiCard({ label, value }: { label: string; value: number | undefined }) {
  return (
    <div className="visit-kpi-card">
      <span className="visit-kpi-label">{label}</span>
      <span className="visit-kpi-value mono">
        {value != null ? formatNumber(value) : '\u2014'}
      </span>
    </div>
  );
}
