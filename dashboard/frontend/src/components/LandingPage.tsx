import { useQuery } from '@tanstack/react-query';
import { getCollectionInfo, getStats } from '../api';
import SiteNav from './SiteNav';
import PipelineOverview from './PipelineOverview';

interface LandingPageProps {
  onViewSummary: () => void;
  explorerUrl?: string;
  onExploreData?: () => void;
}

const valueProps = [
  {
    title: 'Triple Orthogonal Validation',
    description:
      'Every dataset processed with FragPipe, DIA-NN, and Sage independently. We store both consensus and disagreement as first-class scientific data.',
  },
  {
    title: 'Bias-Aware by Design',
    description:
      'Lab identity, organism, gradient length, column type, and acquisition mode are tracked explicitly \u2014 enabling stratified sampling and cross-lab validation.',
  },
  {
    title: 'Full 4D Raw Signal',
    description:
      'Beyond identifications: raw retention time, m/z, ion mobility, and intensity traces extracted directly from timsTOF .d files via rustims.',
  },
  {
    title: 'Versioned Snapshots',
    description:
      'Frozen, reproducible, citable datasets. Always rebuildable from PRIDE + pipeline version. Multiple models trainable from the same snapshot.',
  },
  {
    title: 'Human Quality Gates',
    description:
      'Automation runs the pipeline; humans guard scientific integrity. Five checkpoint types ensure dataset selection, QC, consensus rules, and releases meet standards.',
  },
  {
    title: 'Open Science',
    description:
      'Built on PRIDE, the largest public proteomics repository. Leverages open-source tools (Sage, rustims, Snakemake) and standard formats (Parquet).',
  },
];

const whyNowCards = [
  {
    title: 'Sage',
    subtitle: 'Multi-engine validation at scale',
    description:
      'A fast, high-quality open-source search engine makes triple-engine processing economically feasible for thousands of datasets.',
  },
  {
    title: 'rustims / imspy',
    subtitle: '4D feature extraction now feasible',
    description:
      'Rust-powered raw data access enables extraction of full chromatograms, mobilograms, isotope envelopes, and peak shapes from timsTOF data.',
  },
  {
    title: 'LLM-Assisted Engineering',
    subtitle: 'Reduced infrastructure cost',
    description:
      'Complex bioinformatics pipelines \u2014 orchestration, containerization, QC automation \u2014 can now be built and iterated at unprecedented speed.',
  },
];

function formatNumber(n: number): string {
  if (n >= 1_000_000) return `${(n / 1_000_000).toFixed(1)}M`;
  if (n >= 1_000) return `${(n / 1_000).toFixed(0)}K`;
  return n.toLocaleString();
}

export default function LandingPage({
  onViewSummary,
  explorerUrl,
  onExploreData,
}: LandingPageProps) {
  // Try fetching live stats (graceful if backend unavailable)
  const { data: collectionInfo } = useQuery({
    queryKey: ['collection-info-landing'],
    queryFn: getCollectionInfo,
    retry: false,
    staleTime: 60_000,
  });

  const { data: stats } = useQuery({
    queryKey: ['stats-landing'],
    queryFn: getStats,
    retry: false,
    enabled: !collectionInfo,
    staleTime: 60_000,
  });

  const totalPrecursors = collectionInfo?.n_total_precursors ?? stats?.total_precursors;
  const totalDatasets = collectionInfo?.n_datasets;
  const totalStudies = collectionInfo?.n_studies;

  return (
    <div className="app-shell min-h-screen flex flex-col">
      <SiteNav
        currentPage="landing"
        onNavigateLanding={() => {}}
        onNavigateVisit={onViewSummary}
        explorerUrl={explorerUrl}
        onExploreData={onExploreData}
      />

      <main className="flex-1 overflow-y-auto">
        {/* Hero */}
        <section className="hero-section py-20 md:py-28 px-6 text-center">
          <div className="max-w-3xl mx-auto reveal-up">
            <h1
              className="text-4xl md:text-5xl font-bold tracking-tight mb-4"
              style={{ color: 'var(--accent-1)' }}
            >
              San Jose
            </h1>
            <p className="text-lg md:text-xl mb-3" style={{ color: 'var(--text-2)' }}>
              A Reproducible, Bias-Aware Reference Layer
              <br />
              for timsTOF Data on PRIDE
            </p>
            <p
              className="text-sm md:text-base italic max-w-xl mx-auto mb-8"
              style={{ color: 'var(--text-3)' }}
            >
              &ldquo;We are not collecting peptides. We are collecting peptide observations
              in experimental context.&rdquo;
            </p>

            {/* Live stats badge */}
            {totalPrecursors != null && (
              <div className="flex justify-center gap-3 mb-8 flex-wrap">
                <span className="metric-pill">
                  <span className="metric-value">{formatNumber(totalPrecursors)}</span> precursors
                </span>
                {totalDatasets != null && (
                  <span className="metric-pill">
                    <span className="metric-value">{totalDatasets}</span> datasets
                  </span>
                )}
                {totalStudies != null && (
                  <span className="metric-pill">
                    <span className="metric-value">{totalStudies}</span> studies
                  </span>
                )}
              </div>
            )}

            {/* CTAs */}
            <div className="flex justify-center gap-3 flex-wrap">
              {(explorerUrl || onExploreData) && (
                <button
                  className="btn-primary"
                  style={{ padding: '0.6rem 1.5rem', fontSize: '0.9rem' }}
                  onClick={() => {
                    if (explorerUrl) window.open(explorerUrl, '_blank', 'noopener');
                    else onExploreData?.();
                  }}
                >
                  Explore Data
                </button>
              )}
              <button
                className="btn-secondary"
                style={{ padding: '0.6rem 1.5rem', fontSize: '0.9rem' }}
                onClick={onViewSummary}
              >
                Project Summary
              </button>
            </div>
          </div>
        </section>

        {/* Value Propositions */}
        <section className="px-6 py-12 max-w-6xl mx-auto">
          <h2 className="subtle-label mb-6 text-center">What Makes San Jose Unique</h2>
          <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-3 gap-4">
            {valueProps.map((vp) => (
              <div key={vp.title} className="chrome-panel p-5">
                <h3
                  className="text-sm font-bold mb-2"
                  style={{ color: 'var(--accent-2)' }}
                >
                  {vp.title}
                </h3>
                <p className="text-sm leading-relaxed" style={{ color: 'var(--text-2)' }}>
                  {vp.description}
                </p>
              </div>
            ))}
          </div>
        </section>

        {/* Pipeline Overview */}
        <section className="px-6 py-12 max-w-5xl mx-auto">
          <h2 className="subtle-label mb-2 text-center">Processing Pipeline</h2>
          <p className="text-sm text-center mb-6" style={{ color: 'var(--text-3)' }}>
            Each PRIDE dataset passes through a 6-step runner pipeline with checkpointing and QC validation.
          </p>
          <PipelineOverview />
        </section>

        {/* Why Now */}
        <section className="px-6 py-12 max-w-5xl mx-auto">
          <h2 className="subtle-label mb-2 text-center">Why Now?</h2>
          <p className="text-sm text-center mb-6" style={{ color: 'var(--text-3)' }}>
            Three enabling technologies have converged to make systematic reprocessing at scale feasible.
          </p>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            {whyNowCards.map((card) => (
              <div key={card.title} className="chrome-panel p-5">
                <h3
                  className="mono text-sm font-bold mb-1"
                  style={{ color: 'var(--accent-1)' }}
                >
                  {card.title}
                </h3>
                <p
                  className="text-xs font-semibold mb-2"
                  style={{ color: 'var(--text-2)' }}
                >
                  {card.subtitle}
                </p>
                <p className="text-sm leading-relaxed" style={{ color: 'var(--text-3)' }}>
                  {card.description}
                </p>
              </div>
            ))}
          </div>
        </section>

        {/* Footer */}
        <footer className="px-6 py-10 text-center" style={{ borderTop: '1px solid var(--border-2)' }}>
          <p className="text-sm mb-1" style={{ color: 'var(--text-2)' }}>
            Johannes Gutenberg University Mainz
          </p>
          <p className="text-xs" style={{ color: 'var(--text-3)' }}>
            Built on PRIDE &middot; Powered by Sage, rustims, and Snakemake
          </p>
        </footer>
      </main>
    </div>
  );
}
