import { useState } from 'react';
import { useQuery } from '@tanstack/react-query';
import { getCollectionInfo, getStats } from '../api';
import { useScrollReveal } from '../hooks/useScrollReveal';
import SiteNav from './SiteNav';
import PipelineOverview from './PipelineOverview';
import IonCloudCanvas from './IonCloudCanvas';

interface LandingPageProps {
  onViewSummary: () => void;
  onNavigateBlueprint?: () => void;
  explorerUrl?: string;
  onExploreData?: () => void;
}

const valueProps = [
  {
    title: 'Triple Orthogonal Validation',
    description:
      'Every dataset processed with FragPipe, DIA-NN, and Sage independently. We store both consensus and disagreement as first-class scientific data.',
    details: [
      'FragPipe provides spectrum-centric database search with MSFragger at its core.',
      'DIA-NN runs in library-free DDA mode, offering a peptide-centric perspective on the same data.',
      'Sage adds a fast, fully open-source Rust engine as a third independent opinion.',
      'Precursors are stratified by engine agreement: 3/3, 2/3, or single-engine hits \u2014 all retained for downstream analysis.',
    ],
  },
  {
    title: 'Bias-Aware by Design',
    description:
      'Lab identity, organism, gradient length, column type, and acquisition mode are tracked explicitly \u2014 enabling stratified sampling and cross-lab validation.',
    details: [
      'Every observation is tagged with its experimental context: lab, instrument, column type, gradient length, and acquisition mode.',
      'This enables stratified train/test splits that prevent data leakage across labs or instruments.',
      'Cross-lab validation becomes possible by design, not as an afterthought.',
    ],
  },
  {
    title: 'Full 4D Raw Signal',
    description:
      'Beyond identifications: raw retention time, m/z, ion mobility, and intensity traces extracted directly from timsTOF .d files via rustims.',
    details: [
      'XICs (extracted ion chromatograms), mobilograms, and isotope envelopes stored as variable-length arrays in Parquet.',
      'Fragment spectra include matched b/y ions with ppm-level error annotation.',
      'All raw features extracted via rustims \u2014 a Rust backend that reads Bruker .d files natively without mzML conversion.',
    ],
  },
  {
    title: 'Versioned Snapshots',
    description:
      'Frozen, reproducible, citable datasets. Always rebuildable from PRIDE + pipeline version. Multiple models trainable from the same snapshot.',
    details: [
      'Each snapshot is a frozen, immutable slice of the database at a point in time.',
      'Snapshots record the exact pipeline version, config, and PRIDE accessions used \u2014 fully reproducible.',
      'Multiple prediction models (CCS, RT, MS2) can be trained and compared on the same snapshot.',
    ],
  },
  {
    title: 'Human Quality Gates',
    description:
      'Automation runs the pipeline; humans guard scientific integrity. Five checkpoint types ensure dataset selection, QC, consensus rules, and releases meet standards.',
    details: [
      'Dataset selection: curated PRIDE accessions with metadata normalization.',
      'Automated QC: transparent pass/fail thresholds with detailed logs.',
      'Consensus review: engine agreement rules reviewed before stratification.',
      'Release approval: manual sign-off required before any versioned snapshot is published.',
    ],
  },
  {
    title: 'Open Science',
    description:
      'Built on PRIDE, the largest public proteomics repository. Leverages open-source tools (Sage, rustims, Snakemake) and standard formats (Parquet).',
    details: [
      'All source data comes from PRIDE, ensuring public availability and traceability.',
      'Core tools are open-source: Sage (search), rustims (raw access), Snakemake (orchestration).',
      'Data stored in Apache Parquet for efficient columnar access and broad ecosystem compatibility.',
    ],
  },
];

const qualityGates = [
  'Dataset selection and metadata normalization',
  'Automated QC thresholds with transparent pass/fail logs',
  'Consensus rule review across all three search engines',
  'Manual release approval for each versioned snapshot',
  'Reproducibility checks against frozen pipeline versions',
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
  onNavigateBlueprint,
  explorerUrl,
  onExploreData,
}: LandingPageProps) {
  const [expandedCard, setExpandedCard] = useState<number | null>(null);

  const uniqueSection = useScrollReveal();
  const pipelineSection = useScrollReveal();
  const whyNowSection = useScrollReveal();

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
  const hasExploreAction = Boolean(explorerUrl || onExploreData);

  const heroStats = [
    { label: 'Precursors', value: totalPrecursors != null ? formatNumber(totalPrecursors) : '--' },
    { label: 'Datasets', value: totalDatasets != null ? String(totalDatasets) : '--' },
    { label: 'Studies', value: totalStudies != null ? String(totalStudies) : '--' },
  ];

  const handleExplore = () => {
    if (explorerUrl) window.open(explorerUrl, '_blank', 'noopener');
    else onExploreData?.();
  };

  return (
    <div className="app-shell min-h-screen flex flex-col landing-page">
      <SiteNav
        currentPage="landing"
        onNavigateLanding={() => {}}
        onNavigateVisit={onViewSummary}
        onNavigateBlueprint={onNavigateBlueprint}
        explorerUrl={explorerUrl}
        onExploreData={onExploreData}
      />

      <main className="flex-1 overflow-y-auto">
        <section className="landing-hero px-6 pt-14 pb-12 md:pt-20 md:pb-16">
          <IonCloudCanvas className="ion-cloud-canvas" />
          <div className="max-w-6xl mx-auto landing-hero-grid reveal-up">
            <div className="landing-hero-content">
              <p className="subtle-label mb-3">Reference Layer For PRIDE timsTOF</p>
              <h1 className="landing-hero-title mb-4">
                <span className="landing-hero-title-main">San José</span>
              </h1>
              <p className="text-lg md:text-xl mb-4" style={{ color: 'var(--text-2)' }}>
                A reproducible, bias-aware reference layer for timsTOF data on PRIDE.
              </p>
              <p className="text-sm md:text-base italic mb-7 landing-quote">
                "We are not collecting peptides. We are collecting peptide observations in experimental context."
              </p>

              <div className="landing-hero-actions mt-3">
                {hasExploreAction && (
                  <button
                    className="btn-primary landing-cta-primary"
                    onClick={handleExplore}
                  >
                    Explore Data
                  </button>
                )}
                <button
                  className="btn-secondary landing-cta-secondary"
                  onClick={onViewSummary}
                >
                  Project Summary
                </button>
              </div>

              <div className="flex gap-2 mt-5 flex-wrap">
                <span className="metric-pill metric-pill-shimmer">Triple-engine consensus</span>
                <span className="metric-pill metric-pill-shimmer">4D raw signal extraction</span>
                <span className="metric-pill metric-pill-shimmer">Versioned and citable snapshots</span>
              </div>
            </div>

            <aside className="chrome-panel landing-snapshot p-5 md:p-6">
              <p className="subtle-label mb-3">Live Collection Snapshot</p>
              <div className="landing-stat-grid mb-4">
                {heroStats.map((stat) => (
                  <div key={stat.label} className="landing-stat-card">
                    <span className="landing-stat-label">{stat.label}</span>
                    <span className="landing-stat-value mono">{stat.value}</span>
                  </div>
                ))}
              </div>

              <div className="landing-divider" />
              <h2 className="text-sm font-semibold mb-2" style={{ color: 'var(--text-1)' }}>
                Built for reliable model training
              </h2>
              <ul className="landing-checklist">
                <li>Bias-aware metadata tracked per acquisition context</li>
                <li>Consensus and disagreement both retained as signal</li>
                <li>Frozen snapshots reproducible from source data</li>
              </ul>
            </aside>
          </div>
        </section>

        <section
          ref={uniqueSection.ref}
          className={`px-6 py-10 max-w-6xl mx-auto scroll-section${uniqueSection.visible ? ' scroll-section--visible' : ''}`}
        >
          <div className="flex flex-col items-center text-center mb-6">
            <h2 className="subtle-label mb-2">What Makes San José Unique</h2>
            <p className="text-sm max-w-2xl" style={{ color: 'var(--text-3)' }}>
              Scientific rigor at scale requires both strong automation and traceable, explainable quality control.
            </p>
          </div>
          <div className={`grid grid-cols-1 md:grid-cols-2 xl:grid-cols-3 gap-4 stagger-grid${uniqueSection.visible ? ' stagger-grid--visible' : ''}`}>
            {valueProps.map((vp, index) => {
              const isOpen = expandedCard === index;
              return (
                <article
                  key={vp.title}
                  className={`chrome-panel p-5 landing-value-card landing-value-card--clickable${isOpen ? ' landing-value-card--open' : ''}`}
                  onClick={() => setExpandedCard(isOpen ? null : index)}
                >
                  <div className="flex items-start justify-between gap-2">
                    <span className="landing-value-index mono">{String(index + 1).padStart(2, '0')}</span>
                    <span className={`landing-value-chevron${isOpen ? ' landing-value-chevron--open' : ''}`}>&#x25BE;</span>
                  </div>
                  <h3 className="text-sm font-bold mb-2" style={{ color: 'var(--accent-2)' }}>
                    {vp.title}
                  </h3>
                  <p className="text-sm leading-relaxed" style={{ color: 'var(--text-2)' }}>
                    {vp.description}
                  </p>
                  <div className={`landing-value-details${isOpen ? ' landing-value-details--open' : ''}`}>
                    <ul className="landing-checklist mt-3">
                      {vp.details.map((d) => (
                        <li key={d}>{d}</li>
                      ))}
                    </ul>
                  </div>
                </article>
              );
            })}
          </div>
        </section>

        <section
          ref={pipelineSection.ref}
          className={`px-6 py-10 max-w-6xl mx-auto scroll-section${pipelineSection.visible ? ' scroll-section--visible' : ''}`}
        >
          <div className="grid grid-cols-1 lg:grid-cols-[1.45fr,1fr] gap-4">
            <div className="chrome-panel p-5 md:p-6">
              <h2 className="subtle-label mb-2 text-center">Processing Pipeline</h2>
              <p className="text-sm text-center mb-6" style={{ color: 'var(--text-3)' }}>
                Each PRIDE dataset passes through a 6-step runner pipeline with checkpointing and QC validation.
              </p>
              <PipelineOverview />
            </div>

            <div className="chrome-panel p-5 md:p-6 landing-quality-panel">
              <h2 className="subtle-label mb-2">Human Quality Gates</h2>
              <p className="text-sm mb-4" style={{ color: 'var(--text-3)' }}>
                Automation runs continuously, while release decisions stay auditable and human-reviewed.
              </p>
              <ul className="landing-checklist">
                {qualityGates.map((gate) => (
                  <li key={gate}>{gate}</li>
                ))}
              </ul>
            </div>
          </div>
        </section>

        <section
          ref={whyNowSection.ref}
          className={`px-6 py-10 max-w-6xl mx-auto scroll-section${whyNowSection.visible ? ' scroll-section--visible' : ''}`}
        >
          <div className="flex flex-col items-center text-center mb-6">
            <h2 className="subtle-label mb-2">Why Now?</h2>
            <p className="text-sm max-w-2xl" style={{ color: 'var(--text-3)' }}>
              Three enabling technologies have converged to make systematic reprocessing at scale feasible.
            </p>
          </div>
          <div className={`grid grid-cols-1 md:grid-cols-3 gap-4 stagger-grid${whyNowSection.visible ? ' stagger-grid--visible' : ''}`}>
            {whyNowCards.map((card) => (
              <article key={card.title} className="chrome-panel p-5 landing-why-card">
                <p className="subtle-label mb-1">{card.title}</p>
                <h3 className="text-sm font-semibold mb-2" style={{ color: 'var(--accent-1)' }}>
                  {card.subtitle}
                </h3>
                <p className="text-sm leading-relaxed" style={{ color: 'var(--text-2)' }}>
                  {card.description}
                </p>
              </article>
            ))}
          </div>
        </section>

        <footer className="px-6 py-10 text-center landing-footer">
          <p className="text-sm mb-1" style={{ color: 'var(--text-2)' }}>
            Johannes Gutenberg University Mainz
          </p>
          <p className="text-xs mb-5" style={{ color: 'var(--text-3)' }}>
            Built on PRIDE | Powered by Sage, rustims, and Snakemake
          </p>
          {hasExploreAction && (
            <button className="btn-primary landing-cta-primary" onClick={handleExplore}>
              Explore the Data
            </button>
          )}
          <div className="mt-3">
            <button className="btn-secondary landing-cta-secondary" onClick={onViewSummary}>
              View Project Summary
            </button>
          </div>
        </footer>
      </main>
    </div>
  );
}
