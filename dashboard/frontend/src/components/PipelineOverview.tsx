import { useState, useRef, useEffect } from 'react';

const steps = [
  {
    number: 1,
    name: 'Download',
    description: 'PRIDE raw .d files',
    icon: '\u2193',
    back: {
      input: 'PRIDE accession',
      output: 'raw/{acc}/*.d',
      detail: 'Streams Bruker .d folders from PRIDE. Validates checksums, skips non-timsTOF files.',
      cmd: 'run_dataset.py PXD019086',
    },
  },
  {
    number: 2,
    name: 'Search',
    description: 'FragPipe + DIA-NN + Sage',
    icon: '\u2759',
    back: {
      input: '.d files + FASTA',
      output: 'PSMs per engine',
      detail: 'Triple orthogonal search \u2014 spectrum- and peptide-centric. All results unfiltered.',
      cmd: '--no-fdr-filter',
    },
  },
  {
    number: 3,
    name: 'Stratify',
    description: 'Engine consensus',
    icon: '\u25B3',
    back: {
      input: '3 engine results',
      output: 'precursor_index',
      detail: 'Consensus tiers: 3/3, 2/3, or 1/3 engine agreement. Unified precursor index.',
      cmd: 'build_precursor_index',
    },
  },
  {
    number: 4,
    name: 'Extract',
    description: '4D raw signal',
    icon: '\u25CE',
    back: {
      input: '.d + precursor index',
      output: 'blobs.bin',
      detail: 'Full RT \u00d7 m/z \u00d7 IM \u00d7 intensity per precursor via rustims. Compressed blobs.',
      cmd: 'precursor_store_parquet',
    },
  },
  {
    number: 5,
    name: 'Merge',
    description: 'IDs + features joined',
    icon: '\u29D6',
    back: {
      input: 'Index + 4D blobs',
      output: 'precursor_store',
      detail: 'Joins IDs with raw signal features into one Parquet store. Adds QC metrics.',
      cmd: 'precursor_store --merge',
    },
  },
  {
    number: 6,
    name: 'Package',
    description: 'Versioned archive',
    icon: '\u25A0',
    back: {
      input: 'All outputs',
      output: '{acc}_v{ver}.zip',
      detail: 'Self-contained archive: manifest, engines, consensus, raw features.',
      cmd: '--package --version 1.0',
    },
  },
];

interface PipelineOverviewProps {
  compact?: boolean;
}

function FlipCard({
  step,
  compact,
  index = 0,
}: {
  step: (typeof steps)[number];
  compact?: boolean;
  index?: number;
}) {
  const [face, setFace] = useState<'front' | 'back'>('front');
  const [animating, setAnimating] = useState(false);
  const [entranceFlip, setEntranceFlip] = useState(false);
  const nextFace = useRef<'front' | 'back'>('back');

  // Staggered entrance flip on mount
  useEffect(() => {
    const delay = 300 + index * 150; // stagger each card by 150ms
    const timer = setTimeout(() => setEntranceFlip(true), delay);
    return () => clearTimeout(timer);
  }, [index]);

  function handleClick() {
    if (animating) return;
    nextFace.current = face === 'front' ? 'back' : 'front';
    setAnimating(true);
  }

  function handleHalf() {
    // Mid-point of the animation — swap content
    setFace(nextFace.current);
  }

  function handleDone() {
    setAnimating(false);
  }

  const isFront = face === 'front';
  const classes = [
    'pipeline-flip-card',
    animating ? 'pipeline-flip-card--animating' : '',
    entranceFlip ? 'pipeline-flip-card--entered' : 'pipeline-flip-card--hidden',
  ].filter(Boolean).join(' ');

  return (
    <div
      className={classes}
      onClick={handleClick}
      onAnimationIteration={handleHalf}
      onAnimationEnd={handleDone}
    >
      {isFront ? (
        <div
          className={`pipeline-step ${compact ? 'p-3' : 'p-4'} h-full flex flex-col items-center justify-center text-center`}
        >
          <span
            className="mono text-xs font-bold mb-1"
            style={{ color: 'var(--accent-1)' }}
          >
            {step.icon}
          </span>
          <span
            className="mono text-xs font-bold mb-0.5"
            style={{ color: 'var(--text-3)' }}
          >
            Step {step.number}
          </span>
          <span
            className="text-sm font-semibold"
            style={{ color: 'var(--text-1)' }}
          >
            {step.name}
          </span>
          <span className="text-xs mt-0.5" style={{ color: 'var(--text-3)' }}>
            {step.description}
          </span>
          <span className="pipeline-flip-hint">click to flip</span>
        </div>
      ) : (
        <div className="pipeline-flip-back-card h-full p-3">
          <div className="pipeline-flip-back-content">
            <span
              className="mono text-xs font-bold"
              style={{ color: 'var(--accent-1)' }}
            >
              {step.name}
            </span>

            <div className="pipeline-flip-row">
              <span className="pipeline-flip-label">IN</span>
              <span className="pipeline-flip-value">{step.back.input}</span>
            </div>

            <div className="pipeline-flip-row">
              <span className="pipeline-flip-label">OUT</span>
              <span className="pipeline-flip-value mono" style={{ fontSize: '0.68rem' }}>
                {step.back.output}
              </span>
            </div>

            <p className="pipeline-flip-detail">{step.back.detail}</p>

            <code className="pipeline-flip-cmd">{step.back.cmd}</code>
          </div>
        </div>
      )}
    </div>
  );
}

export default function PipelineOverview({ compact }: PipelineOverviewProps) {
  return (
    <div className={`flex flex-col ${compact ? 'gap-3' : 'gap-4'}`}>
      {/* Desktop: horizontal */}
      <div className="hidden md:flex items-stretch gap-0 pipeline-flow">
        {steps.map((step, i) => (
          <div key={step.number} className="flex items-stretch flex-1 min-w-0">
            <FlipCard step={step} compact={compact} index={i} />
            {i < steps.length - 1 && <div className="pipeline-connector" />}
          </div>
        ))}
      </div>

      {/* Mobile: vertical */}
      <div className="flex md:hidden flex-col gap-2">
        {steps.map((step, i) => (
          <div key={step.number} className="flex flex-col items-stretch">
            <FlipCard step={step} compact={compact} index={i} />
            {i < steps.length - 1 && (
              <div
                className="w-px h-3 mx-auto"
                style={{ background: 'var(--border-1)' }}
              />
            )}
          </div>
        ))}
      </div>
    </div>
  );
}
