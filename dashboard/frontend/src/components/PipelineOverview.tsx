const steps = [
  { number: 1, name: 'Download', description: 'PRIDE raw .d files', icon: '\u2193' },
  { number: 2, name: 'Search', description: 'FragPipe + DIA-NN + Sage', icon: '\u2759' },
  { number: 3, name: 'Stratify', description: 'Engine consensus', icon: '\u25B3' },
  { number: 4, name: 'Extract', description: '4D raw signal', icon: '\u25CE' },
  { number: 5, name: 'Merge', description: 'IDs + features joined', icon: '\u29D6' },
  { number: 6, name: 'Package', description: 'Versioned archive', icon: '\u25A0' },
];

interface PipelineOverviewProps {
  compact?: boolean;
}

export default function PipelineOverview({ compact }: PipelineOverviewProps) {
  return (
    <div className={`flex flex-col ${compact ? 'gap-3' : 'gap-4'}`}>
      {/* Desktop: horizontal */}
      <div className="hidden md:flex items-stretch gap-0 pipeline-flow">
        {steps.map((step, i) => (
          <div key={step.number} className="flex items-stretch flex-1 min-w-0">
            <div className={`pipeline-step flex-1 ${compact ? 'p-3' : 'p-4'} flex flex-col items-center text-center`}>
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
              <span className="text-sm font-semibold" style={{ color: 'var(--text-1)' }}>
                {step.name}
              </span>
              <span className="text-xs mt-0.5" style={{ color: 'var(--text-3)' }}>
                {step.description}
              </span>
            </div>
            {i < steps.length - 1 && <div className="pipeline-connector" />}
          </div>
        ))}
      </div>

      {/* Mobile: vertical */}
      <div className="flex md:hidden flex-col gap-2">
        {steps.map((step, i) => (
          <div key={step.number} className="flex flex-col items-stretch">
            <div className="pipeline-step p-3 flex items-center gap-3">
              <span
                className="mono text-sm font-bold shrink-0"
                style={{ color: 'var(--accent-1)', width: '2rem', textAlign: 'center' }}
              >
                {step.number}
              </span>
              <div className="min-w-0">
                <span className="text-sm font-semibold" style={{ color: 'var(--text-1)' }}>
                  {step.name}
                </span>
                <span className="text-xs ml-2" style={{ color: 'var(--text-3)' }}>
                  {step.description}
                </span>
              </div>
            </div>
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
