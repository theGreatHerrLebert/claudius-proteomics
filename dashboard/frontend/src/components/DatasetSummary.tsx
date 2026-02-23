import { useQuery } from '@tanstack/react-query';
import { getDatasetFullSummary } from '../api';
import type { DatasetFullSummary, QualitySummaryData, EngineOverlapStats, QualityDistBucket } from '../api';

interface DatasetSummaryPageProps {
  activeDataset: string | null;
  onBrowse: () => void;
}

function fmtNum(n: number): string {
  return n.toLocaleString();
}

function fmtPct(n: number): string {
  return `${n.toFixed(1)}%`;
}

function fmtRange(range: [number, number] | null, digits = 1, suffix = ''): string {
  if (!range) return 'N/A';
  return `${range[0].toFixed(digits)}${suffix} – ${range[1].toFixed(digits)}${suffix}`;
}

function qualityColor(value: number | null): string {
  if (value === null) return 'text-slate-400';
  if (value >= 0.9) return 'status-good';
  if (value >= 0.7) return 'status-warn';
  return 'status-bad';
}

// ─── Engine Agreement Stacked Bar + Per-Engine IDs ──────────────────────

function EngineAgreementCard({
  quality,
  engines,
  unidentified,
  overlap,
  uniquePeptides,
}: {
  quality: QualitySummaryData | null;
  engines: Record<string, number> | null;
  unidentified: number;
  overlap: EngineOverlapStats | null;
  uniquePeptides: Record<string, number> | null;
}) {
  const segments = quality ? [
    { count: quality.n_0_engines, pct: quality.pct_0_engines, label: '0', color: '#334155' },
    { count: quality.n_1_engines, pct: quality.pct_1_engines, label: '1', color: '#d97706' },
    { count: quality.n_2_engines, pct: quality.pct_2_engines, label: '2', color: '#0ea5e9' },
    { count: quality.n_3_engines, pct: quality.pct_3_engines, label: '3', color: '#10b981' },
  ] : [];

  const total = segments.reduce((s, seg) => s + seg.count, 0);

  const engineEntries = [
    { key: 'fragpipe', label: 'FragPipe', color: '#06b6d4' },
    { key: 'diann', label: 'DIA-NN', color: '#10b981' },
    { key: 'sage', label: 'Sage', color: '#0ea5e9' },
  ];

  return (
    <div className="space-y-3">
      {/* Stacked bar */}
      {total > 0 && (
        <div>
          <div className="h-5 rounded-md overflow-hidden flex bg-slate-900/70 border border-[#2a4368]">
            {segments.map((seg, i) => {
              const w = Math.max((seg.count / total) * 100, seg.count > 0 ? 0.5 : 0);
              return (
                <div
                  key={i}
                  className="h-full transition-all"
                  style={{ width: `${w}%`, background: seg.color, opacity: 0.85 }}
                  title={`${seg.label} engine${seg.label !== '1' ? 's' : ''}: ${fmtNum(seg.count)} (${fmtPct(seg.pct)})`}
                />
              );
            })}
          </div>
          <div className="flex flex-wrap gap-1.5 mt-1.5">
            {segments.map((seg, i) => (
              <span key={i} className="metric-pill text-[10px]">
                <span className="inline-block w-1.5 h-1.5 rounded-full" style={{ background: seg.color }} />
                {seg.label}
                <span className="metric-value">{fmtNum(seg.count)}</span>
                <span className="text-slate-500">{fmtPct(seg.pct)}</span>
              </span>
            ))}
          </div>
        </div>
      )}

      {/* Per-engine PSMs + unique peptides */}
      {engines && (
        <div className="border-t border-[#25415f] pt-2">
          <div className="grid grid-cols-3 gap-2">
            {engineEntries.map((e) => {
              const psms = engines[e.key] ?? 0;
              const peps = uniquePeptides?.[e.key];
              return (
                <div key={e.key} className="panel-inset p-2">
                  <div className="text-[10px] font-semibold uppercase tracking-wide" style={{ color: e.color }}>{e.label}</div>
                  <div className="mono text-sm text-[#8dd6e8]">{fmtNum(psms)}</div>
                  <div className="text-[10px] text-slate-400">PSMs</div>
                  {peps !== undefined && (
                    <>
                      <div className="mono text-xs text-[#afc0db] mt-0.5">{fmtNum(peps)}</div>
                      <div className="text-[10px] text-slate-500">unique peptides</div>
                    </>
                  )}
                </div>
              );
            })}
          </div>
          {unidentified > 0 && (
            <div className="mt-1.5">
              <span className="metric-pill text-[10px]">Unidentified <span className="metric-value">{fmtNum(unidentified)}</span></span>
            </div>
          )}
        </div>
      )}

      {/* Overlap details */}
      {overlap && (
        <div className="border-t border-[#25415f] pt-2 space-y-1.5">
          <div className="flex flex-wrap gap-1.5">
            <span className="metric-pill text-[10px]">3-way <span className="metric-value">{fmtNum(overlap.n_all_three)}</span> <span className="text-slate-500">{fmtPct(overlap.three_way_rate * 100)}</span></span>
            <span className="metric-pill text-[10px]">2+ <span className="metric-value">{fmtNum(overlap.n_at_least_two)}</span> <span className="text-slate-500">{fmtPct(overlap.at_least_two_rate * 100)}</span></span>
            <span className="metric-pill text-[10px]">union <span className="metric-value">{fmtNum(overlap.n_union)}</span></span>
          </div>
          <div className="flex flex-wrap gap-1 text-[10px]">
            <span className="metric-pill">FP+DN <span className="metric-value">{fmtNum(overlap.n_fp_dn_only)}</span></span>
            <span className="metric-pill">FP+SG <span className="metric-value">{fmtNum(overlap.n_fp_sg_only)}</span></span>
            <span className="metric-pill">DN+SG <span className="metric-value">{fmtNum(overlap.n_dn_sg_only)}</span></span>
            <span className="metric-pill">FP only <span className="metric-value">{fmtNum(overlap.n_fragpipe_only)}</span></span>
            <span className="metric-pill">DN only <span className="metric-value">{fmtNum(overlap.n_diann_only)}</span></span>
            <span className="metric-pill">SG only <span className="metric-value">{fmtNum(overlap.n_sage_only)}</span></span>
          </div>
        </div>
      )}
    </div>
  );
}

// ─── Quality Metrics + Distributions ────────────────────────────────────

function QualityHistogram({ buckets, label }: { buckets: QualityDistBucket[]; label: string }) {
  const maxPct = Math.max(...buckets.map(b => b.pct), 1);
  const barH = 32;

  return (
    <div>
      <div className="text-[10px] text-slate-400 mb-1">{label}</div>
      <div className="flex items-end gap-px h-9">
        {buckets.map((b, i) => {
          const h = Math.max((b.pct / maxPct) * barH, b.count > 0 ? 2 : 0);
          // Color: red for low R², amber for mid, green for high
          const color = b.upper <= 0.5 ? '#ef4444' : b.upper <= 0.7 ? '#f59e0b' : b.upper <= 0.9 ? '#0ea5e9' : '#10b981';
          return (
            <div
              key={i}
              className="flex-1 rounded-t-sm"
              style={{ height: `${h}px`, background: color, opacity: 0.75 }}
              title={`${b.lower.toFixed(2)}–${b.upper.toFixed(2)}: ${fmtNum(b.count)} (${fmtPct(b.pct)})`}
            />
          );
        })}
      </div>
      <div className="flex justify-between text-[9px] text-slate-500 mt-0.5">
        <span>0</span>
        <span>0.5</span>
        <span>1.0</span>
      </div>
    </div>
  );
}

function QualityMetricsCard({
  quality,
  distributions,
}: {
  quality: QualitySummaryData | null;
  distributions: Record<string, QualityDistBucket[]> | null;
}) {
  if (!quality) return <div className="text-xs text-slate-400">No quality metrics available</div>;

  const metrics = [
    { key: 'ms1_rt_r2', label: 'RT R²', mean: quality.ms1_rt_r2_mean, median: quality.ms1_rt_r2_median },
    { key: 'ms1_im_r2', label: 'IM R²', mean: quality.ms1_im_r2_mean, median: quality.ms1_im_r2_median },
    { key: 'isotope_cosim', label: 'Isotope Cosim', mean: quality.isotope_cosim_mean, median: quality.isotope_cosim_median },
  ];

  return (
    <div className="space-y-3">
      {/* Metric mean/median + histogram side by side */}
      {metrics.map((m) => (
        <div key={m.key} className="flex items-start gap-3">
          <div className="flex-none w-28">
            <div className="text-[11px] text-slate-300 font-semibold">{m.label}</div>
            <div className="flex gap-1.5 mt-0.5">
              <span className="metric-pill text-[10px]">
                med <span className={`metric-value ${qualityColor(m.median)}`}>{m.median?.toFixed(3) ?? 'N/A'}</span>
              </span>
            </div>
            <div className="flex gap-1.5 mt-0.5">
              <span className="metric-pill text-[10px]">
                avg <span className={`metric-value ${qualityColor(m.mean)}`}>{m.mean?.toFixed(3) ?? 'N/A'}</span>
              </span>
            </div>
          </div>
          {distributions?.[m.key] && (
            <div className="flex-1 min-w-0">
              <QualityHistogram buckets={distributions[m.key]} label="distribution" />
            </div>
          )}
        </div>
      ))}

      {/* High quality summary */}
      <div className="border-t border-[#25415f] pt-2 flex items-center justify-between">
        <span className="text-[11px] text-slate-300">High Quality Precursors</span>
        <span className="metric-pill text-[11px]">
          <span className={`metric-value ${quality.pct_high_quality >= 5 ? 'status-good' : quality.pct_high_quality >= 1 ? 'status-warn' : 'status-bad'}`}>
            {fmtNum(quality.n_high_quality)}
          </span>
          <span className="text-slate-400">({fmtPct(quality.pct_high_quality)})</span>
        </span>
      </div>
    </div>
  );
}

// ─── Data Space Card ────────────────────────────────────────────────────

function DataSpaceCard({ summary }: { summary: DatasetFullSummary }) {
  const charges = Object.entries(summary.by_charge)
    .map(([k, v]) => ({ charge: parseInt(k), count: v }))
    .filter(c => c.charge > 0)
    .sort((a, b) => a.charge - b.charge);
  const totalCharged = charges.reduce((s, c) => s + c.count, 0) || 1;

  return (
    <div className="space-y-2.5">
      {/* Ranges */}
      <div className="grid grid-cols-2 gap-2">
        <div className="panel-inset p-2">
          <div className="text-[10px] text-slate-400 uppercase tracking-wide font-semibold">m/z</div>
          <div className="mono text-sm text-[#8dd6e8]">{fmtRange(summary.mz_range)}</div>
        </div>
        <div className="panel-inset p-2">
          <div className="text-[10px] text-slate-400 uppercase tracking-wide font-semibold">RT (min)</div>
          <div className="mono text-sm text-[#8dd6e8]">{fmtRange(summary.rt_range_minutes)}</div>
        </div>
        <div className="panel-inset p-2">
          <div className="text-[10px] text-slate-400 uppercase tracking-wide font-semibold">1/K0</div>
          <div className="mono text-sm text-[#8dd6e8]">{fmtRange(summary.mobility_range, 3)}</div>
        </div>
        <div className="panel-inset p-2">
          <div className="text-[10px] text-slate-400 uppercase tracking-wide font-semibold">CE (eV)</div>
          <div className="mono text-sm text-[#8dd6e8]">{fmtRange(summary.collision_energy_range)}</div>
        </div>
      </div>

      {/* Charge distribution - compact inline bars */}
      <div className="border-t border-[#25415f] pt-2">
        <div className="text-[10px] text-slate-400 uppercase tracking-wide font-semibold mb-1.5">Charge State</div>
        <div className="space-y-1">
          {charges.map((c) => {
            const pct = (c.count / totalCharged) * 100;
            return (
              <div key={c.charge} className="flex items-center gap-2">
                <span className="text-[11px] text-slate-300 w-6 text-right mono">{c.charge}+</span>
                <div className="flex-1 h-3 bg-slate-800/60 rounded overflow-hidden">
                  <div className="h-full rounded bg-cyan-600/70" style={{ width: `${pct}%` }} />
                </div>
                <span className="text-[10px] text-slate-400 w-16 text-right mono">{fmtNum(c.count)}</span>
                <span className="text-[10px] text-slate-500 w-10 text-right">{fmtPct(pct)}</span>
              </div>
            );
          })}
        </div>
      </div>
    </div>
  );
}

// ─── Raw Files Card ─────────────────────────────────────────────────────

function RawFilesCard({ summary }: { summary: DatasetFullSummary }) {
  const files = summary.raw_files;
  if (files.length === 0) return <div className="text-xs text-slate-400">No raw file data</div>;

  const sorted = [...files].sort((a, b) => b.count - a.count);
  const maxVal = Math.max(...sorted.map(f => f.count), 1);
  const totalPrecursors = sorted.reduce((s, f) => s + f.count, 0);

  return (
    <div className="space-y-2">
      <div className="flex items-center gap-2">
        <span className="metric-pill text-[10px]">{files.length} file{files.length !== 1 ? 's' : ''}</span>
        <span className="metric-pill text-[10px]">total <span className="metric-value">{fmtNum(totalPrecursors)}</span></span>
      </div>
      <div className="space-y-1">
        {sorted.map((f) => (
          <div key={f.name} className="flex items-center gap-2">
            <span className="text-[10px] text-slate-300 mono truncate flex-none max-w-[180px]" title={f.name}>
              {f.name.replace('.d', '')}
            </span>
            <div className="flex-1 h-2.5 bg-slate-800/60 rounded overflow-hidden">
              <div className="h-full rounded bg-cyan-600/60" style={{ width: `${(f.count / maxVal) * 100}%` }} />
            </div>
            <span className="text-[10px] text-slate-400 mono flex-none">{fmtNum(f.count)}</span>
          </div>
        ))}
      </div>
    </div>
  );
}

// ─── Main Summary Page ──────────────────────────────────────────────────

export default function DatasetSummaryPage({ activeDataset, onBrowse }: DatasetSummaryPageProps) {
  const { data: summary, isLoading, error } = useQuery<DatasetFullSummary>({
    queryKey: ['dataset-summary', activeDataset],
    queryFn: getDatasetFullSummary,
  });

  if (isLoading) {
    return (
      <div className="h-full flex items-center justify-center">
        <div className="metric-pill">Loading dataset summary...</div>
      </div>
    );
  }

  if (error || !summary) {
    return (
      <div className="h-full flex items-center justify-center">
        <div className="metric-pill">
          <span className="status-warn">Failed to load summary</span>
        </div>
      </div>
    );
  }

  const accession = summary.accession ?? activeDataset ?? 'Unknown';
  const dateStr = summary.generated_at
    ? new Date(summary.generated_at).toLocaleDateString('en-US', { year: 'numeric', month: 'short', day: 'numeric' })
    : null;

  return (
    <div className="h-full overflow-y-auto p-4 space-y-3 reveal-up">
      {/* Header */}
      <div className="chrome-panel overflow-hidden">
        <div className="panel-inset-head">
          <div className="flex flex-wrap items-center gap-2">
            <span className="control-label">Dataset Summary</span>
            <h2 className="text-lg font-bold tracking-tight">{accession}</h2>
          </div>
          <button onClick={onBrowse} className="btn-primary">
            Browse Precursors &rarr;
          </button>
        </div>
        <div className="p-3">
          <div className="toolbar-lane">
            <span className="metric-pill text-sm">
              Total Precursors <span className="metric-value text-base">{fmtNum(summary.n_total_precursors)}</span>
            </span>
            {summary.organism && (
              <span className="metric-pill">
                <span className="subtle-label text-[0.6rem]">Organism</span>
                <span className="metric-value">{summary.organism}</span>
              </span>
            )}
            {summary.pipeline_version && (
              <span className="metric-pill">
                <span className="subtle-label text-[0.6rem]">Version</span>
                <span className="metric-value">v{summary.pipeline_version}</span>
              </span>
            )}
            {dateStr && (
              <span className="metric-pill">
                <span className="subtle-label text-[0.6rem]">Generated</span>
                <span className="metric-value">{dateStr}</span>
              </span>
            )}
            {summary.study_id && (
              <span className="metric-pill">
                <span className="subtle-label text-[0.6rem]">Study</span>
                <span className="metric-value">{summary.study_id}</span>
              </span>
            )}
            {summary.group_id && (
              <span className="metric-pill">
                <span className="subtle-label text-[0.6rem]">Group</span>
                <span className="metric-value">{summary.group_id}</span>
              </span>
            )}
          </div>
        </div>
      </div>

      {/* 2x2 Grid */}
      <div className="grid grid-cols-1 xl:grid-cols-2 gap-3">
        {/* Engine Agreement + Per-Engine + Overlap */}
        <div className="chrome-panel overflow-hidden">
          <div className="panel-inset-head">
            <div className="control-label">Engine Agreement & Identifications</div>
          </div>
          <div className="p-3">
            <EngineAgreementCard
              quality={summary.quality_summary}
              engines={summary.n_per_engine}
              unidentified={summary.n_unidentified}
              overlap={summary.overlap_stats}
              uniquePeptides={summary.n_unique_peptides}
            />
          </div>
        </div>

        {/* Quality Metrics + Distributions */}
        <div className="chrome-panel overflow-hidden">
          <div className="panel-inset-head">
            <div className="control-label">Quality Metrics</div>
          </div>
          <div className="p-3">
            <QualityMetricsCard
              quality={summary.quality_summary}
              distributions={summary.quality_distributions}
            />
          </div>
        </div>

        {/* Data Space: ranges + charge */}
        <div className="chrome-panel overflow-hidden">
          <div className="panel-inset-head">
            <div className="control-label">Data Space</div>
          </div>
          <div className="p-3">
            <DataSpaceCard summary={summary} />
          </div>
        </div>

        {/* Raw Files */}
        <div className="chrome-panel overflow-hidden">
          <div className="panel-inset-head">
            <div className="control-label">Raw Files</div>
          </div>
          <div className="p-3">
            <RawFilesCard summary={summary} />
          </div>
        </div>
      </div>
    </div>
  );
}
