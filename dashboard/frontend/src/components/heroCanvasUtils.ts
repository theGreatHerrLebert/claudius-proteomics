/* ── Shared canvas utilities for hero background visualizations ── */

/* ── types ── */

export interface ProjectionResult {
  endXNorm: Float32Array;
  endYNorm: Float32Array;
  xLabel: string;
  yLabel: string;
  xLo: number;
  xHi: number;
  yLo: number;
  yHi: number;
}

/* ── constants ── */

export const PAD = 0.05;

export const AXIS_COLOR = 'rgba(127, 148, 184, 0.35)';
export const TICK_COLOR = 'rgba(127, 148, 184, 0.28)';
export const LABEL_COLOR = 'rgba(175, 192, 219, 0.5)';
export const AXIS_FONT = '10px "IBM Plex Mono", monospace';
export const LABEL_FONT = '11px "IBM Plex Sans", sans-serif';
export const MARGINAL_ALPHA = 0.35;

/* ── helpers ── */

export function bounds(arr: number[]): [number, number] {
  let lo = arr[0];
  let hi = arr[0];
  for (let i = 1; i < arr.length; i++) {
    if (arr[i] < lo) lo = arr[i];
    if (arr[i] > hi) hi = arr[i];
  }
  return [lo, hi];
}

export function easeOutCubic(t: number): number {
  return 1 - (1 - t) ** 3;
}

export function hash(i: number): number {
  let x = ((i + 1) * 2654435761) >>> 0;
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = (x >> 16) ^ x;
  return (x >>> 0) / 0xffffffff;
}

/** Generate nice round tick values between lo and hi */
export function niceTicks(lo: number, hi: number, maxTicks = 5): number[] {
  const range = hi - lo;
  const rough = range / maxTicks;
  const mag = Math.pow(10, Math.floor(Math.log10(rough)));
  let step: number;
  const norm = rough / mag;
  if (norm < 1.5) step = mag;
  else if (norm < 3.5) step = 2 * mag;
  else if (norm < 7.5) step = 5 * mag;
  else step = 10 * mag;

  const ticks: number[] = [];
  const start = Math.ceil(lo / step) * step;
  for (let v = start; v <= hi; v += step) {
    ticks.push(v);
  }
  return ticks;
}

/** Format tick value — drop unnecessary decimals */
export function formatTick(v: number): string {
  if (Math.abs(v) >= 100) return Math.round(v).toString();
  if (Math.abs(v) >= 1) return v.toFixed(1).replace(/\.0$/, '');
  return v.toFixed(2);
}

/** Draw subtle axes with ticks and labels */
export function drawAxes(
  ctx: CanvasRenderingContext2D,
  w: number,
  h: number,
  proj: ProjectionResult,
) {
  const xOff = PAD * w;
  const yOff = PAD * h;
  const xEnd = xOff + w * (1 - 2 * PAD);
  const yEnd = yOff + h * (1 - 2 * PAD);

  ctx.save();
  ctx.globalAlpha = 1;

  // Axis lines
  ctx.strokeStyle = AXIS_COLOR;
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(xOff, yEnd);
  ctx.lineTo(xEnd, yEnd);
  ctx.moveTo(xOff, yOff);
  ctx.lineTo(xOff, yEnd);
  ctx.stroke();

  // X-axis ticks + labels
  const xTicks = niceTicks(proj.xLo, proj.xHi, 5);
  const xRange = proj.xHi - proj.xLo || 1;
  ctx.font = AXIS_FONT;
  ctx.fillStyle = TICK_COLOR;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';

  for (const v of xTicks) {
    const frac = (v - proj.xLo) / xRange;
    const x = xOff + frac * (xEnd - xOff);
    ctx.beginPath();
    ctx.strokeStyle = AXIS_COLOR;
    ctx.moveTo(x, yEnd);
    ctx.lineTo(x, yEnd + 4);
    ctx.stroke();
    ctx.fillText(formatTick(v), x, yEnd + 6);
  }

  // X-axis label
  ctx.font = LABEL_FONT;
  ctx.fillStyle = LABEL_COLOR;
  ctx.fillText(proj.xLabel, (xOff + xEnd) / 2, yEnd + 20);

  // Y-axis ticks + labels (inverted: yHi at top, yLo at bottom)
  const yTicks = niceTicks(proj.yLo, proj.yHi, 4);
  const yRange = proj.yHi - proj.yLo || 1;
  ctx.font = AXIS_FONT;
  ctx.fillStyle = TICK_COLOR;
  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';

  for (const v of yTicks) {
    const frac = (proj.yHi - v) / yRange;
    const y = yOff + frac * (yEnd - yOff);
    ctx.beginPath();
    ctx.strokeStyle = AXIS_COLOR;
    ctx.moveTo(xOff, y);
    ctx.lineTo(xOff - 4, y);
    ctx.stroke();
    ctx.fillText(formatTick(v), xOff - 7, y);
  }

  // Y-axis label (rotated)
  ctx.save();
  ctx.font = LABEL_FONT;
  ctx.fillStyle = LABEL_COLOR;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'bottom';
  const yMid = (yOff + yEnd) / 2;
  ctx.translate(xOff - 30, yMid);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText(proj.yLabel, 0, 0);
  ctx.restore();

  ctx.restore();
}

/** Draw a smooth filled area curve from bin counts */
export function drawMarginalCurve(
  ctx: CanvasRenderingContext2D,
  bins: Float32Array,
  color: [number, number, number],
  axis: 'x' | 'y',
  areaStart: number,
  areaEnd: number,
  rangeStart: number,
  rangeEnd: number,
  stackOffset: Float32Array | null,
) {
  const nBins = bins.length;
  const step = (rangeEnd - rangeStart) / nBins;

  ctx.beginPath();

  if (axis === 'x') {
    const baseline = areaStart;
    ctx.moveTo(rangeStart, baseline);
    for (let b = 0; b < nBins; b++) {
      const cx = rangeStart + (b + 0.5) * step;
      const offset = stackOffset ? stackOffset[b] : 0;
      const cy = baseline - (bins[b] + offset) * (baseline - areaEnd);
      ctx.lineTo(cx, cy);
    }
    ctx.lineTo(rangeEnd, baseline);
    ctx.closePath();
  } else {
    const baseline = areaStart;
    ctx.moveTo(baseline, rangeStart);
    for (let b = 0; b < nBins; b++) {
      const cy = rangeStart + (b + 0.5) * step;
      const offset = stackOffset ? stackOffset[b] : 0;
      const cx = baseline + (bins[b] + offset) * (areaEnd - baseline);
      ctx.lineTo(cx, cy);
    }
    ctx.lineTo(baseline, rangeEnd);
    ctx.closePath();
  }

  ctx.fillStyle = `rgba(${color[0]},${color[1]},${color[2]},${MARGINAL_ALPHA})`;
  ctx.fill();
}
