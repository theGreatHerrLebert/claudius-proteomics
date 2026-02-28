import { useRef, useEffect } from 'react';

/* ── types ── */

interface IonCloudData {
  mz: number[];
  inv_mob: number[];
  rt: number[];
  charge: number[];
}

type Projection = 'mz_vs_im' | 'rt_vs_mz';

/* ── constants ── */

const CHARGE_COLORS: Record<number, [number, number, number]> = {
  2: [82, 186, 216],   // teal
  3: [100, 160, 220],  // soft blue
  4: [130, 200, 230],  // light blue
};
const CHARGES = [2, 3, 4];

const DOT_SIZE = 1.5;
const DOT_ALPHA = 0.55;
const PAD = 0.05;

// Entrance animation
const ANIM_DURATION = 8_000;
const PARTICLE_DURATION = 1_200;
const MAX_STAGGER = 8_000;

// Post-settle
const SETTLED_OPACITY = 0.28;
const FADE_DURATION = 6_000;   // smooth fade from full → settled opacity
const BREATHE_AMPLITUDE = 0.04;
const BREATHE_PERIOD = 6_000;

// Marginals
const N_BINS = 120;
const MARGINAL_X_HEIGHT = 50;  // px – bottom marginal
const MARGINAL_Y_WIDTH = 40;   // px – left marginal
const MARGINAL_ALPHA = 0.35;

/* ── helpers ── */

function bounds(arr: number[]): [number, number] {
  let lo = arr[0];
  let hi = arr[0];
  for (let i = 1; i < arr.length; i++) {
    if (arr[i] < lo) lo = arr[i];
    if (arr[i] > hi) hi = arr[i];
  }
  return [lo, hi];
}

function easeOutCubic(t: number): number {
  return 1 - (1 - t) ** 3;
}

function hash(i: number): number {
  let x = ((i + 1) * 2654435761) >>> 0;
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = (x >> 16) ^ x;
  return (x >>> 0) / 0xffffffff;
}

function pickProjection(): Projection {
  return Math.random() < 0.5 ? 'mz_vs_im' : 'rt_vs_mz';
}

interface ProjectionResult {
  endXNorm: Float32Array;
  endYNorm: Float32Array;
  xLabel: string;
  yLabel: string;
  xLo: number;
  xHi: number;
  yLo: number;
  yHi: number;
}

function projectData(data: IonCloudData, proj: Projection): ProjectionResult {
  const n = data.mz.length;
  let xRaw: number[], yRaw: number[];
  let xLabel: string, yLabel: string;

  if (proj === 'mz_vs_im') {
    xRaw = data.mz;
    yRaw = data.inv_mob;
    xLabel = 'm/z';
    yLabel = '1/K\u2080';
  } else {
    xRaw = data.rt;
    yRaw = data.mz;
    xLabel = 'RT (norm.)';
    yLabel = 'm/z';
  }

  const [xLo, xHi] = bounds(xRaw);
  const [yLo, yHi] = bounds(yRaw);
  const xRange = xHi - xLo || 1;
  const yRange = yHi - yLo || 1;

  const endXNorm = new Float32Array(n);
  const endYNorm = new Float32Array(n);

  for (let i = 0; i < n; i++) {
    endXNorm[i] = (xRaw[i] - xLo) / xRange;
    endYNorm[i] = (yHi - yRaw[i]) / yRange; // invert y
  }

  return { endXNorm, endYNorm, xLabel, yLabel, xLo, xHi, yLo, yHi };
}

/** Generate nice round tick values between lo and hi */
function niceTicks(lo: number, hi: number, maxTicks = 5): number[] {
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
function formatTick(v: number): string {
  if (Math.abs(v) >= 100) return Math.round(v).toString();
  if (Math.abs(v) >= 1) return v.toFixed(1).replace(/\.0$/, '');
  return v.toFixed(2);
}

const AXIS_COLOR = 'rgba(127, 148, 184, 0.35)';  // --text-3 at low alpha
const TICK_COLOR = 'rgba(127, 148, 184, 0.28)';
const LABEL_COLOR = 'rgba(175, 192, 219, 0.5)';   // --text-2 at low alpha
const AXIS_FONT = '10px "IBM Plex Mono", monospace';
const LABEL_FONT = '11px "IBM Plex Sans", sans-serif';

/** Draw subtle axes with ticks and labels */
function drawAxes(
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
  // x-axis (bottom)
  ctx.moveTo(xOff, yEnd);
  ctx.lineTo(xEnd, yEnd);
  // y-axis (left)
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
    // Tick mark
    ctx.beginPath();
    ctx.strokeStyle = AXIS_COLOR;
    ctx.moveTo(x, yEnd);
    ctx.lineTo(x, yEnd + 4);
    ctx.stroke();
    // Label
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
    const frac = (proj.yHi - v) / yRange; // inverted
    const y = yOff + frac * (yEnd - yOff);
    // Tick mark
    ctx.beginPath();
    ctx.strokeStyle = AXIS_COLOR;
    ctx.moveTo(xOff, y);
    ctx.lineTo(xOff - 4, y);
    ctx.stroke();
    // Label
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
function drawMarginalCurve(
  ctx: CanvasRenderingContext2D,
  bins: Float32Array,    // counts per bin (already normalized to 0..1)
  color: [number, number, number],
  axis: 'x' | 'y',
  areaStart: number,     // baseline coordinate
  areaEnd: number,       // max coordinate (peak extends toward here)
  rangeStart: number,    // start along the bins axis
  rangeEnd: number,      // end along the bins axis
  stackOffset: Float32Array | null, // stacked offset from previous charge
) {
  const nBins = bins.length;
  const step = (rangeEnd - rangeStart) / nBins;

  ctx.beginPath();

  if (axis === 'x') {
    // Bottom marginal: bins along x, height grows upward (areaStart=bottom, areaEnd=top)
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
    // Left marginal: bins along y, width grows rightward (areaStart=left, areaEnd=right)
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

/* ── component ── */

export default function IonCloudCanvas({ className = '' }: { className?: string }) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const rafRef = useRef(0);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    let cancelled = false;
    const projection = pickProjection();

    fetch('/data/ion_cloud.json')
      .then((r) => r.json())
      .then((json: IonCloudData) => {
        if (cancelled) return;
        startRendering(json, projection);
      })
      .catch(() => {});

    function startRendering(data: IonCloudData, proj: Projection) {
      const ctx = canvas!.getContext('2d');
      if (!ctx) return;

      const prefersReduced = window.matchMedia(
        '(prefers-reduced-motion: reduce)',
      ).matches;

      const n = data.mz.length;
      const projResult = projectData(data, proj);
      const { endXNorm, endYNorm } = projResult;

      // Per-particle setup
      const startXNorm = new Float32Array(n);
      const startYNorm = new Float32Array(n);
      const delays = new Float32Array(n);
      const colors: string[] = new Array(n);

      // Marginal bin assignments per particle
      const xBin = new Uint8Array(n);
      const yBin = new Uint8Array(n);

      // Final histogram counts per (charge, bin) — for stacked marginals
      const finalXHist: Record<number, Float32Array> = {};
      const finalYHist: Record<number, Float32Array> = {};
      for (const c of CHARGES) {
        finalXHist[c] = new Float32Array(N_BINS);
        finalYHist[c] = new Float32Array(N_BINS);
      }

      for (let i = 0; i < n; i++) {
        startXNorm[i] = -0.05 - hash(i) * 0.15;
        startYNorm[i] = endYNorm[i];
        delays[i] = endXNorm[i] * MAX_STAGGER;

        const c = data.charge[i];
        const [r, g, b] = CHARGE_COLORS[c] ?? CHARGE_COLORS[2];
        colors[i] = `rgba(${r},${g},${b},${DOT_ALPHA})`;

        // Bin assignment
        const bx = Math.min(Math.floor(endXNorm[i] * N_BINS), N_BINS - 1);
        const by = Math.min(Math.floor(endYNorm[i] * N_BINS), N_BINS - 1);
        xBin[i] = bx;
        yBin[i] = by;
        finalXHist[c][bx]++;
        finalYHist[c][by]++;
      }

      // Normalize final histograms: find global max across all charges (stacked)
      const stackedX = new Float32Array(N_BINS);
      const stackedY = new Float32Array(N_BINS);
      for (const c of CHARGES) {
        for (let b = 0; b < N_BINS; b++) {
          stackedX[b] += finalXHist[c][b];
          stackedY[b] += finalYHist[c][b];
        }
      }
      let maxX = 0, maxY = 0;
      for (let b = 0; b < N_BINS; b++) {
        if (stackedX[b] > maxX) maxX = stackedX[b];
        if (stackedY[b] > maxY) maxY = stackedY[b];
      }
      // Normalize to [0..1]
      for (const c of CHARGES) {
        for (let b = 0; b < N_BINS; b++) {
          finalXHist[c][b] /= maxX || 1;
          finalYHist[c][b] /= maxY || 1;
        }
      }

      // Live histogram accumulators (raw counts, normalized each frame)
      const liveXHist: Record<number, Float32Array> = {};
      const liveYHist: Record<number, Float32Array> = {};
      for (const c of CHARGES) {
        liveXHist[c] = new Float32Array(N_BINS);
        liveYHist[c] = new Float32Array(N_BINS);
      }
      const arrived = new Uint8Array(n); // 0 = not arrived, 1 = counted

      // Static settled canvas — bakes in the final settled alpha per dot
      const FINAL_DOT_ALPHA = DOT_ALPHA * SETTLED_OPACITY;
      let staticCanvas: HTMLCanvasElement | null = null;

      // Pre-compute settled colors (with final alpha baked in)
      const settledColors: string[] = new Array(n);
      for (let i = 0; i < n; i++) {
        const c = data.charge[i];
        const [r, g, b] = CHARGE_COLORS[c] ?? CHARGE_COLORS[2];
        settledColors[i] = `rgba(${r},${g},${b},${FINAL_DOT_ALPHA})`;
      }

      function buildStatic(w: number, h: number, dpr: number) {
        const off = document.createElement('canvas');
        off.width = w * dpr;
        off.height = h * dpr;
        const octx = off.getContext('2d')!;
        octx.scale(dpr, dpr);
        const xOff = PAD * w;
        const yOff = PAD * h;
        const xScale = w * (1 - 2 * PAD);
        const yScale = h * (1 - 2 * PAD);
        for (let i = 0; i < n; i++) {
          octx.fillStyle = settledColors[i];
          octx.fillRect(
            xOff + endXNorm[i] * xScale,
            yOff + endYNorm[i] * yScale,
            DOT_SIZE,
            DOT_SIZE,
          );
        }
        return off;
      }

      let w = 0, h = 0, dpr = 1;

      function resize() {
        const rect = canvas!.getBoundingClientRect();
        dpr = Math.min(window.devicePixelRatio || 1, 2);
        w = rect.width;
        h = rect.height;
        canvas!.width = w * dpr;
        canvas!.height = h * dpr;
        ctx!.setTransform(dpr, 0, 0, dpr, 0, 0);
        staticCanvas = buildStatic(w, h, dpr);
      }

      resize();

      /** Draw marginals from the given histogram maps */
      function drawMarginals(
        xHist: Record<number, Float32Array>,
        yHist: Record<number, Float32Array>,
      ) {
        const xOff = PAD * w;
        const yOff = PAD * h;
        const xEnd = xOff + w * (1 - 2 * PAD);
        const yEnd = yOff + h * (1 - 2 * PAD);

        // Bottom marginal (x axis): baseline on axis line, peaks grow upward into plot
        let xStack: Float32Array | null = null;
        for (const c of CHARGES) {
          drawMarginalCurve(
            ctx!, xHist[c], CHARGE_COLORS[c],
            'x',
            yEnd, yEnd - MARGINAL_X_HEIGHT,
            xOff, xEnd,
            xStack,
          );
          if (!xStack) xStack = new Float32Array(N_BINS);
          for (let b = 0; b < N_BINS; b++) xStack[b] += xHist[c][b];
        }

        // Left marginal (y axis): baseline on axis line, peaks grow rightward into plot
        let yStack: Float32Array | null = null;
        for (const c of CHARGES) {
          drawMarginalCurve(
            ctx!, yHist[c], CHARGE_COLORS[c],
            'y',
            xOff, xOff + MARGINAL_Y_WIDTH,
            yOff, yEnd,
            yStack,
          );
          if (!yStack) yStack = new Float32Array(N_BINS);
          for (let b = 0; b < N_BINS; b++) yStack[b] += yHist[c][b];
        }
      }

      if (prefersReduced) {
        if (staticCanvas) {
          ctx!.globalAlpha = SETTLED_OPACITY;
          ctx!.drawImage(staticCanvas, 0, 0, w, h);
          ctx!.globalAlpha = 1;
          drawMarginals(finalXHist, finalYHist);
          drawAxes(ctx!, w, h, projResult);
        }
        return;
      }

      let resizeTimer = 0;
      const ro = new ResizeObserver(() => {
        clearTimeout(resizeTimer);
        resizeTimer = window.setTimeout(resize, 120);
      });
      ro.observe(canvas!);

      const t0 = performance.now();
      let settled = false;

      function frame(now: number) {
        if (cancelled) return;
        if (w === 0) {
          rafRef.current = requestAnimationFrame(frame);
          return;
        }

        const elapsed = now - t0;
        ctx!.clearRect(0, 0, w, h);

        const animEnd = ANIM_DURATION + PARTICLE_DURATION;
        const allDone = elapsed >= animEnd;
        const fading = allDone && elapsed < animEnd + FADE_DURATION;
        const fullySettled = elapsed >= animEnd + FADE_DURATION;

        // Fade progress: 0 = animation just ended, 1 = fully settled
        const fadeFrac = fading
          ? (elapsed - animEnd) / FADE_DURATION
          : fullySettled ? 1 : 0;

        // Per-dot alpha at Phase 1 end = SETTLED_OPACITY (0.28)
        // Final settled per-dot alpha = FINAL_DOT_ALPHA (0.154)
        // Marginals at Phase 1 end = full opacity (1.0)

        const xOff = PAD * w;
        const yOff = PAD * h;
        const xScale = w * (1 - 2 * PAD);
        const yScale = h * (1 - 2 * PAD);

        if (!allDone) {
          // ── Phase 1: Entrance ──
          ctx!.globalAlpha = 1;

          for (let i = 0; i < n; i++) {
            const localT = (elapsed - delays[i]) / PARTICLE_DURATION;
            if (localT < 0) continue;

            // Track arrival for marginals (count once when > 80% settled)
            if (!arrived[i] && localT >= 0.8) {
              arrived[i] = 1;
              const c = data.charge[i];
              liveXHist[c][xBin[i]]++;
              liveYHist[c][yBin[i]]++;
            }

            const progress = localT >= 1 ? 1 : easeOutCubic(localT);

            const nx = startXNorm[i] + (endXNorm[i] - startXNorm[i]) * progress;
            const ny = startYNorm[i] + (endYNorm[i] - startYNorm[i]) * progress;

            const x = xOff + nx * xScale;
            const y = yOff + ny * yScale;

            const fadeIn = Math.min(localT * 3, 1);
            const settledAlpha = SETTLED_OPACITY;
            const flyAlpha = DOT_ALPHA * fadeIn;
            const alpha = settledAlpha + (flyAlpha - settledAlpha) * (1 - progress);
            const [r, g, b] = CHARGE_COLORS[data.charge[i]] ?? CHARGE_COLORS[2];
            ctx!.fillStyle = `rgba(${r},${g},${b},${alpha})`;
            ctx!.fillRect(x, y, DOT_SIZE, DOT_SIZE);
          }

          // Draw live marginals (normalize by same global max)
          const liveNormX: Record<number, Float32Array> = {};
          const liveNormY: Record<number, Float32Array> = {};
          for (const c of CHARGES) {
            liveNormX[c] = new Float32Array(N_BINS);
            liveNormY[c] = new Float32Array(N_BINS);
            for (let b = 0; b < N_BINS; b++) {
              liveNormX[c][b] = liveXHist[c][b] / (maxX || 1);
              liveNormY[c][b] = liveYHist[c][b] / (maxY || 1);
            }
          }
          drawMarginals(liveNormX, liveNormY);
          drawAxes(ctx!, w, h, projResult);

        } else if (fading) {
          // ── Phase 2a: Fade – draw individual particles (same method as Phase 1) ──
          const easedFade = easeOutCubic(fadeFrac);
          ctx!.globalAlpha = 1;

          // Per-dot alpha lerps from SETTLED_OPACITY → FINAL_DOT_ALPHA
          for (let i = 0; i < n; i++) {
            const alpha = SETTLED_OPACITY + (FINAL_DOT_ALPHA - SETTLED_OPACITY) * easedFade;
            const [r, g, b] = CHARGE_COLORS[data.charge[i]] ?? CHARGE_COLORS[2];
            ctx!.fillStyle = `rgba(${r},${g},${b},${alpha})`;
            ctx!.fillRect(
              xOff + endXNorm[i] * xScale,
              yOff + endYNorm[i] * yScale,
              DOT_SIZE,
              DOT_SIZE,
            );
          }

          // Marginals fade from 1.0 → SETTLED_OPACITY
          const overlayMul = 1 - (1 - SETTLED_OPACITY) * easedFade;
          ctx!.globalAlpha = overlayMul;
          drawMarginals(finalXHist, finalYHist);
          drawAxes(ctx!, w, h, projResult);

        } else {
          // ── Phase 2b: Fully settled – use pre-rendered static canvas ──
          const breathe =
            Math.sin(((now % BREATHE_PERIOD) / BREATHE_PERIOD) * Math.PI * 2) * BREATHE_AMPLITUDE;
          ctx!.globalAlpha = 1 + breathe;
          ctx!.drawImage(staticCanvas!, 0, 0, w, h);

          ctx!.globalAlpha = SETTLED_OPACITY + breathe;
          drawMarginals(finalXHist, finalYHist);
          drawAxes(ctx!, w, h, projResult);
        }

        rafRef.current = requestAnimationFrame(frame);
      }

      rafRef.current = requestAnimationFrame(frame);

      return () => {
        ro.disconnect();
      };
    }

    return () => {
      cancelled = true;
      cancelAnimationFrame(rafRef.current);
    };
  }, []);

  return (
    <canvas
      ref={canvasRef}
      className={className}
      style={{ display: 'block', width: '100%', height: '100%' }}
    />
  );
}
