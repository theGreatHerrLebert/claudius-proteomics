import { useRef, useEffect } from 'react';
import {
  easeOutCubic, drawAxes, drawMarginalCurve,
  PAD,
  type ProjectionResult,
} from './heroCanvasUtils';

/* ── types ── */

interface HeatmapData {
  grid: {
    mz_lo: number; mz_hi: number; n_mz_bins: number;
    mob_lo: number; mob_hi: number; n_mob_bins: number;
    counts: number[][];
    mean_engines: number[][];
  };
  marginal_mz: Record<string, number[]>;
  marginal_mob: Record<string, number[]>;
}

/* ── constants ── */

const ENGINE_COLORS: Record<number, [number, number, number]> = {
  3: [16, 185, 129],
  2: [82, 186, 216],
  1: [74, 85, 104],
};

const SWEEP_DURATION = 5_000;
const MARGINAL_DELAY = 2_000;      // Start marginals during sweep
const SETTLE_DELAY = 1_500;        // After sweep ends
const FADE_DURATION = 4_000;
const SETTLED_OPACITY = 0.18;

const MARGINAL_X_HEIGHT = 50;
const MARGINAL_Y_WIDTH = 40;
const N_MARGINAL_BINS = 120;

/* ── color helpers ── */

function engineCountToColor(mean: number): [number, number, number] {
  if (mean <= 1.0) return [42, 63, 85];
  if (mean <= 2.0) {
    const t = mean - 1.0;
    return [
      Math.round(42 + t * (58 - 42)),
      Math.round(63 + t * (143 - 63)),
      Math.round(85 + t * (168 - 85)),
    ];
  }
  const t = Math.min(mean - 2.0, 1.0);
  return [
    Math.round(58 + t * (16 - 58)),
    Math.round(143 + t * (185 - 143)),
    Math.round(168 + t * (129 - 168)),
  ];
}

/* ── component ── */

export default function ConsensusHeatmapCanvas({ className = '' }: { className?: string }) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const rafRef = useRef(0);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    let cancelled = false;

    fetch('/data/consensus_heatmap.json')
      .then((r) => r.json())
      .then((json: HeatmapData) => {
        if (cancelled) return;
        startRendering(json);
      })
      .catch(() => {});

    function startRendering(data: HeatmapData) {
      const ctx = canvas!.getContext('2d');
      if (!ctx) return;

      const prefersReduced = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

      const { grid, marginal_mz, marginal_mob } = data;
      const nMz = grid.n_mz_bins;
      const nMob = grid.n_mob_bins;

      // Find max count for alpha scaling
      let maxCount = 0;
      for (let mx = 0; mx < nMz; mx++) {
        for (let my = 0; my < nMob; my++) {
          if (grid.counts[mx][my] > maxCount) maxCount = grid.counts[mx][my];
        }
      }

      // Pre-compute diagonal distance for each bin (for sweep timing)
      const maxDiag = nMz + nMob;
      const binDiag = new Float32Array(nMz * nMob);
      for (let mx = 0; mx < nMz; mx++) {
        for (let my = 0; my < nMob; my++) {
          binDiag[mx * nMob + my] = (mx + my) / maxDiag;
        }
      }

      // Pre-compute bin colors
      const binColors: Array<[number, number, number]> = new Array(nMz * nMob);
      for (let mx = 0; mx < nMz; mx++) {
        for (let my = 0; my < nMob; my++) {
          const me = grid.mean_engines[mx][my];
          binColors[mx * nMob + my] = me > 0 ? engineCountToColor(me) : [10, 21, 32];
        }
      }

      // Prepare marginal histograms (resample to N_MARGINAL_BINS)
      const finalMzHist: Record<number, Float32Array> = {};
      const finalMobHist: Record<number, Float32Array> = {};

      for (const k of [1, 2, 3]) {
        const srcMz = marginal_mz[String(k)] || [];
        const srcMob = marginal_mob[String(k)] || [];
        finalMzHist[k] = new Float32Array(N_MARGINAL_BINS);
        finalMobHist[k] = new Float32Array(N_MARGINAL_BINS);

        // Resample from grid bins to marginal bins
        for (let b = 0; b < N_MARGINAL_BINS; b++) {
          const srcIdx = Math.floor(b / N_MARGINAL_BINS * srcMz.length);
          finalMzHist[k][b] = srcMz[srcIdx] || 0;
          const srcIdxM = Math.floor(b / N_MARGINAL_BINS * srcMob.length);
          finalMobHist[k][b] = srcMob[srcIdxM] || 0;
        }
      }

      // Normalize marginals
      let maxMzStack = 0, maxMobStack = 0;
      for (let b = 0; b < N_MARGINAL_BINS; b++) {
        let sx = 0, sy = 0;
        for (const k of [1, 2, 3]) {
          sx += finalMzHist[k][b];
          sy += finalMobHist[k][b];
        }
        if (sx > maxMzStack) maxMzStack = sx;
        if (sy > maxMobStack) maxMobStack = sy;
      }
      for (const k of [1, 2, 3]) {
        for (let b = 0; b < N_MARGINAL_BINS; b++) {
          finalMzHist[k][b] /= maxMzStack || 1;
          finalMobHist[k][b] /= maxMobStack || 1;
        }
      }

      const projResult: ProjectionResult = {
        endXNorm: new Float32Array(0),
        endYNorm: new Float32Array(0),
        xLabel: 'm/z',
        yLabel: '1/K\u2080',
        xLo: grid.mz_lo,
        xHi: grid.mz_hi,
        yLo: grid.mob_lo,
        yHi: grid.mob_hi,
      };

      let w = 0, h = 0, dpr = 1;

      function resize() {
        const rect = canvas!.getBoundingClientRect();
        dpr = Math.min(window.devicePixelRatio || 1, 2);
        w = rect.width;
        h = rect.height;
        canvas!.width = w * dpr;
        canvas!.height = h * dpr;
        ctx!.setTransform(dpr, 0, 0, dpr, 0, 0);
      }

      resize();

      function drawHeatmap(sweepFrac: number) {
        const xOff = PAD * w;
        const yOff = PAD * h;
        const plotW = w * (1 - 2 * PAD);
        const plotH = h * (1 - 2 * PAD);
        const binW = plotW / nMz;
        const binH = plotH / nMob;

        for (let mx = 0; mx < nMz; mx++) {
          for (let my = 0; my < nMob; my++) {
            const count = grid.counts[mx][my];
            if (count === 0) continue;

            const diagFrac = binDiag[mx * nMob + my];
            const localT = (sweepFrac - diagFrac * 0.7) / 0.3;
            if (localT <= 0) continue;

            const reveal = Math.min(localT, 1);
            const alpha = easeOutCubic(reveal) * Math.min(count / maxCount * 3, 0.8);

            const [r, g, b] = binColors[mx * nMob + my];
            ctx!.fillStyle = `rgba(${r},${g},${b},${alpha})`;
            // Inverted y: high mobility at top
            ctx!.fillRect(
              xOff + mx * binW,
              yOff + (nMob - 1 - my) * binH,
              binW + 0.5, binH + 0.5,
            );
          }
        }
      }

      function drawMarginals(alpha: number) {
        const xOff = PAD * w;
        const yOff = PAD * h;
        const xEnd = xOff + w * (1 - 2 * PAD);
        const yEnd = yOff + h * (1 - 2 * PAD);

        ctx!.globalAlpha = alpha;

        let xStack: Float32Array | null = null;
        for (const k of [1, 2, 3]) {
          drawMarginalCurve(ctx!, finalMzHist[k], ENGINE_COLORS[k], 'x', yEnd, yEnd - MARGINAL_X_HEIGHT, xOff, xEnd, xStack);
          if (!xStack) xStack = new Float32Array(N_MARGINAL_BINS);
          for (let b = 0; b < N_MARGINAL_BINS; b++) xStack[b] += finalMzHist[k][b];
        }

        let yStack: Float32Array | null = null;
        for (const k of [1, 2, 3]) {
          drawMarginalCurve(ctx!, finalMobHist[k], ENGINE_COLORS[k], 'y', xOff, xOff + MARGINAL_Y_WIDTH, yOff, yEnd, yStack);
          if (!yStack) yStack = new Float32Array(N_MARGINAL_BINS);
          for (let b = 0; b < N_MARGINAL_BINS; b++) yStack[b] += finalMobHist[k][b];
        }

        ctx!.globalAlpha = 1;
      }

      function paintSettled() {
        ctx!.clearRect(0, 0, w, h);
        ctx!.globalAlpha = SETTLED_OPACITY;
        drawHeatmap(1);
        drawMarginals(SETTLED_OPACITY);
        drawAxes(ctx!, w, h, projResult);
        ctx!.globalAlpha = 1;
      }

      if (prefersReduced) {
        paintSettled();
        return;
      }

      let resizeTimer = 0;
      let animDone = false;
      const ro = new ResizeObserver(() => {
        clearTimeout(resizeTimer);
        resizeTimer = window.setTimeout(() => {
          resize();
          if (animDone) paintSettled();
        }, 120);
      });
      ro.observe(canvas!);

      const t0 = performance.now();
      const SETTLE_START = SWEEP_DURATION + SETTLE_DELAY;
      const TOTAL_ANIM = SETTLE_START + FADE_DURATION;

      function frame(now: number) {
        if (cancelled) return;
        if (w === 0) { rafRef.current = requestAnimationFrame(frame); return; }

        const elapsed = now - t0;
        ctx!.clearRect(0, 0, w, h);

        if (elapsed < SETTLE_START) {
          // Sweep + marginals
          const sweepFrac = Math.min(elapsed / SWEEP_DURATION, 1);
          drawHeatmap(sweepFrac);

          // Marginals fade in during sweep
          const margAlpha = Math.max(0, Math.min((elapsed - MARGINAL_DELAY) / (SWEEP_DURATION - MARGINAL_DELAY), 1));
          if (margAlpha > 0) drawMarginals(margAlpha);

          drawAxes(ctx!, w, h, projResult);

        } else if (elapsed < TOTAL_ANIM) {
          // Fade to settled
          const fadeFrac = easeOutCubic((elapsed - SETTLE_START) / FADE_DURATION);
          const alpha = 1 + (SETTLED_OPACITY - 1) * fadeFrac;
          ctx!.globalAlpha = alpha;
          drawHeatmap(1);
          drawMarginals(alpha);
          drawAxes(ctx!, w, h, projResult);
          ctx!.globalAlpha = 1;

        } else {
          animDone = true;
          paintSettled();
          return;
        }

        rafRef.current = requestAnimationFrame(frame);
      }

      rafRef.current = requestAnimationFrame(frame);

      return () => { ro.disconnect(); };
    }

    return () => {
      cancelled = true;
      cancelAnimationFrame(rafRef.current);
    };
  }, []);

  return (
    <canvas ref={canvasRef} className={className} style={{ display: 'block', width: '100%', height: '100%' }} />
  );
}
