import { useRef, useEffect } from 'react';
import {
  bounds, easeOutCubic, hash, drawAxes, drawMarginalCurve,
  PAD,
  type ProjectionResult,
} from './heroCanvasUtils';

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

// Entrance animation
const ANIM_DURATION = 8_000;
const PARTICLE_DURATION = 1_200;
const MAX_STAGGER = 8_000;

// Post-settle
const SETTLED_OPACITY = 0.28;
const FADE_DURATION = 6_000;

// Marginals
const N_BINS = 120;
const MARGINAL_X_HEIGHT = 50;
const MARGINAL_Y_WIDTH = 40;

/* ── helpers ── */

function pickProjection(): Projection {
  return Math.random() < 0.5 ? 'mz_vs_im' : 'rt_vs_mz';
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

      /** Paint the final settled state (scatter + marginals + axes) */
      function paintSettled() {
        ctx!.clearRect(0, 0, w, h);
        ctx!.globalAlpha = 1;
        ctx!.drawImage(staticCanvas!, 0, 0, w, h);
        ctx!.globalAlpha = SETTLED_OPACITY;
        drawMarginals(finalXHist, finalYHist);
        drawAxes(ctx!, w, h, projResult);
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
      const ENTRANCE_END = ANIM_DURATION + PARTICLE_DURATION;
      const FADE_END = ENTRANCE_END + FADE_DURATION;

      function frame(now: number) {
        if (cancelled) return;
        if (w === 0) {
          rafRef.current = requestAnimationFrame(frame);
          return;
        }

        const elapsed = now - t0;
        ctx!.clearRect(0, 0, w, h);

        const xOff = PAD * w;
        const yOff = PAD * h;
        const xScale = w * (1 - 2 * PAD);
        const yScale = h * (1 - 2 * PAD);

        if (elapsed < ENTRANCE_END) {
          ctx!.globalAlpha = 1;

          for (let i = 0; i < n; i++) {
            const localT = (elapsed - delays[i]) / PARTICLE_DURATION;
            if (localT < 0) continue;

            if (!arrived[i] && localT >= 0.8) {
              arrived[i] = 1;
              const c = data.charge[i];
              liveXHist[c][xBin[i]]++;
              liveYHist[c][yBin[i]]++;
            }

            const progress = localT >= 1 ? 1 : easeOutCubic(localT);

            const nx = startXNorm[i] + (endXNorm[i] - startXNorm[i]) * progress;
            const ny = startYNorm[i] + (endYNorm[i] - startYNorm[i]) * progress;

            const fadeIn = Math.min(localT * 3, 1);
            const flyAlpha = DOT_ALPHA * fadeIn;
            const alpha = SETTLED_OPACITY + (flyAlpha - SETTLED_OPACITY) * (1 - progress);
            const [r, g, b] = CHARGE_COLORS[data.charge[i]] ?? CHARGE_COLORS[2];
            ctx!.fillStyle = `rgba(${r},${g},${b},${alpha})`;
            ctx!.fillRect(
              xOff + nx * xScale,
              yOff + ny * yScale,
              DOT_SIZE, DOT_SIZE,
            );
          }

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

        } else if (elapsed < FADE_END) {
          const fadeFrac = easeOutCubic((elapsed - ENTRANCE_END) / FADE_DURATION);
          const dotAlpha = SETTLED_OPACITY + (FINAL_DOT_ALPHA - SETTLED_OPACITY) * fadeFrac;
          const overlayAlpha = 1 + (SETTLED_OPACITY - 1) * fadeFrac;

          ctx!.globalAlpha = 1;
          for (let i = 0; i < n; i++) {
            const [r, g, b] = CHARGE_COLORS[data.charge[i]] ?? CHARGE_COLORS[2];
            ctx!.fillStyle = `rgba(${r},${g},${b},${dotAlpha})`;
            ctx!.fillRect(
              xOff + endXNorm[i] * xScale,
              yOff + endYNorm[i] * yScale,
              DOT_SIZE, DOT_SIZE,
            );
          }

          ctx!.globalAlpha = overlayAlpha;
          drawMarginals(finalXHist, finalYHist);
          drawAxes(ctx!, w, h, projResult);

        } else {
          animDone = true;
          paintSettled();
          return;
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
