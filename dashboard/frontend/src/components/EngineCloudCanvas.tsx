import { useRef, useEffect } from 'react';
import {
  bounds, easeOutCubic, hash, drawAxes, drawMarginalCurve,
  PAD,
  type ProjectionResult,
} from './heroCanvasUtils';

/* ── types ── */

interface EngineCloudData {
  mz: number[];
  rt_norm: number[];
  mobility: number[];
  n_engines: number[];
}

type Projection = 'mz_vs_im' | 'rt_vs_mz';

/* ── constants ── */

const ENGINE_COLORS: Record<number, [number, number, number]> = {
  3: [16, 185, 129],    // emerald green – 3-engine consensus
  2: [82, 186, 216],    // teal/cyan – 2-engine agreement
  1: [74, 85, 104],     // dim gray – single engine
};
const ENGINE_KEYS = [3, 2, 1]; // draw order: 1 first (behind), 3 last (on top)

const DOT_SIZE = 1.5;
const DOT_ALPHA = 0.55;

const ANIM_DURATION = 8_000;
const PARTICLE_DURATION = 1_200;
const MAX_STAGGER = 8_000;

const SETTLED_OPACITY = 0.28;
const FADE_DURATION = 6_000;

const N_BINS = 120;
const MARGINAL_X_HEIGHT = 50;
const MARGINAL_Y_WIDTH = 40;

/* ── helpers ── */

function pickProjection(): Projection {
  return Math.random() < 0.5 ? 'mz_vs_im' : 'rt_vs_mz';
}

function projectData(data: EngineCloudData, proj: Projection): ProjectionResult {
  const n = data.mz.length;
  let xRaw: number[], yRaw: number[];
  let xLabel: string, yLabel: string;

  if (proj === 'mz_vs_im') {
    xRaw = data.mz;
    yRaw = data.mobility;
    xLabel = 'm/z';
    yLabel = '1/K\u2080';
  } else {
    xRaw = data.rt_norm;
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
    endYNorm[i] = (yHi - yRaw[i]) / yRange;
  }

  return { endXNorm, endYNorm, xLabel, yLabel, xLo, xHi, yLo, yHi };
}

/* ── component ── */

export default function EngineCloudCanvas({ className = '' }: { className?: string }) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const rafRef = useRef(0);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    let cancelled = false;
    const projection = pickProjection();

    fetch('/data/engine_cloud.json')
      .then((r) => r.json())
      .then((json: EngineCloudData) => {
        if (cancelled) return;
        startRendering(json, projection);
      })
      .catch(() => {});

    function startRendering(data: EngineCloudData, proj: Projection) {
      const ctx = canvas!.getContext('2d');
      if (!ctx) return;

      const prefersReduced = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

      const n = data.mz.length;
      const projResult = projectData(data, proj);
      const { endXNorm, endYNorm } = projResult;

      const startXNorm = new Float32Array(n);
      const startYNorm = new Float32Array(n);
      const delays = new Float32Array(n);

      const xBin = new Uint8Array(n);
      const yBin = new Uint8Array(n);

      const finalXHist: Record<number, Float32Array> = {};
      const finalYHist: Record<number, Float32Array> = {};
      for (const k of ENGINE_KEYS) {
        finalXHist[k] = new Float32Array(N_BINS);
        finalYHist[k] = new Float32Array(N_BINS);
      }

      for (let i = 0; i < n; i++) {
        startXNorm[i] = -0.05 - hash(i) * 0.15;
        startYNorm[i] = endYNorm[i];
        delays[i] = endXNorm[i] * MAX_STAGGER;

        const bx = Math.min(Math.floor(endXNorm[i] * N_BINS), N_BINS - 1);
        const by = Math.min(Math.floor(endYNorm[i] * N_BINS), N_BINS - 1);
        xBin[i] = bx;
        yBin[i] = by;
        const eng = data.n_engines[i];
        if (finalXHist[eng]) {
          finalXHist[eng][bx]++;
          finalYHist[eng][by]++;
        }
      }

      // Normalize histograms
      const stackedX = new Float32Array(N_BINS);
      const stackedY = new Float32Array(N_BINS);
      for (const k of ENGINE_KEYS) {
        for (let b = 0; b < N_BINS; b++) {
          stackedX[b] += finalXHist[k][b];
          stackedY[b] += finalYHist[k][b];
        }
      }
      let maxX = 0, maxY = 0;
      for (let b = 0; b < N_BINS; b++) {
        if (stackedX[b] > maxX) maxX = stackedX[b];
        if (stackedY[b] > maxY) maxY = stackedY[b];
      }
      for (const k of ENGINE_KEYS) {
        for (let b = 0; b < N_BINS; b++) {
          finalXHist[k][b] /= maxX || 1;
          finalYHist[k][b] /= maxY || 1;
        }
      }

      // Live histograms
      const liveXHist: Record<number, Float32Array> = {};
      const liveYHist: Record<number, Float32Array> = {};
      for (const k of ENGINE_KEYS) {
        liveXHist[k] = new Float32Array(N_BINS);
        liveYHist[k] = new Float32Array(N_BINS);
      }
      const arrived = new Uint8Array(n);

      const FINAL_DOT_ALPHA = DOT_ALPHA * SETTLED_OPACITY;
      let staticCanvas: HTMLCanvasElement | null = null;

      const settledColors: string[] = new Array(n);
      for (let i = 0; i < n; i++) {
        const eng = data.n_engines[i];
        const [r, g, b] = ENGINE_COLORS[eng] ?? ENGINE_COLORS[1];
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
          octx.fillRect(xOff + endXNorm[i] * xScale, yOff + endYNorm[i] * yScale, DOT_SIZE, DOT_SIZE);
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

      function drawMarginals(
        xHist: Record<number, Float32Array>,
        yHist: Record<number, Float32Array>,
      ) {
        const xOff = PAD * w;
        const yOff = PAD * h;
        const xEnd = xOff + w * (1 - 2 * PAD);
        const yEnd = yOff + h * (1 - 2 * PAD);

        // Draw in order: 1 (back), 2, 3 (front)
        let xStack: Float32Array | null = null;
        for (const k of [1, 2, 3]) {
          drawMarginalCurve(ctx!, xHist[k], ENGINE_COLORS[k], 'x', yEnd, yEnd - MARGINAL_X_HEIGHT, xOff, xEnd, xStack);
          if (!xStack) xStack = new Float32Array(N_BINS);
          for (let b = 0; b < N_BINS; b++) xStack[b] += xHist[k][b];
        }

        let yStack: Float32Array | null = null;
        for (const k of [1, 2, 3]) {
          drawMarginalCurve(ctx!, yHist[k], ENGINE_COLORS[k], 'y', xOff, xOff + MARGINAL_Y_WIDTH, yOff, yEnd, yStack);
          if (!yStack) yStack = new Float32Array(N_BINS);
          for (let b = 0; b < N_BINS; b++) yStack[b] += yHist[k][b];
        }
      }

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
        if (w === 0) { rafRef.current = requestAnimationFrame(frame); return; }

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
              const eng = data.n_engines[i];
              if (liveXHist[eng]) {
                liveXHist[eng][xBin[i]]++;
                liveYHist[eng][yBin[i]]++;
              }
            }

            const progress = localT >= 1 ? 1 : easeOutCubic(localT);
            const nx = startXNorm[i] + (endXNorm[i] - startXNorm[i]) * progress;
            const ny = startYNorm[i] + (endYNorm[i] - startYNorm[i]) * progress;

            const fadeIn = Math.min(localT * 3, 1);
            const flyAlpha = DOT_ALPHA * fadeIn;
            const alpha = SETTLED_OPACITY + (flyAlpha - SETTLED_OPACITY) * (1 - progress);
            const eng = data.n_engines[i];
            const [r, g, b] = ENGINE_COLORS[eng] ?? ENGINE_COLORS[1];
            ctx!.fillStyle = `rgba(${r},${g},${b},${alpha})`;
            ctx!.fillRect(xOff + nx * xScale, yOff + ny * yScale, DOT_SIZE, DOT_SIZE);
          }

          const liveNormX: Record<number, Float32Array> = {};
          const liveNormY: Record<number, Float32Array> = {};
          for (const k of ENGINE_KEYS) {
            liveNormX[k] = new Float32Array(N_BINS);
            liveNormY[k] = new Float32Array(N_BINS);
            for (let b = 0; b < N_BINS; b++) {
              liveNormX[k][b] = liveXHist[k][b] / (maxX || 1);
              liveNormY[k][b] = liveYHist[k][b] / (maxY || 1);
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
            const eng = data.n_engines[i];
            const [r, g, b] = ENGINE_COLORS[eng] ?? ENGINE_COLORS[1];
            ctx!.fillStyle = `rgba(${r},${g},${b},${dotAlpha})`;
            ctx!.fillRect(xOff + endXNorm[i] * xScale, yOff + endYNorm[i] * yScale, DOT_SIZE, DOT_SIZE);
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
