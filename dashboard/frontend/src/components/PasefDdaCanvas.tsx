import { useRef, useEffect } from 'react';
import {
  easeOutCubic, hash, drawAxes,
  PAD,
  type ProjectionResult,
} from './heroCanvasUtils';

/* ── types ── */

interface PasefWindow {
  frame_offset: number;
  scan_begin: number;
  scan_end: number;
  mz_center: number;
  mz_width: number;
}

interface PasefCycle {
  n_frames: number;
  windows: PasefWindow[];
}

interface PasefDdaData {
  background_ions: {
    mz: number[];
    scan_norm: number[];
    intensity: number[];
    charge: number[];
  };
  cycles: PasefCycle[];
  axes: { mz_lo: number; mz_hi: number; scan_lo: number; scan_hi: number };
}

/* ── constants ── */

const CHARGE_COLORS: Record<number, [number, number, number]> = {
  1: [100, 130, 170],
  2: [82, 186, 216],
  3: [100, 160, 220],
  4: [130, 200, 230],
};

const DOT_SIZE = 1.5;
const DOT_ALPHA = 0.35;

// Timing
const BG_FADE_IN = 2_000;             // Background ions fade in
const CYCLE_START = 2_500;             // First isolation window appears
const WINDOW_INTERVAL = 120;           // ms between window appearances
const WINDOW_GLOW_DURATION = 400;      // ms window stays bright
const CYCLE_PAUSE = 600;               // ms pause between cycles
const SETTLE_START_AFTER = 2_000;      // ms after last window before settling
const FADE_DURATION = 4_000;
const SETTLED_OPACITY = 0.15;

// Window appearance
const WINDOW_BORDER_COLOR = '#52bad8';
const WINDOW_FILL = 'rgba(82, 186, 216, 0.08)';
const WINDOW_GLOW_FILL = 'rgba(82, 186, 216, 0.25)';
const WINDOW_GLOW_BORDER = '#8de4f4';

/* ── component ── */

export default function PasefDdaCanvas({ className = '' }: { className?: string }) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const rafRef = useRef(0);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    let cancelled = false;

    fetch('/data/pasef_dda.json')
      .then((r) => r.json())
      .then((json: PasefDdaData) => {
        if (cancelled) return;
        startRendering(json);
      })
      .catch(() => {});

    function startRendering(data: PasefDdaData) {
      const ctx = canvas!.getContext('2d');
      if (!ctx) return;

      const prefersReduced = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

      const bg = data.background_ions;
      const n = bg.mz.length;
      const axes = data.axes;
      const mzRange = axes.mz_hi - axes.mz_lo || 1;

      // Pre-compute normalized positions for background ions
      const bgXNorm = new Float32Array(n);
      const bgYNorm = new Float32Array(n);
      for (let i = 0; i < n; i++) {
        bgXNorm[i] = (bg.mz[i] - axes.mz_lo) / mzRange;
        bgYNorm[i] = bg.scan_norm[i]; // no invert: match scatter plot orientation
      }

      // Flatten all windows across all cycles with timing
      interface TimedWindow {
        scanBegin: number;
        scanEnd: number;
        mzNormLeft: number;
        mzNormRight: number;
        startTime: number;
      }

      const allWindows: TimedWindow[] = [];
      let t = CYCLE_START;
      for (const cycle of data.cycles) {
        for (const win of cycle.windows) {
          const mzLeft = (win.mz_center - win.mz_width / 2 - axes.mz_lo) / mzRange;
          const mzRight = (win.mz_center + win.mz_width / 2 - axes.mz_lo) / mzRange;
          allWindows.push({
            scanBegin: win.scan_begin,   // no invert: match scatter plot orientation
            scanEnd: win.scan_end,
            mzNormLeft: mzLeft,
            mzNormRight: mzRight,
            startTime: t,
          });
          t += WINDOW_INTERVAL;
        }
        t += CYCLE_PAUSE;
      }

      const lastWindowEnd = allWindows.length > 0
        ? allWindows[allWindows.length - 1].startTime + WINDOW_GLOW_DURATION
        : CYCLE_START;
      const SETTLE_TIME = lastWindowEnd + SETTLE_START_AFTER;
      const TOTAL_ANIM = SETTLE_TIME + FADE_DURATION;

      // Projection result for axes
      const projResult: ProjectionResult = {
        endXNorm: new Float32Array(0),
        endYNorm: new Float32Array(0),
        xLabel: 'm/z',
        yLabel: 'Ion Mobility',
        xLo: axes.mz_lo,
        xHi: axes.mz_hi,
        yLo: 0,
        yHi: 1,
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

      function drawBackground(alpha: number) {
        const xOff = PAD * w;
        const yOff = PAD * h;
        const xScale = w * (1 - 2 * PAD);
        const yScale = h * (1 - 2 * PAD);

        for (let i = 0; i < n; i++) {
          const c = bg.charge[i];
          const [r, g, b] = CHARGE_COLORS[c] ?? CHARGE_COLORS[2];
          const intAlpha = alpha * (0.3 + 0.7 * bg.intensity[i]);
          ctx!.fillStyle = `rgba(${r},${g},${b},${intAlpha})`;
          ctx!.fillRect(
            xOff + bgXNorm[i] * xScale,
            yOff + bgYNorm[i] * yScale,
            DOT_SIZE, DOT_SIZE,
          );
        }
      }

      function drawWindow(win: TimedWindow, age: number) {
        const xOff = PAD * w;
        const yOff = PAD * h;
        const xScale = w * (1 - 2 * PAD);
        const yScale = h * (1 - 2 * PAD);

        const x1 = xOff + win.mzNormLeft * xScale;
        const x2 = xOff + win.mzNormRight * xScale;
        const y1 = yOff + win.scanBegin * yScale;
        const y2 = yOff + win.scanEnd * yScale;
        const rw = Math.max(x2 - x1, 3);
        const rh = Math.max(y2 - y1, 3);

        // Grow-in animation
        const growT = Math.min(age / 200, 1);
        const grow = easeOutCubic(growT);
        const cx = (x1 + x2) / 2;
        const cy = (y1 + y2) / 2;
        const gw = rw * grow;
        const gh = rh * grow;

        // Glow fades after WINDOW_GLOW_DURATION
        const glowT = Math.max(0, 1 - age / WINDOW_GLOW_DURATION);
        const isGlowing = glowT > 0;

        ctx!.fillStyle = isGlowing ? WINDOW_GLOW_FILL : WINDOW_FILL;
        ctx!.fillRect(cx - gw / 2, cy - gh / 2, gw, gh);

        ctx!.strokeStyle = isGlowing ? WINDOW_GLOW_BORDER : WINDOW_BORDER_COLOR;
        ctx!.lineWidth = isGlowing ? 1.5 : 0.8;
        ctx!.globalAlpha = isGlowing ? 0.8 + 0.2 * glowT : 0.3;
        ctx!.strokeRect(cx - gw / 2, cy - gh / 2, gw, gh);
        ctx!.globalAlpha = 1;
      }

      function paintSettled() {
        ctx!.clearRect(0, 0, w, h);
        drawBackground(SETTLED_OPACITY * DOT_ALPHA);
        ctx!.globalAlpha = SETTLED_OPACITY;
        for (const win of allWindows) {
          drawWindow(win, 10_000);
        }
        ctx!.globalAlpha = SETTLED_OPACITY;
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

      function frame(now: number) {
        if (cancelled) return;
        if (w === 0) { rafRef.current = requestAnimationFrame(frame); return; }

        const elapsed = now - t0;
        ctx!.clearRect(0, 0, w, h);

        if (elapsed < SETTLE_TIME) {
          // Background fade-in
          const bgAlpha = Math.min(elapsed / BG_FADE_IN, 1) * DOT_ALPHA;
          drawBackground(bgAlpha);

          // Draw visible windows
          for (const win of allWindows) {
            const age = elapsed - win.startTime;
            if (age < 0) continue;
            drawWindow(win, age);
          }

          drawAxes(ctx!, w, h, projResult);

        } else if (elapsed < TOTAL_ANIM) {
          // Fade to settled
          const fadeFrac = easeOutCubic((elapsed - SETTLE_TIME) / FADE_DURATION);
          const alpha = DOT_ALPHA + (SETTLED_OPACITY * DOT_ALPHA - DOT_ALPHA) * fadeFrac;
          drawBackground(alpha);

          const winAlpha = 1 + (SETTLED_OPACITY - 1) * fadeFrac;
          ctx!.globalAlpha = winAlpha;
          for (const win of allWindows) {
            drawWindow(win, 10_000);
          }
          ctx!.globalAlpha = winAlpha;
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
