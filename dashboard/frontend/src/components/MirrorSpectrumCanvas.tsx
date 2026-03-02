import { useRef, useEffect } from 'react';
import {
  easeOutCubic, bounds,
  PAD, AXIS_COLOR, TICK_COLOR, LABEL_COLOR, AXIS_FONT, LABEL_FONT,
  niceTicks, formatTick,
} from './heroCanvasUtils';

/* ── types ── */

interface MatchedFragment {
  fragment_type: string;  // "b" or "y"
  ion_number: number;
  charge: number;
  mz_calculated: number;
  mz_experimental: number;
  intensity: number;
}

interface MirrorSpectrumData {
  peptide: string;
  charge: number;
  cosine: number;
  experimental_mz: number[];
  experimental_intensity: number[];
  matched_fragments: MatchedFragment[];
}

/* ── constants ── */

const B_COLOR: [number, number, number] = [95, 179, 255];
const Y_COLOR: [number, number, number] = [247, 143, 100];
const EXP_COLOR: [number, number, number] = [137, 160, 192];

const EXP_BUILD_START = 500;
const EXP_BUILD_DURATION = 3_000;
const MATCH_BUILD_START = 3_000;
const MATCH_BUILD_DURATION = 3_500;
const LABEL_DELAY = 200;        // ms after peak appears before label fades in
const LABEL_FADE_DURATION = 400;
const SEQ_BUILD_START = 5_500;
const SEQ_BUILD_DURATION = 2_000;
const COSINE_REVEAL_START = 7_500;
const COSINE_REVEAL_DURATION = 800;
const SETTLE_START = 8_500;
const FADE_DURATION = 3_000;
const SETTLED_OPACITY = 0.20;

const PEAK_LINE_WIDTH = 1.2;

/* ── component ── */

export default function MirrorSpectrumCanvas({ className = '' }: { className?: string }) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const rafRef = useRef(0);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    let cancelled = false;

    fetch('/data/mirror_spectrum.json')
      .then((r) => r.json())
      .then((json: MirrorSpectrumData) => {
        if (cancelled) return;
        startRendering(json);
      })
      .catch(() => {});

    function startRendering(data: MirrorSpectrumData) {
      const ctx = canvas!.getContext('2d');
      if (!ctx) return;

      const prefersReduced = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

      const expMz = data.experimental_mz;
      const expInt = data.experimental_intensity;
      const nExp = expMz.length;
      const matched = data.matched_fragments;

      // Sort matched by m/z for left-to-right reveal
      const sortedMatched = [...matched].sort((a, b) => a.mz_experimental - b.mz_experimental);

      // Normalize intensities to [0, 1]
      const maxExpInt = Math.max(...expInt, 1);
      const maxMatchInt = Math.max(...matched.map(m => m.intensity), 1);

      // m/z range
      const allMz = [...expMz, ...matched.map(m => m.mz_experimental)];
      const [mzLo, mzHi] = bounds(allMz);
      const mzRange = mzHi - mzLo || 1;

      // Strip mods from peptide for display
      const plainPeptide = data.peptide.replace(/\[[^\]]+\]/g, '').replace(/^-|-$/g, '');

      // Compute b/y coverage for sequence bar
      const seqLen = plainPeptide.length;
      const bCovered = new Set(matched.filter(m => m.fragment_type === 'b').map(m => m.ion_number));
      const yCovered = new Set(matched.filter(m => m.fragment_type === 'y').map(m => m.ion_number));

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

      const SEQ_BAR_H = 28;  // px for sequence bar at top

      function mzToX(mz: number): number {
        return PAD * w + ((mz - mzLo) / mzRange) * w * (1 - 2 * PAD);
      }

      function drawExpPeaks(frac: number, alpha: number) {
        const midY = (SEQ_BAR_H + h) / 2;
        const topY = SEQ_BAR_H + PAD * h;
        const peakH = midY - topY;

        ctx!.lineWidth = PEAK_LINE_WIDTH;

        for (let i = 0; i < nExp; i++) {
          const peakFrac = (frac - i / nExp * 0.7) / 0.3;
          if (peakFrac <= 0) continue;
          const grow = Math.min(easeOutCubic(peakFrac), 1);

          const x = mzToX(expMz[i]);
          const intNorm = expInt[i] / maxExpInt;
          const peakLen = intNorm * peakH * grow;

          const [r, g, b] = EXP_COLOR;
          ctx!.strokeStyle = `rgba(${r},${g},${b},${alpha * 0.6})`;
          ctx!.beginPath();
          ctx!.moveTo(x, midY);
          ctx!.lineTo(x, midY - peakLen);
          ctx!.stroke();
        }
      }

      function drawMatchedPeaks(frac: number, alpha: number, showLabels: boolean, labelAlpha: number) {
        const midY = (SEQ_BAR_H + h) / 2;
        const botY = h - PAD * h;
        const peakH = botY - midY;

        ctx!.lineWidth = PEAK_LINE_WIDTH;

        for (let i = 0; i < sortedMatched.length; i++) {
          const m = sortedMatched[i];
          const peakFrac = (frac - i / sortedMatched.length * 0.7) / 0.3;
          if (peakFrac <= 0) continue;
          const grow = Math.min(easeOutCubic(peakFrac), 1);

          const x = mzToX(m.mz_experimental);
          const intNorm = m.intensity / maxMatchInt;
          const peakLen = intNorm * peakH * 0.85 * grow;

          const color = m.fragment_type === 'b' ? B_COLOR : Y_COLOR;
          const [r, g, b] = color;
          ctx!.strokeStyle = `rgba(${r},${g},${b},${alpha * 0.8})`;
          ctx!.beginPath();
          ctx!.moveTo(x, midY);
          ctx!.lineTo(x, midY + peakLen);
          ctx!.stroke();

          // Ion label
          if (showLabels && grow >= 1) {
            const lAlpha = Math.min(labelAlpha, alpha);
            ctx!.font = '9px "IBM Plex Mono", monospace';
            ctx!.fillStyle = `rgba(${r},${g},${b},${lAlpha * 0.7})`;
            ctx!.textAlign = 'center';
            ctx!.textBaseline = 'top';
            ctx!.fillText(
              `${m.fragment_type}${m.ion_number}`,
              x,
              midY + peakLen + 3,
            );
          }
        }
      }

      function drawSequenceBar(frac: number, alpha: number) {
        if (seqLen === 0) return;

        const xOff = PAD * w + w * 0.05;
        const xEnd = w * (1 - PAD) - w * 0.05;
        const charW = Math.min((xEnd - xOff) / seqLen, 14);
        const startX = (w - charW * seqLen) / 2;

        ctx!.font = '11px "IBM Plex Mono", monospace';
        ctx!.textAlign = 'center';
        ctx!.textBaseline = 'middle';

        const revealedCount = Math.floor(frac * seqLen);

        for (let i = 0; i < seqLen; i++) {
          if (i >= revealedCount) break;

          const x = startX + (i + 0.5) * charW;
          const y = SEQ_BAR_H / 2;

          // Color by coverage
          const bHit = bCovered.has(i + 1); // b-ion number is 1-based from N-term
          const yHit = yCovered.has(seqLen - i); // y-ion number from C-term

          if (bHit && yHit) {
            ctx!.fillStyle = `rgba(16, 185, 129, ${alpha * 0.8})`; // green = both
          } else if (bHit) {
            ctx!.fillStyle = `rgba(${B_COLOR[0]},${B_COLOR[1]},${B_COLOR[2]},${alpha * 0.7})`;
          } else if (yHit) {
            ctx!.fillStyle = `rgba(${Y_COLOR[0]},${Y_COLOR[1]},${Y_COLOR[2]},${alpha * 0.7})`;
          } else {
            ctx!.fillStyle = `rgba(175, 192, 219, ${alpha * 0.35})`;
          }

          ctx!.fillText(plainPeptide[i], x, y);
        }
      }

      function drawCosineScore(alpha: number) {
        ctx!.font = '12px "IBM Plex Mono", monospace';
        ctx!.fillStyle = `rgba(82, 186, 216, ${alpha * 0.6})`;
        ctx!.textAlign = 'right';
        ctx!.textBaseline = 'top';
        ctx!.fillText(
          `cos = ${data.cosine.toFixed(2)}`,
          w - PAD * w,
          SEQ_BAR_H + PAD * h + 4,
        );
      }

      function drawMirrorAxes(alpha: number) {
        const midY = (SEQ_BAR_H + h) / 2;
        const xOff = PAD * w;
        const xEnd = w - PAD * w;

        ctx!.save();
        ctx!.globalAlpha = alpha;

        // Center line
        ctx!.strokeStyle = AXIS_COLOR;
        ctx!.lineWidth = 0.8;
        ctx!.beginPath();
        ctx!.moveTo(xOff, midY);
        ctx!.lineTo(xEnd, midY);
        ctx!.stroke();

        // m/z axis ticks along center line
        const ticks = niceTicks(mzLo, mzHi, 5);
        ctx!.font = AXIS_FONT;
        ctx!.fillStyle = TICK_COLOR;
        ctx!.textAlign = 'center';
        ctx!.textBaseline = 'bottom';

        for (const v of ticks) {
          const x = mzToX(v);
          ctx!.beginPath();
          ctx!.moveTo(x, midY - 3);
          ctx!.lineTo(x, midY + 3);
          ctx!.stroke();
          ctx!.fillText(formatTick(v), x, midY - 5);
        }

        // m/z label
        ctx!.font = LABEL_FONT;
        ctx!.fillStyle = LABEL_COLOR;
        ctx!.textBaseline = 'top';
        ctx!.fillText('m/z', (xOff + xEnd) / 2, midY + 6);

        ctx!.restore();
      }

      function paintSettled() {
        ctx!.clearRect(0, 0, w, h);
        drawExpPeaks(1, SETTLED_OPACITY);
        drawMatchedPeaks(1, SETTLED_OPACITY, true, SETTLED_OPACITY);
        drawSequenceBar(1, SETTLED_OPACITY);
        drawCosineScore(SETTLED_OPACITY);
        drawMirrorAxes(SETTLED_OPACITY);
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
      const TOTAL_ANIM = SETTLE_START + FADE_DURATION;

      function frame(now: number) {
        if (cancelled) return;
        if (w === 0) { rafRef.current = requestAnimationFrame(frame); return; }

        const elapsed = now - t0;
        ctx!.clearRect(0, 0, w, h);

        if (elapsed < SETTLE_START) {
          // Build-up phase
          const expFrac = Math.max(0, Math.min((elapsed - EXP_BUILD_START) / EXP_BUILD_DURATION, 1));
          const matchFrac = Math.max(0, Math.min((elapsed - MATCH_BUILD_START) / MATCH_BUILD_DURATION, 1));
          const seqFrac = Math.max(0, Math.min((elapsed - SEQ_BUILD_START) / SEQ_BUILD_DURATION, 1));
          const cosineFrac = Math.max(0, Math.min((elapsed - COSINE_REVEAL_START) / COSINE_REVEAL_DURATION, 1));

          // Label alpha: starts after a peak has fully appeared + delay
          const labelStartT = MATCH_BUILD_START + MATCH_BUILD_DURATION * 0.3 + LABEL_DELAY;
          const labelAlpha = Math.max(0, Math.min((elapsed - labelStartT) / LABEL_FADE_DURATION, 1));

          drawExpPeaks(expFrac, 1);
          drawMatchedPeaks(matchFrac, 1, labelAlpha > 0, labelAlpha);
          drawSequenceBar(seqFrac, 1);
          if (cosineFrac > 0) drawCosineScore(easeOutCubic(cosineFrac));
          drawMirrorAxes(1);

        } else if (elapsed < TOTAL_ANIM) {
          // Fade to settled
          const fadeFrac = easeOutCubic((elapsed - SETTLE_START) / FADE_DURATION);
          const alpha = 1 + (SETTLED_OPACITY - 1) * fadeFrac;

          drawExpPeaks(1, alpha);
          drawMatchedPeaks(1, alpha, true, alpha);
          drawSequenceBar(1, alpha);
          drawCosineScore(alpha);
          drawMirrorAxes(alpha);

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
