import { useMemo } from 'react';
import type { SageMatchedFragment } from '../api';

interface MirrorSpectrumPlotProps {
  mz: number[];
  intensity: number[];
  sageFragments: SageMatchedFragment[];
  peptide: string | null;
}

function formatAxisValue(val: number): string {
  if (Math.abs(val) >= 1000) return val.toFixed(0);
  if (Math.abs(val) >= 1) return val.toFixed(1);
  return val.toFixed(3);
}

export default function MirrorSpectrumPlot({
  mz,
  intensity,
  sageFragments,
  peptide,
}: MirrorSpectrumPlotProps) {
  // Process data for plotting
  const plotData = useMemo(() => {
    if (!mz.length && !sageFragments.length) {
      return null;
    }

    // Calculate m/z range from experimental spectrum and theoretical fragments
    const allMz = [
      ...mz,
      ...sageFragments.map(f => f.mz_calculated),
    ];
    const minMz = Math.min(...allMz);
    const maxMz = Math.max(...allMz);
    const mzRange = maxMz - minMz || 1;
    const mzPadding = mzRange * 0.05;

    // Calculate intensity ranges
    const maxExpInt = mz.length > 0 ? Math.max(...intensity) : 1;
    // Use Sage intensities for the theoretical panel (normalized display)
    const maxFragInt = sageFragments.length > 0
      ? Math.max(...sageFragments.map(f => f.intensity))
      : 1;

    return {
      minMz: minMz - mzPadding,
      maxMz: maxMz + mzPadding,
      mzRange: mzRange + 2 * mzPadding,
      maxExpInt,
      maxFragInt,
    };
  }, [mz, intensity, sageFragments]);

  if (!plotData) {
    return (
      <div className="h-full flex items-center justify-center text-gray-500 text-sm">
        No data
      </div>
    );
  }

  const width = 600;
  const height = 280;
  const margin = { left: 50, right: 20, top: 25, bottom: 30 };
  const plotW = width - margin.left - margin.right;
  const plotH = height - margin.top - margin.bottom;
  const halfPlotH = plotH / 2;
  const centerY = margin.top + halfPlotH;

  // Colors
  const bIonColor = '#3b82f6';  // blue-500
  const yIonColor = '#ef4444';  // red-500
  const expColor = '#6b7280';   // gray-500

  // Scale functions
  const xScale = (m: number) => margin.left + ((m - plotData.minMz) / plotData.mzRange) * plotW;
  const yScaleExp = (i: number) => centerY - (i / plotData.maxExpInt) * halfPlotH * 0.9;
  const yScaleFrag = (i: number) => centerY + (i / plotData.maxFragInt) * halfPlotH * 0.9;

  // Build ion label
  const getIonLabel = (frag: SageMatchedFragment) => {
    const base = `${frag.fragment_type}${frag.ion_number}`;
    return frag.charge > 1 ? `${base}+${frag.charge}` : base;
  };

  // Count b and y ions
  const nBIons = sageFragments.filter(f => f.fragment_type === 'b').length;
  const nYIons = sageFragments.filter(f => f.fragment_type === 'y').length;

  return (
    <div className="h-full relative">
      {/* Title and legend */}
      <div className="absolute top-1 left-2 text-xs text-gray-400 z-10 flex items-center gap-4">
        <span>Spectrum ({mz.length} peaks)</span>
        <span className="flex items-center gap-1">
          <span className="inline-block w-3 h-0.5" style={{ backgroundColor: bIonColor }}></span>
          <span>b ({nBIons})</span>
        </span>
        <span className="flex items-center gap-1">
          <span className="inline-block w-3 h-0.5" style={{ backgroundColor: yIonColor }}></span>
          <span>y ({nYIons})</span>
        </span>
        <span className="text-yellow-500">(theoretical)</span>
      </div>
      {peptide && (
        <div className="absolute top-1 right-2 text-xs text-green-400 z-10 font-mono">
          {peptide}
        </div>
      )}

      <svg className="w-full h-full" viewBox={`0 0 ${width} ${height}`} preserveAspectRatio="none">
        {/* Background */}
        <rect x={margin.left} y={margin.top} width={plotW} height={plotH} fill="#1f2937" />

        {/* Center axis line (m/z axis) */}
        <line
          x1={margin.left}
          y1={centerY}
          x2={margin.left + plotW}
          y2={centerY}
          stroke="#4b5563"
          strokeWidth={1}
        />

        {/* Y-axis line */}
        <line
          x1={margin.left}
          y1={margin.top}
          x2={margin.left}
          y2={margin.top + plotH}
          stroke="#4b5563"
          strokeWidth={1}
        />

        {/* Y-axis labels */}
        <text x={margin.left - 5} y={margin.top + 10} fill="#9ca3af" fontSize={8} textAnchor="end">
          {formatAxisValue(plotData.maxExpInt)}
        </text>
        <text x={margin.left - 5} y={centerY - 5} fill="#9ca3af" fontSize={8} textAnchor="end">
          0
        </text>
        <text x={margin.left - 5} y={margin.top + plotH - 5} fill="#9ca3af" fontSize={8} textAnchor="end">
          {formatAxisValue(plotData.maxFragInt)}
        </text>

        {/* Y-axis title - Experimental */}
        <text
          x={15}
          y={margin.top + halfPlotH / 2}
          fill="#9ca3af"
          fontSize={8}
          textAnchor="middle"
          transform={`rotate(-90, 15, ${margin.top + halfPlotH / 2})`}
        >
          Exp
        </text>

        {/* Y-axis title - Theoretical */}
        <text
          x={15}
          y={centerY + halfPlotH / 2}
          fill="#9ca3af"
          fontSize={8}
          textAnchor="middle"
          transform={`rotate(-90, 15, ${centerY + halfPlotH / 2})`}
        >
          b/y
        </text>

        {/* X-axis labels */}
        <text x={margin.left} y={height - 5} fill="#9ca3af" fontSize={9} textAnchor="start">
          {plotData.minMz.toFixed(0)}
        </text>
        <text x={margin.left + plotW} y={height - 5} fill="#9ca3af" fontSize={9} textAnchor="end">
          {plotData.maxMz.toFixed(0)}
        </text>
        <text x={margin.left + plotW / 2} y={height - 5} fill="#9ca3af" fontSize={9} textAnchor="middle">
          m/z
        </text>

        {/* Experimental spectrum (top panel, pointing up) */}
        {mz.map((m, i) => {
          const x = xScale(m);
          const y1 = centerY;
          const y2 = yScaleExp(intensity[i]);
          return (
            <line
              key={`exp-${i}`}
              x1={x}
              y1={y1}
              x2={x}
              y2={y2}
              stroke={expColor}
              strokeWidth={1}
            />
          );
        })}

        {/* Theoretical b/y ions (bottom panel, pointing down, with labels) */}
        {sageFragments.map((frag, i) => {
          const x = xScale(frag.mz_calculated);
          const y1 = centerY;
          const y2 = yScaleFrag(frag.intensity);
          const color = frag.fragment_type === 'b' ? bIonColor : yIonColor;
          const label = getIonLabel(frag);

          return (
            <g key={`frag-${i}`}>
              <line
                x1={x}
                y1={y1}
                x2={x}
                y2={y2}
                stroke={color}
                strokeWidth={1.5}
              />
              {/* Ion label */}
              <text
                x={x}
                y={y2 + 10}
                fill={color}
                fontSize={7}
                textAnchor="middle"
              >
                {label}
              </text>
            </g>
          );
        })}
      </svg>
    </div>
  );
}
