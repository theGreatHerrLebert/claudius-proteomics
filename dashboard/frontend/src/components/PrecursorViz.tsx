import { useMemo, useRef, useEffect } from 'react';
import type { PrecursorDetail } from '../api';
import MirrorSpectrumPlot from './MirrorSpectrumPlot';

interface PrecursorVizProps {
  precursor: PrecursorDetail | null;
  isLoading?: boolean;
}

// Color scale: inferno (black -> purple -> red -> orange -> yellow)
// Using color stops from matplotlib's inferno colormap
const infernoStops: [number, number, number][] = [
  [0, 0, 4],       // 0.0 - near black
  [40, 11, 84],    // 0.15 - dark purple
  [101, 21, 110],  // 0.30 - purple
  [159, 42, 99],   // 0.45 - magenta
  [212, 72, 66],   // 0.60 - red-orange
  [245, 125, 21], // 0.75 - orange
  [250, 193, 39], // 0.90 - yellow-orange
  [252, 255, 164], // 1.0 - bright yellow
];

function intensityToColor(normalizedIntensity: number): [number, number, number, number] {
  const t = Math.max(0, Math.min(1, normalizedIntensity));
  const nStops = infernoStops.length;
  const idx = t * (nStops - 1);
  const i0 = Math.floor(idx);
  const i1 = Math.min(i0 + 1, nStops - 1);
  const frac = idx - i0;

  const c0 = infernoStops[i0];
  const c1 = infernoStops[i1];

  return [
    Math.round(c0[0] + frac * (c1[0] - c0[0])),
    Math.round(c0[1] + frac * (c1[1] - c0[1])),
    Math.round(c0[2] + frac * (c1[2] - c0[2])),
    255
  ];
}

function formatAxisValue(val: number): string {
  if (Math.abs(val) >= 1000) return val.toFixed(0);
  if (Math.abs(val) >= 1) return val.toFixed(1);
  return val.toFixed(3);
}

function SpectrumPlot({
  mz,
  intensity,
  title,
  yLabel = 'Intensity',
}: {
  mz: number[];
  intensity: number[];
  title: string;
  yLabel?: string;
}) {
  if (!mz.length) {
    return (
      <div className="h-full flex items-center justify-center text-gray-500 text-sm">
        No data
      </div>
    );
  }

  const minMz = Math.min(...mz);
  const maxMz = Math.max(...mz);
  const mzRange = maxMz - minMz || 1;
  const maxInt = Math.max(...intensity);

  const width = 600;
  const height = 180;
  const margin = { left: 50, right: 20, top: 25, bottom: 30 };
  const plotW = width - margin.left - margin.right;
  const plotH = height - margin.top - margin.bottom;

  return (
    <div className="h-full relative">
      <div className="absolute top-1 left-2 text-xs text-gray-400 z-10">{title}</div>
      <svg className="w-full h-full" viewBox={`0 0 ${width} ${height}`} preserveAspectRatio="none">
        <rect x={margin.left} y={margin.top} width={plotW} height={plotH} fill="#1f2937" />
        <line x1={margin.left} y1={margin.top} x2={margin.left} y2={margin.top + plotH} stroke="#4b5563" strokeWidth={1} />
        <text x={margin.left - 5} y={margin.top + 5} fill="#9ca3af" fontSize={9} textAnchor="end">
          {formatAxisValue(maxInt)}
        </text>
        <text x={margin.left - 5} y={margin.top + plotH} fill="#9ca3af" fontSize={9} textAnchor="end">
          0
        </text>
        <text x={15} y={margin.top + plotH / 2} fill="#9ca3af" fontSize={9} textAnchor="middle"
          transform={`rotate(-90, 15, ${margin.top + plotH / 2})`}>
          {yLabel}
        </text>
        <line x1={margin.left} y1={margin.top + plotH} x2={margin.left + plotW} y2={margin.top + plotH} stroke="#4b5563" strokeWidth={1} />
        <text x={margin.left} y={height - 5} fill="#9ca3af" fontSize={9} textAnchor="start">
          {minMz.toFixed(0)}
        </text>
        <text x={margin.left + plotW} y={height - 5} fill="#9ca3af" fontSize={9} textAnchor="end">
          {maxMz.toFixed(0)}
        </text>
        <text x={margin.left + plotW / 2} y={height - 5} fill="#9ca3af" fontSize={9} textAnchor="middle">
          m/z
        </text>
        {mz.map((m, i) => {
          const x = margin.left + ((m - minMz) / mzRange) * plotW;
          const h = (intensity[i] / maxInt) * plotH;
          return (
            <line key={i} x1={x} y1={margin.top + plotH} x2={x} y2={margin.top + plotH - h}
              stroke="#3b82f6" strokeWidth={1} />
          );
        })}
      </svg>
    </div>
  );
}

function ProfilePlot({
  x, y, title, xLabel, yLabel, color,
}: {
  x: number[];
  y: number[];
  title: string;
  xLabel: string;
  yLabel: string;
  color: string;
}) {
  if (!x.length) {
    return (
      <div className="h-full flex items-center justify-center text-gray-500 text-sm">
        No data
      </div>
    );
  }

  const maxY = Math.max(...y);
  const minX = Math.min(...x);
  const maxX = Math.max(...x);
  const xRange = maxX - minX || 1;

  const width = 250;
  const height = 180;
  const margin = { left: 40, right: 10, top: 20, bottom: 25 };
  const plotW = width - margin.left - margin.right;
  const plotH = height - margin.top - margin.bottom;

  const points = x.map((xv, i) => {
    const px = margin.left + ((xv - minX) / xRange) * plotW;
    const py = margin.top + plotH - (y[i] / maxY) * plotH;
    return `${px},${py}`;
  });
  const pathD = `M ${points.join(' L ')}`;
  const areaD = `${pathD} L ${margin.left + plotW},${margin.top + plotH} L ${margin.left},${margin.top + plotH} Z`;

  return (
    <div className="h-full relative">
      <svg className="w-full h-full" viewBox={`0 0 ${width} ${height}`} preserveAspectRatio="xMidYMid meet">
        <text x={margin.left + plotW / 2} y={12} fill="#9ca3af" fontSize={10} textAnchor="middle">
          {title}
        </text>
        <rect x={margin.left} y={margin.top} width={plotW} height={plotH} fill="#1f2937" />
        <line x1={margin.left} y1={margin.top} x2={margin.left} y2={margin.top + plotH} stroke="#4b5563" strokeWidth={1} />
        <text x={margin.left - 3} y={margin.top + 8} fill="#6b7280" fontSize={7} textAnchor="end">
          {formatAxisValue(maxY)}
        </text>
        <text x={10} y={margin.top + plotH / 2} fill="#6b7280" fontSize={7} textAnchor="middle"
          transform={`rotate(-90, 10, ${margin.top + plotH / 2})`}>
          {yLabel}
        </text>
        <line x1={margin.left} y1={margin.top + plotH} x2={margin.left + plotW} y2={margin.top + plotH} stroke="#4b5563" strokeWidth={1} />
        <text x={margin.left} y={height - 5} fill="#6b7280" fontSize={7} textAnchor="start">
          {formatAxisValue(minX)}
        </text>
        <text x={margin.left + plotW} y={height - 5} fill="#6b7280" fontSize={7} textAnchor="end">
          {formatAxisValue(maxX)}
        </text>
        <text x={margin.left + plotW / 2} y={height - 5} fill="#6b7280" fontSize={7} textAnchor="middle">
          {xLabel}
        </text>
        <path d={areaD} fill={color} fillOpacity={0.3} />
        <path d={pathD} fill="none" stroke={color} strokeWidth={1.5} />
      </svg>
    </div>
  );
}

// Canvas-based heatmap component for 2D binned data (no gaps between cells)
function HeatmapPlot({
  xData,
  yData,
  intensity,
  title,
  xLabel,
  yLabel,
  xFormat = (v: number) => v.toFixed(0),
  yFormat = (v: number) => v.toFixed(3),
  useDiscreteX = false,
  useDiscreteY = false,
}: {
  xData: number[];
  yData: number[];
  intensity: number[];
  title: string;
  xLabel: string;
  yLabel: string;
  xFormat?: (v: number) => string;
  yFormat?: (v: number) => string;
  useDiscreteX?: boolean;
  useDiscreteY?: boolean;
}) {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  const { grid, bounds, maxInt, nBinsX, nBinsY } = useMemo(() => {
    if (!xData.length || !yData.length) {
      return { grid: null, bounds: null, maxInt: 0, nBinsX: 0, nBinsY: 0 };
    }

    const minX = Math.min(...xData);
    const maxX = Math.max(...xData);
    const minY = Math.min(...yData);
    const maxY = Math.max(...yData);

    // Y binning
    let nBinsY: number;
    let yBinFn: (y: number) => number;

    if (useDiscreteY) {
      // Use actual discrete values as Y bins (e.g., scan indices or mobility values)
      // Round to reduce floating point noise
      const roundedY = yData.map(y => Math.round(y * 10000) / 10000);
      const uniqueY = [...new Set(roundedY)].sort((a, b) => a - b);
      const yToIdx = new Map(uniqueY.map((y, i) => [y, i]));
      nBinsY = uniqueY.length;
      yBinFn = (y: number) => yToIdx.get(Math.round(y * 10000) / 10000) ?? 0;
    } else {
      // Bin Y values into ~100 bins for finer resolution
      nBinsY = Math.min(100, Math.max(20, Math.floor(Math.sqrt(yData.length) * 1.5)));
      const yRange = maxY - minY || 1;
      const yBinWidth = yRange / nBinsY;
      yBinFn = (y: number) => Math.min(nBinsY - 1, Math.floor((y - minY) / yBinWidth));
    }

    // X binning
    let nBinsX: number;
    let xBinFn: (x: number) => number;

    if (useDiscreteX) {
      // Use actual discrete values as X bins (e.g., unique RT frames)
      // Round to reduce floating point noise
      const roundedX = xData.map(x => Math.round(x * 1000) / 1000);
      const uniqueX = [...new Set(roundedX)].sort((a, b) => a - b);
      const xToIdx = new Map(uniqueX.map((x, i) => [x, i]));
      nBinsX = uniqueX.length;
      xBinFn = (x: number) => xToIdx.get(Math.round(x * 1000) / 1000) ?? 0;
    } else {
      // Bin X values into ~120 bins for finer resolution
      nBinsX = Math.min(120, Math.max(20, Math.floor(Math.sqrt(xData.length) * 1.5)));
      const xRange = maxX - minX || 1;
      const xBinWidth = xRange / nBinsX;
      xBinFn = (x: number) => Math.min(nBinsX - 1, Math.floor((x - minX) / xBinWidth));
    }

    // Create grid and accumulate intensities
    const grid: number[][] = Array(nBinsY).fill(null).map(() => Array(nBinsX).fill(0));

    for (let i = 0; i < xData.length; i++) {
      const xBin = xBinFn(xData[i]);
      const yBin = yBinFn(yData[i]);
      grid[yBin][xBin] += intensity[i];
    }

    const maxInt = Math.max(...grid.flat());

    return {
      grid,
      bounds: { minX, maxX, minY, maxY },
      maxInt,
      nBinsX,
      nBinsY,
    };
  }, [xData, yData, intensity, useDiscreteX, useDiscreteY]);

  // Render heatmap to canvas
  useEffect(() => {
    if (!grid || !canvasRef.current) return;

    const canvas = canvasRef.current;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    // Set canvas size to match grid dimensions (1:1 pixel mapping)
    canvas.width = nBinsX;
    canvas.height = nBinsY;

    // Create ImageData for direct pixel manipulation
    const imageData = ctx.createImageData(nBinsX, nBinsY);
    const data = imageData.data;

    const logMax = Math.log10(maxInt + 1);

    for (let yi = 0; yi < nBinsY; yi++) {
      for (let xi = 0; xi < nBinsX; xi++) {
        // Flip Y axis so higher values are at top (canvas row 0)
        const val = grid[nBinsY - 1 - yi][xi];
        const idx = (yi * nBinsX + xi) * 4;

        if (val === 0) {
          // Background color (dark gray)
          data[idx] = 31;     // R
          data[idx + 1] = 41; // G
          data[idx + 2] = 55; // B
          data[idx + 3] = 255; // A
        } else {
          const norm = Math.log10(val + 1) / logMax;
          const [r, g, b, a] = intensityToColor(norm);
          data[idx] = r;
          data[idx + 1] = g;
          data[idx + 2] = b;
          data[idx + 3] = a;
        }
      }
    }

    ctx.putImageData(imageData, 0, 0);
  }, [grid, nBinsX, nBinsY, maxInt]);

  if (!grid || !bounds) {
    return (
      <div className="h-full flex items-center justify-center text-gray-500 text-sm">
        No data
      </div>
    );
  }

  return (
    <div className="h-full relative flex flex-col">
      {/* Title bar */}
      <div className="flex-none flex justify-between px-2 py-1">
        <span className="text-xs text-gray-400">{title}</span>
        <span className="text-xs text-gray-500">{xData.length.toLocaleString()} pts</span>
      </div>

      {/* Main plot area */}
      <div className="flex-1 flex min-h-0">
        {/* Y axis label (top = max, bottom = min) */}
        <div className="flex-none w-12 flex flex-col justify-between items-end pr-1 py-1">
          <span className="text-[10px] text-gray-400">{yFormat(bounds.maxY)}</span>
          <span
            className="text-[10px] text-gray-400 origin-center"
            style={{ writingMode: 'vertical-rl', transform: 'rotate(180deg)' }}
          >
            {yLabel}
          </span>
          <span className="text-[10px] text-gray-400">{yFormat(bounds.minY)}</span>
        </div>

        {/* Canvas container */}
        <div className="flex-1 min-w-0">
          <canvas
            ref={canvasRef}
            className="w-full h-full"
            style={{ imageRendering: 'pixelated' }}
          />
        </div>
      </div>

      {/* X axis */}
      <div className="flex-none h-5 flex pl-12 pr-2">
        <div className="flex-1 flex justify-between items-start">
          <span className="text-[10px] text-gray-400">{xFormat(bounds.minX)}</span>
          <span className="text-[10px] text-gray-400">{xLabel}</span>
          <span className="text-[10px] text-gray-400">{xFormat(bounds.maxX)}</span>
        </div>
      </div>
    </div>
  );
}

export default function PrecursorViz({ precursor, isLoading }: PrecursorVizProps) {
  if (isLoading) {
    return (
      <div className="h-full flex items-center justify-center text-gray-400">
        Loading...
      </div>
    );
  }

  if (!precursor) {
    return (
      <div className="h-full flex items-center justify-center text-gray-500">
        Select a precursor to view details
      </div>
    );
  }

  // Check if we have scan data for heatmap
  const hasFragmentData = precursor.fragment_mobility && precursor.fragment_mobility.length > 0;

  return (
    <div className="h-full flex flex-col gap-2 p-2">
      {/* Header */}
      <div className="flex-none bg-gray-800 rounded p-3">
        <div className="flex items-center gap-4">
          <span className="text-lg font-semibold text-white">
            Precursor {precursor.precursor_id}
          </span>
          <span className="text-gray-400">
            m/z: <span className="text-white">{precursor.mz.toFixed(4)}</span>
          </span>
          <span className="text-gray-400">
            z: <span className="text-white">{precursor.charge}+</span>
          </span>
          <span className="text-gray-400">
            RT: <span className="text-white">{(precursor.rt_seconds / 60).toFixed(2)} min</span>
          </span>
          <span className="text-gray-400">
            1/K0: <span className="text-white">{precursor.mobility.toFixed(3)}</span>
          </span>
        </div>
        <div className="mt-2 flex gap-4 text-sm">
          <div>
            <span className="text-gray-500">FragPipe: </span>
            <span className="font-mono text-blue-400">{precursor.fragpipe_peptide || '-'}</span>
          </div>
          <div>
            <span className="text-gray-500">Sage: </span>
            <span className="font-mono text-green-400">{precursor.sage_modified || precursor.sage_peptide || '-'}</span>
          </div>
          <div>
            <span className="text-gray-500">DIA-NN: </span>
            <span className="font-mono text-purple-400">{precursor.diann_peptide || '-'}</span>
          </div>
        </div>
      </div>

      {/* Fragment spectrum - full width */}
      {/* Show mirror plot if Sage fragments available, otherwise simple spectrum */}
      <div className={`flex-none bg-gray-800 rounded overflow-hidden ${precursor.sage_matched_fragments?.length ? 'h-56' : 'h-44'}`}>
        {precursor.sage_matched_fragments?.length ? (
          <MirrorSpectrumPlot
            mz={precursor.fragment_mz}
            intensity={precursor.fragment_intensity}
            sageFragments={precursor.sage_matched_fragments}
            peptide={precursor.sage_modified || precursor.sage_peptide}
          />
        ) : (
          <SpectrumPlot
            mz={precursor.fragment_mz}
            intensity={precursor.fragment_intensity}
            title={`Fragment Spectrum (${precursor.fragment_mz.length} peaks)`}
          />
        )}
      </div>

      {/* Second row: heatmaps */}
      <div className="flex-1 grid grid-cols-2 gap-2 min-h-0">
        {/* Fragment m/z vs 1/K0 heatmap */}
        <div className="bg-gray-800 rounded overflow-hidden">
          {hasFragmentData ? (
            <HeatmapPlot
              xData={precursor.fragment_mz}
              yData={precursor.fragment_mobility}
              intensity={precursor.fragment_intensity}
              title="Fragments: 1/K0 vs m/z"
              xLabel="m/z"
              yLabel="1/K0"
              xFormat={(v) => v.toFixed(0)}
              yFormat={(v) => v.toFixed(3)}
              useDiscreteY={true}
            />
          ) : (
            <div className="h-full flex items-center justify-center text-gray-500 text-sm">
              No fragment data
            </div>
          )}
        </div>

        {/* Raw 4D: 1/K0 vs RT heatmap */}
        <div className="bg-gray-800 rounded overflow-hidden">
          {precursor.raw_rt.length > 0 ? (
            <HeatmapPlot
              xData={precursor.raw_rt.map(rt => rt / 60)}
              yData={precursor.raw_mobility}
              intensity={precursor.raw_intensity}
              title="Raw 4D: 1/K0 vs RT"
              xLabel="RT (min)"
              yLabel="1/K0"
              xFormat={(v) => v.toFixed(2)}
              yFormat={(v) => v.toFixed(3)}
              useDiscreteX={true}
              useDiscreteY={false}
            />
          ) : (
            <div className="h-full flex items-center justify-center text-gray-500 text-sm">
              No MS1 data
            </div>
          )}
        </div>
      </div>

      {/* Third row: profiles (taller) */}
      <div className="flex-none h-56 grid grid-cols-3 gap-2">
        <div className="bg-gray-800 rounded overflow-hidden">
          <ProfilePlot
            x={precursor.xic_rt.map((rt) => rt / 60)}
            y={precursor.xic_intensity}
            title="XIC"
            xLabel="RT (min)"
            yLabel="Int"
            color="#22c55e"
          />
        </div>
        <div className="bg-gray-800 rounded overflow-hidden">
          <ProfilePlot
            x={precursor.mobilogram_im}
            y={precursor.mobilogram_intensity}
            title="Mobilogram"
            xLabel="1/K0"
            yLabel="Int"
            color="#a855f7"
          />
        </div>
        <div className="bg-gray-800 rounded overflow-hidden">
          <SpectrumPlot
            mz={precursor.isotope_mz}
            intensity={precursor.isotope_intensity}
            title="Isotopes"
          />
        </div>
      </div>
    </div>
  );
}
