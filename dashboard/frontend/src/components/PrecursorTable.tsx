import { useMemo } from 'react';
import {
  createColumnHelper,
  flexRender,
  getCoreRowModel,
  useReactTable,
} from '@tanstack/react-table';
import type { PrecursorSummary } from '../api';

const columnHelper = createColumnHelper<PrecursorSummary>();

// Columns that support server-side sorting (must match backend sort_by pattern)
const SORTABLE_COLUMNS: Record<string, string> = {
  precursor_id: 'precursor_id',
  mz: 'mz',
  rt_seconds: 'rt_seconds',
  n_engines: 'n_engines',
  raw_intensity_meta: 'raw_intensity_meta',
  mobility: 'mobility',
  fragpipe_probability: 'fragpipe_probability',
  fragpipe_hyperscore: 'fragpipe_hyperscore',
  sage_hyperscore: 'sage_hyperscore',
  sage_qvalue: 'sage_qvalue',
  diann_qvalue: 'diann_qvalue',
  diann_match_tier: 'diann_match_tier',
  ms1_rt_r2: 'ms1_rt_r2',
  ms1_im_r2: 'ms1_im_r2',
  isotope_cosim: 'isotope_cosim',
  sage_cosine: 'sage_cosine',
};

interface PrecursorTableProps {
  data: PrecursorSummary[];
  selectedId: number | null;
  onSelect: (id: number, rawFile?: string) => void;
  isLoading?: boolean;
  sortBy?: string;
  sortDesc?: boolean;
  onSort?: (column: string) => void;
}

function SortIndicator({ column, currentSort, sortDesc }: { column: string; currentSort?: string; sortDesc?: boolean }) {
  const isSortable = column in SORTABLE_COLUMNS;
  const isActive = currentSort === SORTABLE_COLUMNS[column];

  if (!isSortable) return null;

  return (
    <span className={`ml-1 inline-block ${isActive ? 'text-cyan-300' : 'text-slate-500'}`}>
      {isActive ? (sortDesc ? '▼' : '▲') : '⇅'}
    </span>
  );
}

export default function PrecursorTable({
  data,
  selectedId,
  onSelect,
  isLoading,
  sortBy,
  sortDesc,
  onSort,
}: PrecursorTableProps) {
  const columns = useMemo(
    () => [
      columnHelper.accessor('precursor_id', {
        header: 'ID',
        cell: (info) => info.getValue(),
        size: 70,
      }),
      columnHelper.accessor('raw_file', {
        header: 'Raw File',
        cell: (info) => {
          const val = info.getValue();
          // Show last portion for readability
          const short = val && val.length > 20 ? '...' + val.slice(-20) : val;
          return <span className="text-xs truncate max-w-28 block text-slate-300" title={val}>{short || '-'}</span>;
        },
        size: 130,
      }),
      columnHelper.accessor('mz', {
        header: 'm/z',
        cell: (info) => info.getValue().toFixed(4),
        size: 90,
      }),
      columnHelper.accessor('charge', {
        header: 'z',
        cell: (info) => `${info.getValue()}+`,
        size: 40,
      }),
      columnHelper.accessor('rt_seconds', {
        header: 'RT',
        cell: (info) => (info.getValue() / 60).toFixed(1),
        size: 50,
      }),
      columnHelper.accessor('mobility', {
        header: '1/K0',
        cell: (info) => info.getValue().toFixed(3),
        size: 60,
      }),
      columnHelper.accessor('ccs', {
        header: 'CCS',
        cell: (info) => {
          const val = info.getValue();
          return val !== null ? val.toFixed(1) : '-';
        },
        size: 55,
      }),
      columnHelper.accessor('n_engines', {
        header: '#',
        cell: (info) => {
          const n = info.getValue();
          const colors = ['bg-slate-600', 'bg-amber-600', 'bg-sky-700', 'bg-emerald-700'];
          return (
            <span className={`px-1.5 py-0.5 rounded text-xs font-semibold text-slate-50 ${colors[n]}`}>
              {n}
            </span>
          );
        },
        size: 40,
      }),
      // FragPipe columns
      columnHelper.accessor('fragpipe_modified', {
        header: 'FragPipe',
        cell: (info) => {
          const modified = info.getValue();
          const plain = info.row.original.fragpipe_peptide;
          return (
            <span className="mono text-xs truncate max-w-32 block text-cyan-100" title={modified || plain || ''}>
              {modified || plain || '-'}
            </span>
          );
        },
        size: 130,
      }),
      columnHelper.accessor('fragpipe_probability', {
        header: 'FP Prob',
        cell: (info) => {
          const val = info.getValue();
          return val !== null ? val.toFixed(3) : '-';
        },
        size: 65,
      }),
      columnHelper.accessor('fragpipe_pep', {
        header: 'FP PEP',
        cell: (info) => {
          const val = info.getValue();
          return val !== null ? val.toExponential(1) : '-';
        },
        size: 60,
      }),
      columnHelper.accessor('fragpipe_hyperscore', {
        header: 'FP Score',
        cell: (info) => {
          const val = info.getValue();
          return val !== null ? val.toFixed(1) : '-';
        },
        size: 65,
      }),
      columnHelper.accessor('fragpipe_rt', {
        header: 'FP RT',
        cell: (info) => {
          const val = info.getValue();
          return val !== null ? (val / 60).toFixed(1) : '-';
        },
        size: 55,
      }),
      // Sage columns
      columnHelper.accessor('sage_modified', {
        header: 'Sage',
        cell: (info) => {
          const modified = info.getValue();
          const plain = info.row.original.sage_peptide;
          const tier = info.row.original.sage_match_tier;
          const isCoord = tier?.includes('COORDINATE');
          const seq = modified || plain || '-';
          return (
            <span
              className={`mono text-xs truncate max-w-32 block ${isCoord ? 'text-orange-300/65 line-through decoration-orange-300/40' : 'text-emerald-100'}`}
              title={isCoord ? `${seq} (coordinate match - unreliable)` : seq}
            >
              {seq}
            </span>
          );
        },
        size: 130,
      }),
      columnHelper.accessor('sage_qvalue', {
        header: 'S qval',
        cell: (info) => {
          const val = info.getValue();
          return val !== null ? val.toExponential(1) : '-';
        },
        size: 60,
      }),
      columnHelper.accessor('sage_pep', {
        header: 'S PEP',
        cell: (info) => {
          const val = info.getValue();
          return val !== null ? val.toExponential(1) : '-';
        },
        size: 60,
      }),
      columnHelper.accessor('sage_hyperscore', {
        header: 'S Score',
        cell: (info) => {
          const val = info.getValue();
          return val !== null ? val.toFixed(1) : '-';
        },
        size: 65,
      }),
      columnHelper.accessor('sage_rt', {
        header: 'S RT',
        cell: (info) => {
          const val = info.getValue();
          return val !== null ? (val / 60).toFixed(1) : '-';
        },
        size: 55,
      }),
      columnHelper.accessor('sage_match_tier', {
        header: 'S Tier',
        cell: (info) => {
          const tier = info.getValue();
          if (!tier) return '-';
          const short = tier.replace('SEQUENCE_', 'SEQ_').replace('COORDINATE_', 'COORD_');
          const color = tier.includes('SEQUENCE') ? 'text-emerald-300' : tier.includes('COORDINATE') ? 'text-orange-300' : 'text-amber-300';
          return <span className={`text-xs ${color}`} title={tier}>{short}</span>;
        },
        size: 80,
      }),
      // DIA-NN columns
      columnHelper.accessor('diann_modified', {
        header: 'DIA-NN',
        cell: (info) => {
          const modified = info.getValue();
          const plain = info.row.original.diann_peptide;
          const tier = info.row.original.diann_match_tier;
          const isCoord = tier?.includes('COORDINATE');
          const seq = modified || plain || '-';
          return (
            <span
              className={`mono text-xs truncate max-w-32 block ${isCoord ? 'text-orange-300/65 line-through decoration-orange-300/40' : 'text-sky-100'}`}
              title={isCoord ? `${seq} (coordinate match - unreliable)` : seq}
            >
              {seq}
            </span>
          );
        },
        size: 130,
      }),
      columnHelper.accessor('diann_qvalue', {
        header: 'D qval',
        cell: (info) => {
          const val = info.getValue();
          return val !== null ? val.toExponential(1) : '-';
        },
        size: 60,
      }),
      columnHelper.accessor('diann_pep', {
        header: 'D PEP',
        cell: (info) => {
          const val = info.getValue();
          return val !== null ? val.toExponential(1) : '-';
        },
        size: 60,
      }),
      columnHelper.accessor('diann_rt', {
        header: 'D RT',
        cell: (info) => {
          const val = info.getValue();
          return val !== null ? (val / 60).toFixed(1) : '-';
        },
        size: 55,
      }),
      columnHelper.accessor('diann_match_tier', {
        header: 'D Tier',
        cell: (info) => {
          const tier = info.getValue();
          if (!tier) return '-';
          const short = tier.replace('SEQUENCE_', 'SEQ_').replace('COORDINATE_', 'COORD_');
          const color = tier.includes('SEQUENCE') ? 'text-emerald-300' : tier.includes('COORDINATE') ? 'text-orange-300' : 'text-amber-300';
          return <span className={`text-xs ${color}`} title={tier}>{short}</span>;
        },
        size: 80,
      }),
      // Raw data columns
      columnHelper.accessor('raw_intensity_meta', {
        header: 'Int',
        cell: (info) => {
          const val = info.getValue();
          return val ? val.toExponential(1) : '-';
        },
        size: 70,
      }),
      columnHelper.accessor('collision_energy', {
        header: 'CE',
        cell: (info) => {
          const val = info.getValue();
          return val !== null ? val.toFixed(1) : '-';
        },
        size: 50,
      }),
      // Quality metrics
      columnHelper.accessor('ms1_rt_r2', {
        header: 'RT R²',
        cell: (info) => {
          const val = info.getValue();
          if (val === null) return '-';
          const color = val >= 0.9 ? 'status-good' : val >= 0.7 ? 'status-warn' : 'status-bad';
          return <span className={color}>{val.toFixed(2)}</span>;
        },
        size: 55,
      }),
      columnHelper.accessor('ms1_im_r2', {
        header: 'IM R²',
        cell: (info) => {
          const val = info.getValue();
          if (val === null) return '-';
          const color = val >= 0.9 ? 'status-good' : val >= 0.7 ? 'status-warn' : 'status-bad';
          return <span className={color}>{val.toFixed(2)}</span>;
        },
        size: 55,
      }),
      columnHelper.accessor('isotope_cosim', {
        header: 'Iso',
        cell: (info) => {
          const val = info.getValue();
          if (val === null) return '-';
          const color = val >= 0.95 ? 'status-good' : val >= 0.85 ? 'status-warn' : 'status-bad';
          return <span className={color}>{val.toFixed(2)}</span>;
        },
        size: 50,
      }),
      columnHelper.accessor('sage_cosine', {
        header: 'S Cos',
        cell: (info) => {
          const val = info.getValue();
          if (val === null || val === undefined) return '-';
          const color = val >= 0.9 ? 'status-good' : val >= 0.7 ? 'status-warn' : 'status-bad';
          return <span className={color}>{val.toFixed(3)}</span>;
        },
        size: 55,
      }),
    ],
    []
  );

  const table = useReactTable({
    data,
    columns,
    getCoreRowModel: getCoreRowModel(),
  });

  const handleHeaderClick = (columnId: string) => {
    if (onSort && columnId in SORTABLE_COLUMNS) {
      onSort(SORTABLE_COLUMNS[columnId]);
    }
  };

  if (isLoading) {
    return (
      <div className="flex items-center justify-center h-64">
        <div className="metric-pill">Loading precursors...</div>
      </div>
    );
  }

  return (
    <div className="overflow-auto h-full">
      <table className="data-table text-sm">
        <thead className="sticky top-0">
          {table.getHeaderGroups().map((headerGroup) => (
            <tr key={headerGroup.id}>
              {headerGroup.headers.map((header) => {
                const columnId = header.column.id;
                const isSortable = columnId in SORTABLE_COLUMNS;
                return (
                  <th
                    key={header.id}
                    className={`px-3 py-2 text-left ${
                      isSortable ? 'cursor-pointer hover:bg-slate-700/50 select-none transition-colors' : ''
                    }`}
                    style={{ width: header.getSize() }}
                    onClick={() => handleHeaderClick(columnId)}
                  >
                    {flexRender(header.column.columnDef.header, header.getContext())}
                    <SortIndicator column={columnId} currentSort={sortBy} sortDesc={sortDesc} />
                  </th>
                );
              })}
            </tr>
          ))}
        </thead>
        <tbody>
          {table.getRowModel().rows.map((row) => {
            const isSelected = row.original.precursor_id === selectedId;
            return (
              <tr
                key={row.id}
                onClick={() => onSelect(row.original.precursor_id, row.original.raw_file)}
                className={`
                  cursor-pointer border-b border-[#24395633] transition-colors
                  ${isSelected ? 'bg-cyan-800/25' : 'hover:bg-slate-700/20'}
                `}
              >
                {row.getVisibleCells().map((cell) => (
                  <td key={cell.id} className="px-3 py-2 text-slate-200">
                    {flexRender(cell.column.columnDef.cell, cell.getContext())}
                  </td>
                ))}
              </tr>
            );
          })}
        </tbody>
      </table>
    </div>
  );
}
