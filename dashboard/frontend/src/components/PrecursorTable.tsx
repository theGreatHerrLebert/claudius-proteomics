import { useMemo } from 'react';
import {
  createColumnHelper,
  flexRender,
  getCoreRowModel,
  useReactTable,
} from '@tanstack/react-table';
import type { PrecursorSummary } from '../api';

const columnHelper = createColumnHelper<PrecursorSummary>();

interface PrecursorTableProps {
  data: PrecursorSummary[];
  selectedId: number | null;
  onSelect: (id: number) => void;
  isLoading?: boolean;
}

export default function PrecursorTable({
  data,
  selectedId,
  onSelect,
  isLoading,
}: PrecursorTableProps) {
  const columns = useMemo(
    () => [
      columnHelper.accessor('precursor_id', {
        header: 'ID',
        cell: (info) => info.getValue(),
        size: 70,
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
      columnHelper.accessor('n_engines', {
        header: '#',
        cell: (info) => {
          const n = info.getValue();
          const colors = ['bg-gray-600', 'bg-yellow-600', 'bg-blue-600', 'bg-green-600'];
          return (
            <span className={`px-1.5 py-0.5 rounded text-xs font-medium ${colors[n]}`}>
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
            <span className="font-mono text-xs truncate max-w-32 block" title={modified || plain || ''}>
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
      // Sage columns
      columnHelper.accessor('sage_modified', {
        header: 'Sage',
        cell: (info) => {
          const modified = info.getValue();
          const plain = info.row.original.sage_peptide;
          return (
            <span className="font-mono text-xs truncate max-w-32 block" title={modified || plain || ''}>
              {modified || plain || '-'}
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
      columnHelper.accessor('sage_match_tier', {
        header: 'S Tier',
        cell: (info) => {
          const tier = info.getValue();
          if (!tier) return '-';
          const short = tier.replace('SEQUENCE_', 'SEQ_').replace('COORDINATE_', 'COORD_');
          const color = tier.includes('SEQUENCE') ? 'text-green-400' : 'text-yellow-400';
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
          return (
            <span className="font-mono text-xs truncate max-w-32 block" title={modified || plain || ''}>
              {modified || plain || '-'}
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
      columnHelper.accessor('diann_match_tier', {
        header: 'D Tier',
        cell: (info) => {
          const tier = info.getValue();
          if (!tier) return '-';
          const short = tier.replace('SEQUENCE_', 'SEQ_').replace('COORDINATE_', 'COORD_');
          const color = tier.includes('SEQUENCE') ? 'text-green-400' : 'text-yellow-400';
          return <span className={`text-xs ${color}`} title={tier}>{short}</span>;
        },
        size: 80,
      }),
      columnHelper.accessor('diann_ccs', {
        header: 'CCS',
        cell: (info) => {
          const val = info.getValue();
          return val !== null ? val.toFixed(1) : '-';
        },
        size: 55,
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
      columnHelper.accessor('frame_id', {
        header: 'Frame',
        cell: (info) => info.getValue() ?? '-',
        size: 60,
      }),
      // Quality metrics
      columnHelper.accessor('ms1_rt_r2', {
        header: 'RT R²',
        cell: (info) => {
          const val = info.getValue();
          if (val === null) return '-';
          const color = val >= 0.9 ? 'text-green-400' : val >= 0.7 ? 'text-yellow-400' : 'text-red-400';
          return <span className={color}>{val.toFixed(2)}</span>;
        },
        size: 55,
      }),
      columnHelper.accessor('ms1_im_r2', {
        header: 'IM R²',
        cell: (info) => {
          const val = info.getValue();
          if (val === null) return '-';
          const color = val >= 0.9 ? 'text-green-400' : val >= 0.7 ? 'text-yellow-400' : 'text-red-400';
          return <span className={color}>{val.toFixed(2)}</span>;
        },
        size: 55,
      }),
      columnHelper.accessor('isotope_cosim', {
        header: 'Iso',
        cell: (info) => {
          const val = info.getValue();
          if (val === null) return '-';
          const color = val >= 0.95 ? 'text-green-400' : val >= 0.85 ? 'text-yellow-400' : 'text-red-400';
          return <span className={color}>{val.toFixed(2)}</span>;
        },
        size: 50,
      }),
    ],
    []
  );

  const table = useReactTable({
    data,
    columns,
    getCoreRowModel: getCoreRowModel(),
  });

  if (isLoading) {
    return (
      <div className="flex items-center justify-center h-64 text-gray-400">
        Loading precursors...
      </div>
    );
  }

  return (
    <div className="overflow-auto h-full">
      <table className="w-full text-sm">
        <thead className="sticky top-0 bg-gray-800 text-gray-300">
          {table.getHeaderGroups().map((headerGroup) => (
            <tr key={headerGroup.id}>
              {headerGroup.headers.map((header) => (
                <th
                  key={header.id}
                  className="px-3 py-2 text-left font-medium border-b border-gray-700"
                  style={{ width: header.getSize() }}
                >
                  {flexRender(header.column.columnDef.header, header.getContext())}
                </th>
              ))}
            </tr>
          ))}
        </thead>
        <tbody>
          {table.getRowModel().rows.map((row) => {
            const isSelected = row.original.precursor_id === selectedId;
            return (
              <tr
                key={row.id}
                onClick={() => onSelect(row.original.precursor_id)}
                className={`
                  cursor-pointer border-b border-gray-800 transition-colors
                  ${isSelected ? 'bg-blue-900/50' : 'hover:bg-gray-800/50'}
                `}
              >
                {row.getVisibleCells().map((cell) => (
                  <td key={cell.id} className="px-3 py-2 text-gray-300">
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
