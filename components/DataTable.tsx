import React, { useState } from 'react';
import { ProteinData } from '../types';
import { Search } from 'lucide-react';
import { SPECIES_CONFIG } from '../constants';

interface DataTableProps {
  data: ProteinData[];
  onSelect: (protein: ProteinData) => void;
}

const DataTable: React.FC<DataTableProps> = ({ data, onSelect }) => {
  const [filter, setFilter] = useState('');
  
  const filteredData = data.filter(d => 
    d.gene.toLowerCase().includes(filter.toLowerCase()) ||
    d.description.toLowerCase().includes(filter.toLowerCase())
  ).slice(0, 100); // Limit display for performance

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden flex flex-col h-[500px]">
      <div className="p-4 border-b border-slate-200 flex justify-between items-center bg-slate-50">
        <h3 className="text-lg font-semibold text-slate-800">Top Data (List)</h3>
        <div className="relative">
          <Search className="absolute left-3 top-1/2 transform -translate-y-1/2 text-slate-400" size={16} />
          <input 
            type="text" 
            placeholder="Search gene..." 
            className="pl-9 pr-4 py-2 border border-slate-300 rounded-lg text-sm focus:outline-none focus:ring-2 focus:ring-indigo-500"
            value={filter}
            onChange={(e) => setFilter(e.target.value)}
          />
        </div>
      </div>
      <div className="overflow-auto flex-1">
        <table className="w-full text-sm text-left">
          <thead className="text-xs text-slate-500 uppercase bg-slate-50 sticky top-0 z-10">
            <tr>
              <th className="px-6 py-3">Gene</th>
              <th className="px-6 py-3">Species</th>
              <th className="px-6 py-3">Log2 FC</th>
              <th className="px-6 py-3">Intensity</th>
              <th className="px-6 py-3">Action</th>
            </tr>
          </thead>
          <tbody>
            {filteredData.map((d) => (
              <tr key={d.id} className="bg-white border-b hover:bg-slate-50 transition-colors">
                <td className="px-6 py-4 font-medium text-slate-900">{d.gene}</td>
                <td className="px-6 py-4">
                    <span 
                        className="text-xs px-2 py-1 rounded-full font-medium border"
                        style={{ 
                            backgroundColor: `${SPECIES_CONFIG[d.species].color}20`, // 20 hex = low opacity
                            borderColor: `${SPECIES_CONFIG[d.species].color}50`,
                            color: SPECIES_CONFIG[d.species].color 
                        }}
                    >
                        {d.species}
                    </span>
                </td>
                <td className={`px-6 py-4 font-mono ${
                    Math.abs(d.log2FoldChange) > 1 ? 'font-bold' : ''
                }`}>
                    {d.log2FoldChange.toFixed(2)}
                </td>
                <td className="px-6 py-4 font-mono text-slate-600">{d.averageIntensityLog10.toFixed(2)}</td>
                <td className="px-6 py-4">
                  <button 
                    onClick={() => onSelect(d)}
                    className="text-indigo-600 hover:text-indigo-900 font-medium text-xs border border-indigo-200 bg-indigo-50 px-3 py-1 rounded-full hover:bg-indigo-100"
                  >
                    Details
                  </button>
                </td>
              </tr>
            ))}
            {filteredData.length === 0 && (
              <tr>
                <td colSpan={5} className="text-center py-8 text-slate-400">
                  No proteins found matching "{filter}"
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>
    </div>
  );
};

export default DataTable;
