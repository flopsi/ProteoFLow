import React from 'react';
import {
  ScatterChart,
  Scatter,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  ReferenceLine,
  Legend
} from 'recharts';
import { ProteinData, Species } from '../types';
import { SPECIES_CONFIG } from '../constants';

interface LFQScatterPlotProps {
  data: ProteinData[];
  onPointClick: (point: ProteinData) => void;
}

const LFQScatterPlot: React.FC<LFQScatterPlotProps> = ({ data, onPointClick }) => {
  
  // Separate data by species for the chart series
  const humanData = data.filter(d => d.species === 'Human');
  const yeastData = data.filter(d => d.species === 'Yeast');
  const ecoliData = data.filter(d => d.species === 'E.coli');

  const CustomTooltip = ({ active, payload }: any) => {
    if (active && payload && payload.length) {
      const d = payload[0].payload as ProteinData;
      return (
        <div className="bg-white border border-slate-200 p-3 rounded shadow-lg text-sm z-50">
          <p className="font-bold text-slate-800">{d.gene} <span className="text-xs font-normal text-slate-500">({d.species})</span></p>
          <p className="text-slate-600">Log2 FC: <span className="font-mono text-slate-800">{d.log2FoldChange.toFixed(2)}</span></p>
          <p className="text-slate-600">Intensity (Log10): <span className="font-mono text-slate-800">{d.averageIntensityLog10.toFixed(2)}</span></p>
        </div>
      );
    }
    return null;
  };

  return (
    <div className="w-full h-[500px] bg-white rounded-xl shadow-sm border border-slate-200 p-4">
      <div className="flex justify-between items-center mb-2">
         <h3 className="text-lg font-semibold text-slate-800">Ratio vs Intensity (MA Plot)</h3>
         <div className="flex gap-4 text-xs text-slate-500">
            {Object.entries(SPECIES_CONFIG).map(([species, conf]) => (
                <div key={species} className="flex items-center">
                    <span className="w-3 h-3 rounded-full mr-1" style={{ backgroundColor: conf.color }}></span>
                    {species}: Exp. {conf.expectedLog2FC}
                </div>
            ))}
         </div>
      </div>
      <ResponsiveContainer width="100%" height="100%">
        <ScatterChart
          margin={{ top: 20, right: 20, bottom: 20, left: 20 }}
          onClick={(e: any) => {
             if (e && e.activePayload && e.activePayload[0]) {
                 onPointClick(e.activePayload[0].payload as ProteinData);
             }
          }}
        >
          <CartesianGrid strokeDasharray="3 3" stroke="#e2e8f0" />
          <XAxis 
            type="number" 
            dataKey="averageIntensityLog10" 
            name="Log10 Intensity" 
            domain={['auto', 'auto']}
            label={{ value: 'Log10 Average Intensity', position: 'bottom', offset: 0, fill: '#64748b' }}
            stroke="#94a3b8"
            tick={{ fill: '#64748b' }}
          />
          <YAxis 
            type="number" 
            dataKey="log2FoldChange" 
            name="Log2 Fold Change" 
            label={{ value: 'Log2 Fold Change', angle: -90, position: 'insideLeft', fill: '#64748b' }} 
            stroke="#94a3b8"
            tick={{ fill: '#64748b' }}
          />
          <Tooltip content={<CustomTooltip />} cursor={{ strokeDasharray: '3 3' }} />
          
          <ReferenceLine y={0} stroke="#cbd5e1" strokeWidth={2} />
          {/* Expected Ratio Lines */}
          <ReferenceLine y={SPECIES_CONFIG['Yeast'].expectedLog2FC} stroke={SPECIES_CONFIG['Yeast'].color} strokeDasharray="3 3" opacity={0.5} />
          <ReferenceLine y={SPECIES_CONFIG['E.coli'].expectedLog2FC} stroke={SPECIES_CONFIG['E.coli'].color} strokeDasharray="3 3" opacity={0.5} />

          <Scatter name="Human" data={humanData} fill={SPECIES_CONFIG['Human'].color} fillOpacity={0.4} />
          <Scatter name="Yeast" data={yeastData} fill={SPECIES_CONFIG['Yeast'].color} fillOpacity={0.6} />
          <Scatter name="E.coli" data={ecoliData} fill={SPECIES_CONFIG['E.coli'].color} fillOpacity={0.6} />
          
          <Legend />
        </ScatterChart>
      </ResponsiveContainer>
    </div>
  );
};

export default LFQScatterPlot;