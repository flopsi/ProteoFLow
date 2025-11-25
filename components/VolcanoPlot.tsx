import React, { useMemo } from 'react';
import {
  ScatterChart,
  Scatter,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  ReferenceLine,
  Cell
} from 'recharts';
import { ProteinData, AnalysisConfig } from '../types';

interface VolcanoPlotProps {
  data: ProteinData[];
  config: AnalysisConfig;
  onPointClick: (point: ProteinData) => void;
}

const VolcanoPlot: React.FC<VolcanoPlotProps> = ({ data, config, onPointClick }) => {
  
  // Memoize chart data to improve performance
  const chartData = useMemo(() => {
    return data.map(d => ({
      x: d.log2FoldChange,
      y: d.negLog10PValue,
      ...d
    }));
  }, [data]);

  const CustomTooltip = ({ active, payload }: any) => {
    if (active && payload && payload.length) {
      const d = payload[0].payload as ProteinData;
      return (
        <div className="bg-white border border-slate-200 p-3 rounded shadow-lg text-sm z-50">
          <p className="font-bold text-slate-800">{d.gene}</p>
          <p className="text-slate-500 text-xs mb-2">{d.description}</p>
          <p className="text-slate-600">Log2 FC: <span className="font-mono text-slate-800">{d.log2FoldChange.toFixed(2)}</span></p>
          <p className="text-slate-600">-Log10 P: <span className="font-mono text-slate-800">{d.negLog10PValue.toFixed(2)}</span></p>
          <p className={`text-xs font-semibold mt-1 ${
             d.significance === 'UP' ? 'text-red-500' : 
             d.significance === 'DOWN' ? 'text-blue-500' : 'text-slate-400'
          }`}>
            {d.significance === 'UP' ? 'Upregulated' : d.significance === 'DOWN' ? 'Downregulated' : 'Not Significant'}
          </p>
        </div>
      );
    }
    return null;
  };

  return (
    <div className="w-full h-[500px] bg-white rounded-xl shadow-sm border border-slate-200 p-4">
      <h3 className="text-lg font-semibold text-slate-800 mb-2">Volcano Plot</h3>
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
            dataKey="x" 
            name="Log2 Fold Change" 
            label={{ value: 'Log2 Fold Change', position: 'bottom', offset: 0, fill: '#64748b' }}
            stroke="#94a3b8"
            tick={{ fill: '#64748b' }}
          />
          <YAxis 
            type="number" 
            dataKey="y" 
            name="-Log10 P-value" 
            label={{ value: '-Log10 P-value', angle: -90, position: 'insideLeft', fill: '#64748b' }} 
            stroke="#94a3b8"
            tick={{ fill: '#64748b' }}
          />
          <Tooltip content={<CustomTooltip />} cursor={{ strokeDasharray: '3 3' }} />
          
          <ReferenceLine x={config.fcCutoff} stroke="#94a3b8" strokeDasharray="3 3" />
          <ReferenceLine x={-config.fcCutoff} stroke="#94a3b8" strokeDasharray="3 3" />
          <ReferenceLine y={config.pValCutoff} stroke="#94a3b8" strokeDasharray="3 3" />

          <Scatter name="Proteins" data={chartData} fill="#8884d8">
            {chartData.map((entry, index) => {
               let fill = '#94a3b8'; // Neutral
               if (entry.log2FoldChange >= config.fcCutoff && entry.negLog10PValue >= config.pValCutoff) {
                   fill = '#ef4444'; // Red (Up)
               } else if (entry.log2FoldChange <= -config.fcCutoff && entry.negLog10PValue >= config.pValCutoff) {
                   fill = '#3b82f6'; // Blue (Down)
               }
               return <Cell key={`cell-${index}`} fill={fill} fillOpacity={fill === '#94a3b8' ? 0.3 : 0.8} />;
            })}
          </Scatter>
        </ScatterChart>
      </ResponsiveContainer>
    </div>
  );
};

export default VolcanoPlot;