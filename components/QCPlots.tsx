import React, { useMemo } from 'react';
import { ProteinData, BoxPlotStat, PCAPoint } from '../types';
import { calculateBoxPlotStats, calculatePCA } from '../utils/stats';
import {
  BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer,
  ScatterChart, Scatter, Cell, LineChart, Line, ReferenceLine
} from 'recharts';

interface QCPlotsProps {
  data: ProteinData[];
  isTransformed: boolean;
}

const QCPlots: React.FC<QCPlotsProps> = ({ data, isTransformed }) => {

  // --- 1. Box Plots (Distribution of Intensities) ---
  const boxPlotData = useMemo(() => {
    // Flatten data for each replicate
    const samples = ['A1', 'A2', 'A3', 'B1', 'B2', 'B3'];
    return samples.map((sample, idx) => {
      const group = idx < 3 ? 'A' : 'B';
      const repIdx = idx % 3;
      const values = data.map(p => p.replicates[group as 'A' | 'B'][repIdx]).filter(v => v > 0);
      return calculateBoxPlotStats(values, sample);
    });
  }, [data]);

  // --- 2. CV Analysis (Coefficient of Variation) ---
  const cvHistogram = useMemo(() => {
    // Create bins for histogram
    const bins = Array(20).fill(0); // 0-5%, 5-10%... up to 100%
    data.forEach(p => {
       const binIdxA = Math.min(Math.floor(p.cvA / 5), 19);
       const binIdxB = Math.min(Math.floor(p.cvB / 5), 19);
       bins[binIdxA]++;
       bins[binIdxB]++;
    });
    return bins.map((count, i) => ({ bin: `${i*5}-${(i+1)*5}%`, count }));
  }, [data]);

  // --- 3. PCA Analysis ---
  const pcaPoints = useMemo(() => calculatePCA(data), [data]);

  // --- 4. Rank Plots (Dynamic Range) ---
  const rankData = useMemo(() => {
    // Sort by intensity
    const sortedA = [...data].sort((a, b) => b.averageIntensityLog10 - a.averageIntensityLog10); // High to Low
    
    // Sample down for performance if huge
    const step = Math.ceil(sortedA.length / 500);
    return sortedA.filter((_, i) => i % step === 0).map((p, i) => ({
      rank: i * step,
      intensity: p.averageIntensityLog10
    }));
  }, [data]);

  // --- 5. Missing Values Heatmap (Simplified SVG) ---
  // Create a 50x30 grid (approx 1500 points)
  const heatmapGrid = useMemo(() => {
      const gridSize = Math.ceil(Math.sqrt(data.length));
      return data.slice(0, gridSize*gridSize).map((p, i) => {
          const hasMissing = p.replicates.A.some(v => v === 0) || p.replicates.B.some(v => v === 0);
          return { i, hasMissing };
      });
  }, [data]);

  const CustomBoxPlot = ({ stats }: { stats: BoxPlotStat }) => (
      <div className="flex flex-col items-center mx-2 w-12 group relative">
          {/* Whiskers */}
          <div className="w-px bg-slate-400 absolute left-1/2 -translate-x-1/2" 
               style={{ top: '0%', bottom: '0%' }}></div> {/* Ideally scaled to container */}
          
          {/* Box */}
          <div className="w-full bg-indigo-100 border border-indigo-500 relative h-32 rounded-sm hover:bg-indigo-200 transition-colors">
             {/* This is a visual dummy. Real boxplot needs absolute scaling based on Y-axis. 
                 Since Recharts boxplot is tricky, we are using the stats to render a summary below. */}
          </div>
          <span className="mt-2 text-xs font-bold text-slate-600">{stats.label}</span>
      </div>
  );

  return (
    <div className="space-y-6 pb-12">
      
      {/* Row 1: Normality & Boxplots */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <div className="bg-white p-4 rounded-xl shadow-sm border border-slate-200">
              <h3 className="font-bold text-slate-700 mb-4">Replicate Distribution (Boxplot Stats)</h3>
              <p className="text-xs text-slate-500 mb-4">
                  Displays Median, Q1, Q3. {isTransformed ? '(Log2 Transformed)' : '(Raw Intensity)'}
              </p>
              <div className="h-64 flex items-end justify-around pb-6 border-b border-slate-200 relative">
                  {/* Custom CSS implementation of Boxplot visualizer */}
                  {boxPlotData.map((stat) => {
                       // Normalize height for visualization
                       const globalMax = Math.max(...boxPlotData.map(s => s.max));
                       const globalMin = Math.min(...boxPlotData.map(s => s.min));
                       const range = globalMax - globalMin;
                       
                       const getH = (val: number) => ((val - globalMin) / range) * 100;
                       
                       return (
                           <div key={stat.label} className="relative h-full w-12 group">
                               {/* Whiskers Line */}
                               <div className="absolute left-1/2 w-px bg-slate-800 -translate-x-1/2" 
                                    style={{ bottom: `${getH(stat.min)}%`, height: `${getH(stat.max) - getH(stat.min)}%` }} />
                               {/* Whiskers Caps */}
                               <div className="absolute left-1/4 w-1/2 h-px bg-slate-800" style={{ bottom: `${getH(stat.min)}%` }} />
                               <div className="absolute left-1/4 w-1/2 h-px bg-slate-800" style={{ bottom: `${getH(stat.max)}%` }} />
                               
                               {/* Box */}
                               <div className="absolute w-full bg-indigo-200 border border-indigo-600 opacity-80"
                                    style={{ bottom: `${getH(stat.q1)}%`, height: `${getH(stat.q3) - getH(stat.q1)}%` }} />
                               
                               {/* Median */}
                               <div className="absolute w-full h-0.5 bg-indigo-900 z-10" style={{ bottom: `${getH(stat.median)}%` }} />
                               
                               <div className="absolute -bottom-6 w-full text-center text-xs font-bold text-slate-600">
                                   {stat.label}
                               </div>
                               
                               {/* Tooltip */}
                               <div className="opacity-0 group-hover:opacity-100 absolute -top-12 left-1/2 -translate-x-1/2 bg-slate-800 text-white text-xs p-2 rounded z-20 pointer-events-none whitespace-nowrap">
                                   Median: {stat.median.toFixed(2)}
                               </div>
                           </div>
                       )
                  })}
              </div>
          </div>

          <div className="bg-white p-4 rounded-xl shadow-sm border border-slate-200">
              <h3 className="font-bold text-slate-700 mb-2">CV Analysis</h3>
              <p className="text-xs text-slate-500 mb-4">Coefficient of Variation Distribution across all proteins</p>
              <ResponsiveContainer width="100%" height={250}>
                  <BarChart data={cvHistogram}>
                      <CartesianGrid strokeDasharray="3 3" vertical={false} />
                      <XAxis dataKey="bin" tick={{fontSize: 10}} label={{ value: '% CV', position: 'bottom', offset: 0 }} />
                      <YAxis tick={{fontSize: 10}} />
                      <Tooltip />
                      <Bar dataKey="count" fill="#8884d8" name="Count" />
                  </BarChart>
              </ResponsiveContainer>
              <div className="mt-2 text-center">
                  <span className="text-xs text-slate-500 bg-slate-100 px-2 py-1 rounded">
                      Median CV: {
                          (data.reduce((acc, p) => acc + p.cvA + p.cvB, 0) / (data.length * 2)).toFixed(1)
                      }%
                  </span>
              </div>
          </div>
      </div>

      {/* Row 2: PCA & Rank */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
           <div className="bg-white p-4 rounded-xl shadow-sm border border-slate-200">
               <h3 className="font-bold text-slate-700 mb-2">Principal Component Analysis</h3>
               <ResponsiveContainer width="100%" height={300}>
                   <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 20 }}>
                       <CartesianGrid />
                       <XAxis type="number" dataKey="x" name="PC1" />
                       <YAxis type="number" dataKey="y" name="PC2" />
                       <Tooltip cursor={{ strokeDasharray: '3 3' }} />
                       <Scatter name="Samples" data={pcaPoints}>
                           {pcaPoints.map((entry, index) => (
                               <Cell key={`cell-${index}`} fill={entry.group === 'A' ? '#ef4444' : '#3b82f6'} />
                           ))}
                       </Scatter>
                   </ScatterChart>
               </ResponsiveContainer>
               <div className="flex justify-center gap-4 text-xs mt-2">
                   <div className="flex items-center"><div className="w-3 h-3 bg-red-500 rounded-full mr-1"></div> Condition A</div>
                   <div className="flex items-center"><div className="w-3 h-3 bg-blue-500 rounded-full mr-1"></div> Condition B</div>
               </div>
           </div>

           <div className="bg-white p-4 rounded-xl shadow-sm border border-slate-200">
               <h3 className="font-bold text-slate-700 mb-2">Rank Plot (Dynamic Range)</h3>
               <ResponsiveContainer width="100%" height={300}>
                   <LineChart data={rankData}>
                       <CartesianGrid strokeDasharray="3 3" />
                       <XAxis dataKey="rank" label={{ value: 'Protein Rank', position: 'bottom' }} />
                       <YAxis label={{ value: 'Log10 Intensity', angle: -90, position: 'insideLeft' }} />
                       <Tooltip />
                       <Line type="monotone" dataKey="intensity" stroke="#82ca9d" dot={false} strokeWidth={2} />
                   </LineChart>
               </ResponsiveContainer>
           </div>
      </div>

      {/* Row 3: Missing Values Heatmap */}
      <div className="bg-white p-4 rounded-xl shadow-sm border border-slate-200">
           <h3 className="font-bold text-slate-700 mb-2">Missing Value Map</h3>
           <p className="text-xs text-slate-500 mb-4">Yellow = Missing Value (NaN or 0) in at least one replicate.</p>
           <div className="w-full overflow-hidden flex flex-wrap gap-0.5">
               {heatmapGrid.map((cell) => (
                   <div 
                      key={cell.i} 
                      className={`w-3 h-3 ${cell.hasMissing ? 'bg-amber-300' : 'bg-slate-800'} rounded-sm`}
                      title={`Protein ${cell.i} ${cell.hasMissing ? '- Missing Data' : ''}`}
                   ></div>
               ))}
           </div>
      </div>

    </div>
  );
};

export default QCPlots;