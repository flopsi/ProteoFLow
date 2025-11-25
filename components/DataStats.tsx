import React from 'react';
import { ProteinData } from '../types';
import { SPECIES_CONFIG } from '../constants';
import { Calculator, CheckCircle2, AlertCircle } from 'lucide-react';

interface DataStatsProps {
  data: ProteinData[];
}

const DataStats: React.FC<DataStatsProps> = ({ data }) => {
  const getMedian = (values: number[]) => {
    if (values.length === 0) return 0;
    values.sort((a, b) => a - b);
    const half = Math.floor(values.length / 2);
    if (values.length % 2) return values[half];
    return (values[half - 1] + values[half]) / 2.0;
  };

  const speciesStats = Object.keys(SPECIES_CONFIG).map(speciesKey => {
      const speciesName = speciesKey as keyof typeof SPECIES_CONFIG;
      const subset = data.filter(d => d.species === speciesName);
      const medianFC = getMedian(subset.map(d => d.log2FoldChange));
      const expected = SPECIES_CONFIG[speciesName].expectedLog2FC;
      const count = subset.length;
      
      // Simple accuracy check
      const diff = Math.abs(medianFC - expected);
      const isGood = diff < 0.3; // Arbitrary tolerance for demo

      return {
          species: speciesName,
          count,
          medianFC,
          expected,
          isGood,
          config: SPECIES_CONFIG[speciesName]
      };
  });

  return (
    <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-6">
      {speciesStats.map((stat) => (
        <div key={stat.species} className="bg-white p-4 rounded-xl shadow-sm border border-slate-200 flex flex-col relative overflow-hidden">
          <div className="absolute top-0 left-0 w-1 h-full" style={{ backgroundColor: stat.config.color }}></div>
          
          <div className="flex justify-between items-start mb-2">
              <div>
                  <h4 className="text-sm font-semibold text-slate-600 uppercase tracking-wide">{stat.species}</h4>
                  <p className="text-xs text-slate-400">{stat.count} Proteins</p>
              </div>
              {stat.isGood ? (
                  <CheckCircle2 size={18} className="text-emerald-500" />
              ) : (
                  <AlertCircle size={18} className="text-amber-500" />
              )}
          </div>
          
          <div className="mt-2 flex items-baseline justify-between">
              <div>
                  <span className="text-2xl font-bold text-slate-800">{stat.medianFC.toFixed(2)}</span>
                  <span className="text-xs text-slate-500 ml-1">Observed Log2FC</span>
              </div>
              <div className="text-right">
                  <span className="text-sm font-mono text-slate-400">Exp: {stat.expected}</span>
              </div>
          </div>
          
          <div className="mt-3 w-full bg-slate-100 rounded-full h-1.5 overflow-hidden">
             <div 
                className="h-full rounded-full transition-all duration-500" 
                style={{ 
                    width: '100%', 
                    backgroundColor: stat.config.color,
                    opacity: 0.5
                }}
             ></div>
          </div>
        </div>
      ))}
    </div>
  );
};

export default DataStats;
