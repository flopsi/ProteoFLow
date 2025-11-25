import React, { useState, useMemo } from 'react';
import { 
  LayoutDashboard, 
  UploadCloud, 
  FileText, 
  Settings, 
  MessageSquare, 
  Zap,
  Download,
  TestTube,
  BarChart,
  Activity,
  Microscope,
  CheckCircle2,
  AlertTriangle,
  ArrowRight
} from 'lucide-react';
import ReactMarkdown from 'react-markdown';

import { ProteinData, AnalysisConfig, AppView, ChatMessage, Species } from './types';
import { DEFAULT_CONFIG, MOCK_DATA_COUNT, SAMPLE_DESCRIPTIONS, SPECIES_CONFIG } from './constants';
import VolcanoPlot from './components/VolcanoPlot';
import LFQScatterPlot from './components/LFQScatterPlot';
import DataStats from './components/DataStats';
import DataTable from './components/DataTable';
import QCPlots from './components/QCPlots';
import { analyzeProteins, chatWithData } from './services/geminiService';
import { checkSkewness, calculateCV, mean } from './utils/stats';

// -- Mock Data Generator with Raw Replicates --
const generateMockData = (count: number): ProteinData[] => {
  const proteins: ProteinData[] = [];
  
  // Ratios: Human 65%, Yeast 20%, E.coli 15%
  const humanCount = Math.floor(count * 0.65);
  const yeastCount = Math.floor(count * 0.20);
  const ecoliCount = count - humanCount - yeastCount;

  const createProtein = (i: number, species: Species, startIdx: number): ProteinData => {
    // Generate base intensity (Log Normal)
    // We generate RAW intensity first (skewed distribution)
    const baseMean = Math.pow(10, 4 + Math.random() * 4); // 10^4 to 10^8
    
    const expectedFC = SPECIES_CONFIG[species].expectedLog2FC;
    // For raw data generation, we convert log2 fc back to linear ratio
    // If expectedLog2FC is 1 (2x), then ratio is 2.
    const foldChangeMultiplier = Math.pow(2, expectedFC);
    
    // Condition B is base, Condition A is changed
    const meanB = baseMean;
    const meanA = baseMean * foldChangeMultiplier;

    // Generate Replicates with noise and occasional missing values
    const genReps = (mu: number) => {
        const reps = [];
        for(let k=0; k<3; k++) {
            // Add multiplicative noise (log-normal noise) to simulate raw instrument data
            const noise = Math.pow(2, (Math.random() - 0.5) * 0.8);
            let val = mu * noise;
            
            // 5% chance of missing value (set to 0)
            if (Math.random() < 0.05) val = 0;
            reps.push(val);
        }
        return reps;
    };

    const repsA = genReps(meanA);
    const repsB = genReps(meanB);
    
    // Initial calculation (on raw data, will be updated during transform)
    const validA = repsA.filter(r => r > 0);
    const validB = repsB.filter(r => r > 0);
    const avgA = validA.length ? validA.reduce((a,b)=>a+b,0)/validA.length : 0;
    const avgB = validB.length ? validB.reduce((a,b)=>a+b,0)/validB.length : 0;
    
    const fc = (avgB > 0) ? avgA / avgB : 0;
    const log2fc = (fc > 0) ? Math.log2(fc) : 0;
    
    // P-value sim (dummy)
    let pVal = 0.05 / (Math.abs(log2fc) + 0.1) * Math.random();
    if (species === 'Human') pVal = Math.random();

    const genePrefix = species === 'Human' ? 'HUM' : species === 'Yeast' ? 'YEA' : 'ECO';

    return {
      id: `ID_${startIdx + i}`,
      gene: `${genePrefix}_${Math.floor(Math.random()*9000)+1000}`,
      species: species,
      description: `${species} protein ${i}`,
      replicates: { A: repsA, B: repsB },
      foldChange: fc,
      log2FoldChange: log2fc,
      pValue: pVal,
      negLog10PValue: -Math.log10(pVal),
      averageIntensityLog10: Math.log10((avgA + avgB) / 2 || 1),
      significance: 'NS',
      cvA: calculateCV(validA),
      cvB: calculateCV(validB)
    };
  };

  for (let i = 0; i < humanCount; i++) proteins.push(createProtein(i, 'Human', 0));
  for (let i = 0; i < yeastCount; i++) proteins.push(createProtein(i, 'Yeast', humanCount));
  for (let i = 0; i < ecoliCount; i++) proteins.push(createProtein(i, 'E.coli', humanCount + yeastCount));

  return proteins.sort(() => Math.random() - 0.5);
};

const App: React.FC = () => {
  // State
  const [view, setView] = useState<AppView>(AppView.UPLOAD);
  const [data, setData] = useState<ProteinData[]>([]);
  const [config, setConfig] = useState<AnalysisConfig>({...DEFAULT_CONFIG, imputationEnabled: false});
  const [selectedProtein, setSelectedProtein] = useState<ProteinData | null>(null);
  const [activePlot, setActivePlot] = useState<'MA' | 'VOLCANO'>('MA');
  
  // QC State
  const [normalityCheck, setNormalityCheck] = useState<{skewness: number, needsTransform: boolean} | null>(null);
  const [isTransformed, setIsTransformed] = useState(false);
  
  // AI & Chat
  const [aiReport, setAiReport] = useState<string>('');
  const [isAnalyzing, setIsAnalyzing] = useState(false);
  const [experimentContext, setExperimentContext] = useState(SAMPLE_DESCRIPTIONS[0]);
  const [chatHistory, setChatHistory] = useState<ChatMessage[]>([]);
  const [chatInput, setChatInput] = useState('');
  const [isChatting, setIsChatting] = useState(false);

  // Apply Significance Thresholds
  const processedData = useMemo(() => {
    return data.map(d => {
      let sig: 'UP' | 'DOWN' | 'NS' = 'NS';
      if (d.negLog10PValue >= config.pValCutoff) {
        if (d.log2FoldChange >= config.fcCutoff) sig = 'UP';
        else if (d.log2FoldChange <= -config.fcCutoff) sig = 'DOWN';
      }
      return { ...d, significance: sig };
    });
  }, [data, config]);

  // Handlers
  const handleLoadData = () => {
    const mock = generateMockData(MOCK_DATA_COUNT);
    setData(mock);
    setNormalityCheck(null);
    setIsTransformed(false);
    setView(AppView.QC);
  };

  const handleNormalityCheck = () => {
    // Check a subset of proteins for skewness in raw data
    // We check the first replicate of Condition A for the first 100 proteins
    const sampleValues = data.slice(0, 100).map(d => d.replicates.A[0]).filter(v => v > 0);
    const check = checkSkewness(sampleValues);
    setNormalityCheck({ skewness: check.skewness, needsTransform: check.isSkewed });
  };

  const handleTransform = () => {
    // Apply Log2 Transformation to replicates and recalculate stats
    const transformedData = data.map(d => {
        // Log2 transform replicates
        const log2RepsA = d.replicates.A.map(v => v > 0 ? Math.log2(v) : 0);
        const log2RepsB = d.replicates.B.map(v => v > 0 ? Math.log2(v) : 0);
        
        // Recalculate means of Log2 values
        const validA = log2RepsA.filter(v => v > 0);
        const validB = log2RepsB.filter(v => v > 0);
        
        const meanA = validA.length ? mean(validA) : 0;
        const meanB = validB.length ? mean(validB) : 0;
        
        // Log2 Fold Change is now simply Difference of Means (MeanA - MeanB)
        const log2fc = (meanA > 0 && meanB > 0) ? meanA - meanB : 0;
        
        // Intensity is sum of means (or avg)
        const avgIntensity = (meanA + meanB) / 2;

        return {
            ...d,
            replicates: { A: log2RepsA, B: log2RepsB },
            log2FoldChange: log2fc,
            averageIntensityLog10: avgIntensity, // Used as x-axis for MA plot
            // Recalculate CV on the Log values (though CV is typically for linear, we update it for display)
            cvA: calculateCV(validA),
            cvB: calculateCV(validB)
        };
    });
    
    setData(transformedData);
    setIsTransformed(true);
  };

  const handleGenerateReport = async () => {
    setIsAnalyzing(true);
    const report = await analyzeProteins(processedData, experimentContext);
    setAiReport(report);
    setIsAnalyzing(false);
    setView(AppView.REPORT);
  };

  const handleSendMessage = async () => {
    if (!chatInput.trim()) return;
    const userMsg: ChatMessage = { role: 'user', text: chatInput, timestamp: Date.now() };
    setChatHistory(prev => [...prev, userMsg]);
    setChatInput('');
    setIsChatting(true);

    const contextStats = `
      Total Proteins: ${data.length}. 
      Significant Up: ${processedData.filter(d => d.significance === 'UP').length}.
      Significant Down: ${processedData.filter(d => d.significance === 'DOWN').length}.
    `;
    
    const responseText = await chatWithData(
        chatHistory.map(h => ({role: h.role, text: h.text})), 
        userMsg.text, 
        contextStats
    );

    const modelMsg: ChatMessage = { role: 'model', text: responseText, timestamp: Date.now() };
    setChatHistory(prev => [...prev, modelMsg]);
    setIsChatting(false);
  };

  // -- Views --

  const renderUpload = () => (
    <div className="flex flex-col items-center justify-center h-full p-8 bg-white rounded-2xl shadow-sm border border-slate-200">
      <div className="w-20 h-20 bg-indigo-50 rounded-full flex items-center justify-center mb-6">
        <UploadCloud size={40} className="text-indigo-600" />
      </div>
      <h2 className="text-2xl font-bold text-slate-800 mb-2">Upload Proteomics Data</h2>
      <p className="text-slate-500 text-center max-w-md mb-8">
        Upload your MaxQuant output or CSV file containing replicate intensities.
        Compatible with LFQ-Bench datasets.
      </p>
      
      <div className="flex gap-4">
        <button className="px-6 py-3 bg-white border border-slate-300 text-slate-700 font-medium rounded-lg hover:bg-slate-50 transition-colors">
          Select File
        </button>
        <button 
          onClick={handleLoadData}
          className="px-6 py-3 bg-indigo-600 text-white font-medium rounded-lg hover:bg-indigo-700 transition-colors shadow-lg shadow-indigo-200 flex items-center gap-2"
        >
          <Zap size={18} /> Load Demo HYE Data
        </button>
      </div>
    </div>
  );

  const renderQC = () => (
    <div className="space-y-6">
        <div className="flex justify-between items-center">
            <h2 className="text-2xl font-bold text-slate-800 flex items-center gap-2">
                <Microscope className="text-indigo-600" /> Quality Control & Normalization
            </h2>
            <div className="flex gap-3">
               {!normalityCheck ? (
                   <button 
                     onClick={handleNormalityCheck}
                     className="px-4 py-2 bg-amber-100 text-amber-800 rounded-lg hover:bg-amber-200 font-medium flex items-center gap-2"
                   >
                       <Activity size={18} /> Check Normality
                   </button>
               ) : (
                   <div className={`px-4 py-2 rounded-lg font-medium flex items-center gap-2 ${normalityCheck.needsTransform ? 'bg-red-50 text-red-700' : 'bg-green-50 text-green-700'}`}>
                       {normalityCheck.needsTransform ? <AlertTriangle size={18}/> : <CheckCircle2 size={18} />}
                       Skewness: {normalityCheck.skewness.toFixed(2)}
                   </div>
               )}
               
               {normalityCheck?.needsTransform && !isTransformed && (
                   <button 
                     onClick={handleTransform}
                     className="px-4 py-2 bg-indigo-600 text-white rounded-lg hover:bg-indigo-700 font-medium flex items-center gap-2"
                   >
                       <Zap size={18} /> Apply Log2 Transform
                   </button>
               )}

               <button 
                 onClick={() => setView(AppView.DASHBOARD)}
                 disabled={!isTransformed}
                 className={`px-4 py-2 rounded-lg font-medium flex items-center gap-2 transition-colors ${
                     isTransformed 
                     ? 'bg-emerald-600 text-white hover:bg-emerald-700' 
                     : 'bg-slate-100 text-slate-400 cursor-not-allowed'
                 }`}
               >
                   Proceed to Analysis <ArrowRight size={18} />
               </button>
            </div>
        </div>

        <QCPlots data={data} isTransformed={isTransformed} />
    </div>
  );

  const renderDashboard = () => (
    <div className="space-y-6 animate-in fade-in duration-500">
      {/* Top Stats */}
      <DataStats data={processedData} />

      {/* Main Charts Area */}
      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        <div className="lg:col-span-2 space-y-6">
            <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-1">
                <div className="flex border-b border-slate-100">
                    <button 
                        onClick={() => setActivePlot('MA')}
                        className={`flex-1 py-3 text-sm font-medium transition-colors ${activePlot === 'MA' ? 'text-indigo-600 border-b-2 border-indigo-600' : 'text-slate-500 hover:text-slate-700'}`}
                    >
                        Ratio Plot (MA)
                    </button>
                    <button 
                        onClick={() => setActivePlot('VOLCANO')}
                        className={`flex-1 py-3 text-sm font-medium transition-colors ${activePlot === 'VOLCANO' ? 'text-indigo-600 border-b-2 border-indigo-600' : 'text-slate-500 hover:text-slate-700'}`}
                    >
                        Volcano Plot
                    </button>
                </div>
                <div className="p-4">
                     {activePlot === 'MA' ? (
                         <LFQScatterPlot data={processedData} onPointClick={setSelectedProtein} />
                     ) : (
                         <VolcanoPlot data={processedData} config={config} onPointClick={setSelectedProtein} />
                     )}
                </div>
            </div>
            
            <DataTable data={processedData} onSelect={setSelectedProtein} />
        </div>

        <div className="space-y-6">
          {/* Selected Protein Details */}
          <div className="bg-white p-5 rounded-xl shadow-sm border border-slate-200">
            <h3 className="text-lg font-bold text-slate-800 mb-4 flex items-center gap-2">
              <Microscope size={20} className="text-indigo-600" />
              Protein Details
            </h3>
            {selectedProtein ? (
              <div className="space-y-3">
                <div className="p-3 bg-slate-50 rounded-lg border border-slate-100">
                  <span className="text-xs font-semibold text-slate-500 uppercase">Gene ID</span>
                  <div className="text-lg font-bold text-slate-800">{selectedProtein.gene}</div>
                  <div className="text-xs text-slate-500">{selectedProtein.id}</div>
                </div>
                
                <div className="grid grid-cols-2 gap-3">
                   <div className="p-3 bg-slate-50 rounded-lg">
                      <span className="text-xs text-slate-500">Log2 FC</span>
                      <div className={`text-xl font-bold ${selectedProtein.log2FoldChange > 0 ? 'text-emerald-600' : 'text-blue-600'}`}>
                          {selectedProtein.log2FoldChange.toFixed(2)}
                      </div>
                   </div>
                   <div className="p-3 bg-slate-50 rounded-lg">
                      <span className="text-xs text-slate-500">-Log10 P</span>
                      <div className="text-xl font-bold text-slate-800">
                          {selectedProtein.negLog10PValue.toFixed(2)}
                      </div>
                   </div>
                </div>

                <div className="p-3 bg-slate-50 rounded-lg">
                     <span className="text-xs text-slate-500 mb-1 block">Replicate CVs</span>
                     <div className="flex justify-between text-sm">
                         <span>Cond A: <span className="font-mono font-bold">{selectedProtein.cvA.toFixed(1)}%</span></span>
                         <span>Cond B: <span className="font-mono font-bold">{selectedProtein.cvB.toFixed(1)}%</span></span>
                     </div>
                </div>
              </div>
            ) : (
              <div className="text-center py-8 text-slate-400 text-sm">
                Select a data point from the plots or table to view details.
              </div>
            )}
          </div>

          {/* AI Chat Assistant */}
          <div className="bg-white rounded-xl shadow-sm border border-slate-200 flex flex-col h-[400px]">
            <div className="p-4 border-b border-slate-200 bg-slate-50 rounded-t-xl">
               <h3 className="font-bold text-slate-800 flex items-center gap-2">
                   <MessageSquare size={18} className="text-indigo-600" /> AI Assistant
               </h3>
            </div>
            <div className="flex-1 overflow-y-auto p-4 space-y-3">
               {chatHistory.length === 0 && (
                   <p className="text-xs text-slate-400 text-center mt-10">
                       Ask me about normalization, outliers, or specific proteins.
                   </p>
               )}
               {chatHistory.map((msg, i) => (
                   <div key={i} className={`flex ${msg.role === 'user' ? 'justify-end' : 'justify-start'}`}>
                       <div className={`max-w-[85%] rounded-lg p-3 text-sm ${
                           msg.role === 'user' 
                           ? 'bg-indigo-600 text-white' 
                           : 'bg-slate-100 text-slate-800'
                       }`}>
                           {msg.text}
                       </div>
                   </div>
               ))}
               {isChatting && <div className="text-xs text-slate-400 animate-pulse ml-2">Thinking...</div>}
            </div>
            <div className="p-3 border-t border-slate-200">
                <div className="flex gap-2">
                    <input 
                        className="flex-1 text-sm border border-slate-300 rounded-md px-3 py-2 focus:outline-none focus:ring-2 focus:ring-indigo-500"
                        placeholder="Ask a question..."
                        value={chatInput}
                        onChange={(e) => setChatInput(e.target.value)}
                        onKeyDown={(e) => e.key === 'Enter' && handleSendMessage()}
                    />
                    <button 
                        onClick={handleSendMessage}
                        disabled={!chatInput.trim() || isChatting}
                        className="bg-indigo-600 text-white p-2 rounded-md hover:bg-indigo-700 disabled:opacity-50"
                    >
                        <ArrowRight size={16} />
                    </button>
                </div>
            </div>
          </div>
        </div>
      </div>
      
      <div className="flex justify-end pt-8">
           <button 
             onClick={handleGenerateReport}
             className="px-6 py-3 bg-slate-800 text-white rounded-xl hover:bg-slate-900 shadow-lg flex items-center gap-2 font-medium"
           >
              {isAnalyzing ? <span className="animate-spin">⌛</span> : <FileText size={20} />}
              Generate Full Report
           </button>
      </div>
    </div>
  );

  const renderReport = () => (
      <div className="max-w-4xl mx-auto bg-white p-10 rounded-xl shadow-sm border border-slate-200 min-h-[600px]">
          <div className="flex justify-between items-center mb-8 pb-4 border-b border-slate-100">
              <h1 className="text-3xl font-bold text-slate-800">Proteomics Analysis Report</h1>
              <button onClick={() => setView(AppView.DASHBOARD)} className="text-slate-500 hover:text-indigo-600">
                  Back to Dashboard
              </button>
          </div>
          <div className="prose prose-slate max-w-none">
              <ReactMarkdown>{aiReport}</ReactMarkdown>
          </div>
      </div>
  );

  return (
    <div className="min-h-screen flex flex-col">
      {/* Header */}
      <header className="bg-white border-b border-slate-200 sticky top-0 z-30">
        <div className="max-w-7xl mx-auto px-6 h-16 flex items-center justify-between">
          <div className="flex items-center gap-3">
            <div className="w-8 h-8 bg-gradient-to-br from-indigo-500 to-purple-600 rounded-lg flex items-center justify-center text-white font-bold shadow-md">
              P
            </div>
            <span className="text-xl font-bold bg-clip-text text-transparent bg-gradient-to-r from-indigo-700 to-purple-600">
              ProteoFlow <span className="text-xs font-normal text-slate-400 ml-1">LFQ Bench</span>
            </span>
          </div>
          
          <nav className="flex gap-1 bg-slate-100 p-1 rounded-lg">
             {[AppView.UPLOAD, AppView.QC, AppView.DASHBOARD, AppView.REPORT].map((v) => (
                 <button
                    key={v}
                    onClick={() => { if (data.length > 0) setView(v); }}
                    disabled={data.length === 0}
                    className={`px-4 py-1.5 text-xs font-medium rounded-md transition-all ${
                        view === v 
                        ? 'bg-white text-indigo-700 shadow-sm' 
                        : 'text-slate-500 hover:text-slate-700'
                    } ${data.length === 0 ? 'opacity-50 cursor-not-allowed' : ''}`}
                 >
                     {v}
                 </button>
             ))}
          </nav>

          <div className="flex items-center gap-4">
             <div className="text-xs text-slate-500 hidden md:block">
                 {data.length > 0 ? `${data.length} Proteins Loaded` : 'No Data'}
             </div>
             <Settings className="text-slate-400 hover:text-slate-600 cursor-pointer" size={20} />
          </div>
        </div>
      </header>

      {/* Main Content */}
      <main className="flex-1 bg-slate-50 p-6 overflow-auto">
        <div className="max-w-7xl mx-auto h-full">
          {view === AppView.UPLOAD && renderUpload()}
          {view === AppView.QC && renderQC()}
          {view === AppView.DASHBOARD && renderDashboard()}
          {view === AppView.REPORT && renderReport()}
        </div>
      </main>
    </div>
  );
};

export default App;