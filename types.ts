export type Species = 'Human' | 'Yeast' | 'E.coli';

export interface ReplicateData {
  A: number[]; // e.g., [A1, A2, A3]
  B: number[]; // e.g., [B1, B2, B3]
}

export interface ProteinData {
  id: string;
  gene: string;
  species: Species;
  description: string;
  replicates: ReplicateData; // Raw intensities
  foldChange: number; // Calculated
  log2FoldChange: number; // Calculated
  pValue: number;
  negLog10PValue: number; // Calculated
  averageIntensityLog10: number; // Average of all valid replicates
  significance: 'UP' | 'DOWN' | 'NS';
  cvA: number; // Coefficient of Variation Condition A
  cvB: number; // Coefficient of Variation Condition B
}

export interface AnalysisConfig {
  fcCutoff: number; // log2 fold change cutoff (absolute)
  pValCutoff: number; // -log10 p-value cutoff
  imputationEnabled: boolean;
}

export enum AppView {
  UPLOAD = 'UPLOAD',
  QC = 'QC',
  DASHBOARD = 'DASHBOARD',
  REPORT = 'REPORT'
}

export interface ChatMessage {
  role: 'user' | 'model';
  text: string;
  timestamp: number;
}

export interface PCAPoint {
  sample: string;
  x: number;
  y: number;
  group: 'A' | 'B';
}

export interface BoxPlotStat {
  label: string;
  min: number;
  q1: number;
  median: number;
  q3: number;
  max: number;
}
