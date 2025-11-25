export const MOCK_DATA_COUNT = 1500;

export const DEFAULT_CONFIG = {
  fcCutoff: 0.58, // ~1.5 fold change
  pValCutoff: 1.3, // p < 0.05
};

export const SAMPLE_DESCRIPTIONS = [
    "LFQ Bench: Mix A vs Mix B (HYE)",
    "Custom 3-Proteome Mix",
];

export const SPECIES_CONFIG = {
    'Human': { color: '#94a3b8', expectedLog2FC: 0, label: 'Human (1:1)' }, // Slate-400
    'Yeast': { color: '#f59e0b', expectedLog2FC: 1, label: 'Yeast (2:1)' }, // Amber-500
    'E.coli': { color: '#3b82f6', expectedLog2FC: -2, label: 'E.coli (1:4)' } // Blue-500
};
