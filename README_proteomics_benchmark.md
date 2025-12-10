# 🔬 Proteomics Benchmark Analyzer

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Streamlit](https://img.shields.io/badge/streamlit-1.24+-red.svg)](https://streamlit.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**A Streamlit application for multi-category protein classification in 3-species benchmark experiments.**

This tool goes beyond traditional binary significance testing to classify proteins into actionable categories based on both **direction** and **magnitude** of fold changes compared to expected values.

![App Screenshot](https://via.placeholder.com/800x400?text=Proteomics+Benchmark+Analyzer)

---

## 📋 Table of Contents

1. [Overview](#overview)
2. [Features](#features)
3. [Installation](#installation)
4. [Quick Start](#quick-start)
5. [Classification System](#classification-system)
6. [Benchmark Metrics](#benchmark-metrics)
7. [Input Data Format](#input-data-format)
8. [Example Data](#example-data)
9. [Interpreting Results](#interpreting-results)
10. [Troubleshooting](#troubleshooting)
11. [Scientific Background](#scientific-background)
12. [Citation](#citation)

---

## Overview

Traditional proteomics analysis tools classify proteins as simply "significant" or "not significant" based on p-value and fold-change thresholds. This loses critical information in benchmark experiments where you **know the expected fold changes**.

For a 3-species benchmark mix (e.g., Human background, Yeast +2 log2FC, E.coli -2 log2FC):

| Observation | Traditional Tool Says | Reality |
|-------------|----------------------|---------|
| E.coli at -2.1 log2FC | ✅ Significant | ✅ Correct |
| E.coli at -4.0 log2FC | ✅ Significant | ⚠️ Wrong magnitude (over-regulated) |
| E.coli at +1.5 log2FC | ✅ Significant | ❌ Wrong direction! |
| Human at +0.8 log2FC | ✅ Significant | 🔵 False positive |

**This tool provides the complete picture.**

---

## Features

### 🎯 Multi-Category Classification
- **TOST equivalence testing** against expected fold changes
- Distinguishes correct results from magnitude and direction errors
- Identifies true false positives in background species

### 📊 Comprehensive Metrics
- **deFDR**: Directional Error False Discovery Rate
- **Asymmetry Factor**: Detects ratio compression/expansion
- **Cpk**: Process Capability Index (from industrial statistics)
- **RMSE**: Root Mean Square Error by species

### 📈 Interactive Visualizations
- Volcano plots with species coloring
- Fold-change distribution violins
- Expected vs Observed scatter plots
- Classification sunburst charts
- Process capability gauges

### 🔧 Flexible Configuration
- Customizable expected fold changes per species
- Adjustable equivalence tolerance
- Configurable significance thresholds

---

## Installation

### Prerequisites
- Python 3.8 or higher
- pip package manager

### Step 1: Clone or Download
```bash
# Create project directory
mkdir proteomics_benchmark
cd proteomics_benchmark

# Download the files (or copy from your source)
# - proteomics_benchmark_app.py
# - requirements.txt
# - run_app.sh
```

### Step 2: Install Dependencies
```bash
pip install -r requirements.txt
```

Or install manually:
```bash
pip install streamlit pandas numpy scipy plotly
```

### Step 3: Verify Installation
```bash
python -c "import streamlit; import pandas; import plotly; print('✅ All dependencies installed')"
```

---

## Quick Start

### Option 1: Using the Launcher Script
```bash
chmod +x run_app.sh
./run_app.sh
```

### Option 2: Direct Streamlit Command
```bash
streamlit run proteomics_benchmark_app.py
```

### Option 3: With Custom Port
```bash
streamlit run proteomics_benchmark_app.py --server.port 8080
```

Then open your browser to **http://localhost:8501** (or your custom port).

---

## Classification System

### The Four Categories

The app classifies each protein into one of these categories:

```
                    ┌─────────────────────────────────────────┐
                    │         PROTEIN CLASSIFICATION          │
                    └─────────────────────────────────────────┘
                                        │
                    ┌───────────────────┴───────────────────┐
                    │                                       │
            Expected to Change                    Expected Unchanged
            (Yeast, E.coli)                         (Human/Background)
                    │                                       │
            ┌───────┴───────┐                       ┌───────┴───────┐
            │               │                       │               │
      Correct          Wrong                   Not Detected      Detected
      Direction        Direction                    │               │
            │               │                       │               │
      ┌─────┴─────┐         │                TRUE_NEGATIVE   FALSE_POSITIVE
      │           │         │                    (✅)            (🔵)
   Within      Outside      │
   Tolerance   Tolerance    │
      │           │         │
   CORRECT   WRONG_MAG   WRONG_DIR
    (✅)        (⚠️)        (❌)
```

### Category Definitions

| Category | Icon | Description | Example |
|----------|------|-------------|---------|
| **CORRECT** | ✅ | Correct direction AND magnitude (passes TOST) | E.coli expected -2, observed -1.9 |
| **WRONG_MAGNITUDE** | ⚠️ | Correct direction but outside tolerance | E.coli expected -2, observed -3.5 |
| **WRONG_DIRECTION** | ❌ | Significant change in opposite direction | E.coli expected -2, observed +1.5 |
| **FALSE_POSITIVE** | 🔵 | Background species showing unexpected change | Human expected 0, observed +1.2 |
| **TRUE_NEGATIVE** | ⚪ | Background correctly shows no change | Human expected 0, observed +0.1 |
| **NOT_SIGNIFICANT** | ⚫ | Not statistically significant | Any protein with p > α |

### TOST Equivalence Testing

The app uses **Two One-Sided Tests (TOST)** to determine if an observed fold change is equivalent to the expected value:

```
H₀: |observed - expected| ≥ tolerance  (NOT equivalent)
H₁: |observed - expected| < tolerance  (EQUIVALENT)

Equivalence established when 90% CI falls entirely within:
    [expected - tolerance, expected + tolerance]
```

This is fundamentally different from traditional difference testing:
- **Difference test**: "Is it different from expected?" (failure to reject ≠ equivalence)
- **Equivalence test**: "Is it the same as expected?" (positive confirmation)

---

## Benchmark Metrics

### Directional Error FDR (deFDR)

Measures the proportion of "significant" calls that have incorrect direction:

```
deFDR = (Wrong Direction + False Positives) / (All Significant Proteins)
```

| deFDR | Interpretation | Action |
|-------|----------------|--------|
| < 1% | ✅ Excellent | Workflow validated |
| 1-5% | ⚠️ Acceptable | Review edge cases |
| > 5% | ❌ Poor | Investigate workflow |

### Asymmetry Factor

Quantifies ratio compression or expansion:

```
Asymmetry = (Left distance from peak) / (Right distance from peak)
            measured at 10% of maximum density height
```

| Asymmetry | Interpretation |
|-----------|----------------|
| 0.8 - 1.2 | ✅ Excellent - symmetric distribution |
| 0.5 - 2.0 | ⚠️ Acceptable |
| < 0.5 | ❌ Ratio compression (common in ToF) |
| > 2.0 | ❌ Ratio expansion |

### Process Capability Index (Cpk)

Borrowed from industrial Six Sigma methodology:

```
Cpk = min[(USL - μ)/(3σ), (μ - LSL)/(3σ)]

Where:
  USL = expected + tolerance (Upper Specification Limit)
  LSL = expected - tolerance (Lower Specification Limit)
  μ = observed mean
  σ = observed standard deviation
```

| Cpk | Sigma Level | Interpretation |
|-----|-------------|----------------|
| < 1.0 | < 3σ | ❌ Process not capable |
| 1.0 - 1.33 | 3-4σ | ⚠️ Marginally capable |
| 1.33 - 1.67 | 4-5σ | ✅ Capable |
| > 1.67 | > 5σ | ✅✅ Highly capable |

### RMSE (Root Mean Square Error)

```
RMSE = √[Σ(observed - expected)² / n]
```

| RMSE (log2) | Interpretation |
|-------------|----------------|
| < 0.3 | ✅ Excellent accuracy |
| 0.3 - 0.5 | ⚠️ Acceptable |
| > 0.5 | ❌ Poor accuracy |

---

## Input Data Format

### Required Columns

| Column | Description | Example |
|--------|-------------|---------|
| `R.Condition` | Experimental condition | "Control", "Treatment" |
| `R.FileName` | Run/sample identifier | "sample_01.raw" |
| `PG.ProteinGroups` | Protein group IDs | "P12345;Q67890" |
| `PG.Organisms` | Species annotation | "Homo sapiens" |
| `FG.MS2RawQuantity` | Raw MS2 intensity | 12345.67 |

### Recommended Columns

| Column | Description |
|--------|-------------|
| `EG.IsDecoy` | Decoy indicator (TRUE/FALSE) |
| `EG.Cscore` | Confidence score |
| `EG.Qvalue` | Q-value |
| `FG.Charge` | Precursor charge |

### Supported Formats
- Tab-separated values (`.tsv`, `.txt`)
- Comma-separated values (`.csv`)

---

## Example Data

### 🌟 Optimal Benchmark Data (Excellent Workflow)

This simulated data represents a **high-quality workflow** with accurate fold-change recovery:

```tsv
R.Condition	R.FileName	PG.Organisms	PG.ProteinGroups	EG.IsDecoy	EG.ModifiedSequence	EG.PEP	EG.Qvalue	EG.Cscore	FG.Charge	FG.Id	FG.MS2RawQuantity
Control	sample_ctrl_01.raw	Homo sapiens	P04406_HUMAN	FALSE	_VGVNGFGR_	0.0001	0.00001	35.5	2	_VGVNGFGR_.2	50000
Control	sample_ctrl_01.raw	Homo sapiens	P04406_HUMAN	FALSE	_LISWYDNEFGYSNR_	0.0002	0.00002	34.2	2	_LISWYDNEFGYSNR_.2	48000
Control	sample_ctrl_01.raw	Homo sapiens	P68871_HUMAN	FALSE	_VVAGVANALAHK_	0.0001	0.00001	36.1	2	_VVAGVANALAHK_.2	62000
Control	sample_ctrl_01.raw	Homo sapiens	P68871_HUMAN	FALSE	_FLASVSTVLTSK_	0.0003	0.00003	33.8	2	_FLASVSTVLTSK_.2	58000
Control	sample_ctrl_01.raw	Homo sapiens	P00558_HUMAN	FALSE	_GCITIIGGGDTATCCAK_	0.0002	0.00002	34.5	3	_GCITIIGGGDTATCCAK_.3	45000
Control	sample_ctrl_01.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_VNQIGTLSESIK_	0.0001	0.00001	35.8	2	_VNQIGTLSESIK_.2	25000
Control	sample_ctrl_01.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_TAGIQIVADDLTVTNPK_	0.0002	0.00002	34.1	2	_TAGIQIVADDLTVTNPK_.2	24000
Control	sample_ctrl_01.raw	Saccharomyces cerevisiae	P00925_YEAST	FALSE	_SIVPSGASTGVHEALEMR_	0.0001	0.00001	35.2	3	_SIVPSGASTGVHEALEMR_.3	26000
Control	sample_ctrl_01.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_VGGTSDVEVNEK_	0.0001	0.00001	36.0	2	_VGGTSDVEVNEK_.2	40000
Control	sample_ctrl_01.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_AAVEEGVVPGGGVALIR_	0.0002	0.00002	34.8	2	_AAVEEGVVPGGGVALIR_.2	38000
Control	sample_ctrl_01.raw	Escherichia coli	P0A6Y8_ECOLI	FALSE	_AIDAGFVGSVLCK_	0.0001	0.00001	35.5	2	_AIDAGFVGSVLCK_.2	42000
Control	sample_ctrl_02.raw	Homo sapiens	P04406_HUMAN	FALSE	_VGVNGFGR_	0.0001	0.00001	35.3	2	_VGVNGFGR_.2	51000
Control	sample_ctrl_02.raw	Homo sapiens	P04406_HUMAN	FALSE	_LISWYDNEFGYSNR_	0.0002	0.00002	34.0	2	_LISWYDNEFGYSNR_.2	49000
Control	sample_ctrl_02.raw	Homo sapiens	P68871_HUMAN	FALSE	_VVAGVANALAHK_	0.0001	0.00001	35.9	2	_VVAGVANALAHK_.2	63000
Control	sample_ctrl_02.raw	Homo sapiens	P68871_HUMAN	FALSE	_FLASVSTVLTSK_	0.0003	0.00003	33.5	2	_FLASVSTVLTSK_.2	59000
Control	sample_ctrl_02.raw	Homo sapiens	P00558_HUMAN	FALSE	_GCITIIGGGDTATCCAK_	0.0002	0.00002	34.3	3	_GCITIIGGGDTATCCAK_.3	46000
Control	sample_ctrl_02.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_VNQIGTLSESIK_	0.0001	0.00001	35.6	2	_VNQIGTLSESIK_.2	26000
Control	sample_ctrl_02.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_TAGIQIVADDLTVTNPK_	0.0002	0.00002	33.9	2	_TAGIQIVADDLTVTNPK_.2	25000
Control	sample_ctrl_02.raw	Saccharomyces cerevisiae	P00925_YEAST	FALSE	_SIVPSGASTGVHEALEMR_	0.0001	0.00001	35.0	3	_SIVPSGASTGVHEALEMR_.3	27000
Control	sample_ctrl_02.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_VGGTSDVEVNEK_	0.0001	0.00001	35.8	2	_VGGTSDVEVNEK_.2	41000
Control	sample_ctrl_02.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_AAVEEGVVPGGGVALIR_	0.0002	0.00002	34.6	2	_AAVEEGVVPGGGVALIR_.2	39000
Control	sample_ctrl_02.raw	Escherichia coli	P0A6Y8_ECOLI	FALSE	_AIDAGFVGSVLCK_	0.0001	0.00001	35.3	2	_AIDAGFVGSVLCK_.2	43000
Control	sample_ctrl_03.raw	Homo sapiens	P04406_HUMAN	FALSE	_VGVNGFGR_	0.0001	0.00001	35.7	2	_VGVNGFGR_.2	49000
Control	sample_ctrl_03.raw	Homo sapiens	P04406_HUMAN	FALSE	_LISWYDNEFGYSNR_	0.0002	0.00002	34.4	2	_LISWYDNEFGYSNR_.2	47000
Control	sample_ctrl_03.raw	Homo sapiens	P68871_HUMAN	FALSE	_VVAGVANALAHK_	0.0001	0.00001	36.3	2	_VVAGVANALAHK_.2	61000
Control	sample_ctrl_03.raw	Homo sapiens	P68871_HUMAN	FALSE	_FLASVSTVLTSK_	0.0003	0.00003	34.0	2	_FLASVSTVLTSK_.2	57000
Control	sample_ctrl_03.raw	Homo sapiens	P00558_HUMAN	FALSE	_GCITIIGGGDTATCCAK_	0.0002	0.00002	34.7	3	_GCITIIGGGDTATCCAK_.3	44000
Control	sample_ctrl_03.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_VNQIGTLSESIK_	0.0001	0.00001	36.0	2	_VNQIGTLSESIK_.2	24000
Control	sample_ctrl_03.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_TAGIQIVADDLTVTNPK_	0.0002	0.00002	34.3	2	_TAGIQIVADDLTVTNPK_.2	23000
Control	sample_ctrl_03.raw	Saccharomyces cerevisiae	P00925_YEAST	FALSE	_SIVPSGASTGVHEALEMR_	0.0001	0.00001	35.4	3	_SIVPSGASTGVHEALEMR_.3	25000
Control	sample_ctrl_03.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_VGGTSDVEVNEK_	0.0001	0.00001	36.2	2	_VGGTSDVEVNEK_.2	39000
Control	sample_ctrl_03.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_AAVEEGVVPGGGVALIR_	0.0002	0.00002	35.0	2	_AAVEEGVVPGGGVALIR_.2	37000
Control	sample_ctrl_03.raw	Escherichia coli	P0A6Y8_ECOLI	FALSE	_AIDAGFVGSVLCK_	0.0001	0.00001	35.7	2	_AIDAGFVGSVLCK_.2	41000
Treatment	sample_treat_01.raw	Homo sapiens	P04406_HUMAN	FALSE	_VGVNGFGR_	0.0001	0.00001	35.4	2	_VGVNGFGR_.2	50500
Treatment	sample_treat_01.raw	Homo sapiens	P04406_HUMAN	FALSE	_LISWYDNEFGYSNR_	0.0002	0.00002	34.1	2	_LISWYDNEFGYSNR_.2	48500
Treatment	sample_treat_01.raw	Homo sapiens	P68871_HUMAN	FALSE	_VVAGVANALAHK_	0.0001	0.00001	36.0	2	_VVAGVANALAHK_.2	62500
Treatment	sample_treat_01.raw	Homo sapiens	P68871_HUMAN	FALSE	_FLASVSTVLTSK_	0.0003	0.00003	33.7	2	_FLASVSTVLTSK_.2	58500
Treatment	sample_treat_01.raw	Homo sapiens	P00558_HUMAN	FALSE	_GCITIIGGGDTATCCAK_	0.0002	0.00002	34.4	3	_GCITIIGGGDTATCCAK_.3	45500
Treatment	sample_treat_01.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_VNQIGTLSESIK_	0.0001	0.00001	35.7	2	_VNQIGTLSESIK_.2	50000
Treatment	sample_treat_01.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_TAGIQIVADDLTVTNPK_	0.0002	0.00002	34.0	2	_TAGIQIVADDLTVTNPK_.2	48000
Treatment	sample_treat_01.raw	Saccharomyces cerevisiae	P00925_YEAST	FALSE	_SIVPSGASTGVHEALEMR_	0.0001	0.00001	35.1	3	_SIVPSGASTGVHEALEMR_.3	52000
Treatment	sample_treat_01.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_VGGTSDVEVNEK_	0.0001	0.00001	35.9	2	_VGGTSDVEVNEK_.2	10000
Treatment	sample_treat_01.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_AAVEEGVVPGGGVALIR_	0.0002	0.00002	34.7	2	_AAVEEGVVPGGGVALIR_.2	9500
Treatment	sample_treat_01.raw	Escherichia coli	P0A6Y8_ECOLI	FALSE	_AIDAGFVGSVLCK_	0.0001	0.00001	35.4	2	_AIDAGFVGSVLCK_.2	10500
Treatment	sample_treat_02.raw	Homo sapiens	P04406_HUMAN	FALSE	_VGVNGFGR_	0.0001	0.00001	35.2	2	_VGVNGFGR_.2	51500
Treatment	sample_treat_02.raw	Homo sapiens	P04406_HUMAN	FALSE	_LISWYDNEFGYSNR_	0.0002	0.00002	33.9	2	_LISWYDNEFGYSNR_.2	49500
Treatment	sample_treat_02.raw	Homo sapiens	P68871_HUMAN	FALSE	_VVAGVANALAHK_	0.0001	0.00001	35.8	2	_VVAGVANALAHK_.2	63500
Treatment	sample_treat_02.raw	Homo sapiens	P68871_HUMAN	FALSE	_FLASVSTVLTSK_	0.0003	0.00003	33.4	2	_FLASVSTVLTSK_.2	59500
Treatment	sample_treat_02.raw	Homo sapiens	P00558_HUMAN	FALSE	_GCITIIGGGDTATCCAK_	0.0002	0.00002	34.2	3	_GCITIIGGGDTATCCAK_.3	46500
Treatment	sample_treat_02.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_VNQIGTLSESIK_	0.0001	0.00001	35.5	2	_VNQIGTLSESIK_.2	52000
Treatment	sample_treat_02.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_TAGIQIVADDLTVTNPK_	0.0002	0.00002	33.8	2	_TAGIQIVADDLTVTNPK_.2	50000
Treatment	sample_treat_02.raw	Saccharomyces cerevisiae	P00925_YEAST	FALSE	_SIVPSGASTGVHEALEMR_	0.0001	0.00001	34.9	3	_SIVPSGASTGVHEALEMR_.3	54000
Treatment	sample_treat_02.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_VGGTSDVEVNEK_	0.0001	0.00001	35.7	2	_VGGTSDVEVNEK_.2	10200
Treatment	sample_treat_02.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_AAVEEGVVPGGGVALIR_	0.0002	0.00002	34.5	2	_AAVEEGVVPGGGVALIR_.2	9700
Treatment	sample_treat_02.raw	Escherichia coli	P0A6Y8_ECOLI	FALSE	_AIDAGFVGSVLCK_	0.0001	0.00001	35.2	2	_AIDAGFVGSVLCK_.2	10700
Treatment	sample_treat_03.raw	Homo sapiens	P04406_HUMAN	FALSE	_VGVNGFGR_	0.0001	0.00001	35.6	2	_VGVNGFGR_.2	49500
Treatment	sample_treat_03.raw	Homo sapiens	P04406_HUMAN	FALSE	_LISWYDNEFGYSNR_	0.0002	0.00002	34.3	2	_LISWYDNEFGYSNR_.2	47500
Treatment	sample_treat_03.raw	Homo sapiens	P68871_HUMAN	FALSE	_VVAGVANALAHK_	0.0001	0.00001	36.2	2	_VVAGVANALAHK_.2	61500
Treatment	sample_treat_03.raw	Homo sapiens	P68871_HUMAN	FALSE	_FLASVSTVLTSK_	0.0003	0.00003	33.9	2	_FLASVSTVLTSK_.2	57500
Treatment	sample_treat_03.raw	Homo sapiens	P00558_HUMAN	FALSE	_GCITIIGGGDTATCCAK_	0.0002	0.00002	34.6	3	_GCITIIGGGDTATCCAK_.3	44500
Treatment	sample_treat_03.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_VNQIGTLSESIK_	0.0001	0.00001	35.9	2	_VNQIGTLSESIK_.2	48000
Treatment	sample_treat_03.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_TAGIQIVADDLTVTNPK_	0.0002	0.00002	34.2	2	_TAGIQIVADDLTVTNPK_.2	46000
Treatment	sample_treat_03.raw	Saccharomyces cerevisiae	P00925_YEAST	FALSE	_SIVPSGASTGVHEALEMR_	0.0001	0.00001	35.3	3	_SIVPSGASTGVHEALEMR_.3	50000
Treatment	sample_treat_03.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_VGGTSDVEVNEK_	0.0001	0.00001	36.1	2	_VGGTSDVEVNEK_.2	9800
Treatment	sample_treat_03.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_AAVEEGVVPGGGVALIR_	0.0002	0.00002	34.9	2	_AAVEEGVVPGGGVALIR_.2	9300
Treatment	sample_treat_03.raw	Escherichia coli	P0A6Y8_ECOLI	FALSE	_AIDAGFVGSVLCK_	0.0001	0.00001	35.6	2	_AIDAGFVGSVLCK_.2	10300
```

**Expected Results with Optimal Data:**
- **Species setup**: Human = 0 (unchanged), Yeast = +1 log2FC (2x up), E.coli = -2 log2FC (4x down)
- **deFDR**: ~0% (no directional errors)
- **Asymmetry**: ~1.0 (symmetric)
- **Cpk**: >1.5 for all species
- **Classification**: All CORRECT or TRUE_NEGATIVE

---

### 💀 Bad Benchmark Data (Poor Workflow)

This simulated data represents a **problematic workflow** with ratio compression, wrong directions, and false positives:

```tsv
R.Condition	R.FileName	PG.Organisms	PG.ProteinGroups	EG.IsDecoy	EG.ModifiedSequence	EG.PEP	EG.Qvalue	EG.Cscore	FG.Charge	FG.Id	FG.MS2RawQuantity
Control	bad_ctrl_01.raw	Homo sapiens	P04406_HUMAN	FALSE	_VGVNGFGR_	0.001	0.0001	30.5	2	_VGVNGFGR_.2	50000
Control	bad_ctrl_01.raw	Homo sapiens	P04406_HUMAN	FALSE	_LISWYDNEFGYSNR_	0.002	0.0002	29.2	2	_LISWYDNEFGYSNR_.2	48000
Control	bad_ctrl_01.raw	Homo sapiens	P68871_HUMAN	FALSE	_VVAGVANALAHK_	0.001	0.0001	31.1	2	_VVAGVANALAHK_.2	62000
Control	bad_ctrl_01.raw	Homo sapiens	P68871_HUMAN	FALSE	_FLASVSTVLTSK_	0.003	0.0003	28.8	2	_FLASVSTVLTSK_.2	58000
Control	bad_ctrl_01.raw	Homo sapiens	P00558_HUMAN	FALSE	_GCITIIGGGDTATCCAK_	0.002	0.0002	29.5	3	_GCITIIGGGDTATCCAK_.3	45000
Control	bad_ctrl_01.raw	Homo sapiens	P06733_HUMAN	FALSE	_AAVPSGASTGIYEALELR_	0.001	0.0001	30.8	2	_AAVPSGASTGIYEALELR_.2	55000
Control	bad_ctrl_01.raw	Homo sapiens	P07195_HUMAN	FALSE	_VIGSGCNLDSAR_	0.002	0.0002	29.9	2	_VIGSGCNLDSAR_.2	47000
Control	bad_ctrl_01.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_VNQIGTLSESIK_	0.001	0.0001	30.8	2	_VNQIGTLSESIK_.2	25000
Control	bad_ctrl_01.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_TAGIQIVADDLTVTNPK_	0.002	0.0002	29.1	2	_TAGIQIVADDLTVTNPK_.2	24000
Control	bad_ctrl_01.raw	Saccharomyces cerevisiae	P00925_YEAST	FALSE	_SIVPSGASTGVHEALEMR_	0.001	0.0001	30.2	3	_SIVPSGASTGVHEALEMR_.3	26000
Control	bad_ctrl_01.raw	Saccharomyces cerevisiae	P00560_YEAST	FALSE	_TGQAPGFTYTDANK_	0.002	0.0002	29.4	2	_TGQAPGFTYTDANK_.2	28000
Control	bad_ctrl_01.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_VGGTSDVEVNEK_	0.001	0.0001	31.0	2	_VGGTSDVEVNEK_.2	40000
Control	bad_ctrl_01.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_AAVEEGVVPGGGVALIR_	0.002	0.0002	29.8	2	_AAVEEGVVPGGGVALIR_.2	38000
Control	bad_ctrl_01.raw	Escherichia coli	P0A6Y8_ECOLI	FALSE	_AIDAGFVGSVLCK_	0.001	0.0001	30.5	2	_AIDAGFVGSVLCK_.2	42000
Control	bad_ctrl_01.raw	Escherichia coli	P0AG67_ECOLI	FALSE	_GKPQIAVIGGGFIGLELIAR_	0.002	0.0002	29.7	3	_GKPQIAVIGGGFIGLELIAR_.3	36000
Control	bad_ctrl_02.raw	Homo sapiens	P04406_HUMAN	FALSE	_VGVNGFGR_	0.001	0.0001	30.3	2	_VGVNGFGR_.2	51000
Control	bad_ctrl_02.raw	Homo sapiens	P04406_HUMAN	FALSE	_LISWYDNEFGYSNR_	0.002	0.0002	29.0	2	_LISWYDNEFGYSNR_.2	49000
Control	bad_ctrl_02.raw	Homo sapiens	P68871_HUMAN	FALSE	_VVAGVANALAHK_	0.001	0.0001	30.9	2	_VVAGVANALAHK_.2	63000
Control	bad_ctrl_02.raw	Homo sapiens	P68871_HUMAN	FALSE	_FLASVSTVLTSK_	0.003	0.0003	28.5	2	_FLASVSTVLTSK_.2	59000
Control	bad_ctrl_02.raw	Homo sapiens	P00558_HUMAN	FALSE	_GCITIIGGGDTATCCAK_	0.002	0.0002	29.3	3	_GCITIIGGGDTATCCAK_.3	46000
Control	bad_ctrl_02.raw	Homo sapiens	P06733_HUMAN	FALSE	_AAVPSGASTGIYEALELR_	0.001	0.0001	30.6	2	_AAVPSGASTGIYEALELR_.2	56000
Control	bad_ctrl_02.raw	Homo sapiens	P07195_HUMAN	FALSE	_VIGSGCNLDSAR_	0.002	0.0002	29.7	2	_VIGSGCNLDSAR_.2	48000
Control	bad_ctrl_02.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_VNQIGTLSESIK_	0.001	0.0001	30.6	2	_VNQIGTLSESIK_.2	26000
Control	bad_ctrl_02.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_TAGIQIVADDLTVTNPK_	0.002	0.0002	28.9	2	_TAGIQIVADDLTVTNPK_.2	25000
Control	bad_ctrl_02.raw	Saccharomyces cerevisiae	P00925_YEAST	FALSE	_SIVPSGASTGVHEALEMR_	0.001	0.0001	30.0	3	_SIVPSGASTGVHEALEMR_.3	27000
Control	bad_ctrl_02.raw	Saccharomyces cerevisiae	P00560_YEAST	FALSE	_TGQAPGFTYTDANK_	0.002	0.0002	29.2	2	_TGQAPGFTYTDANK_.2	29000
Control	bad_ctrl_02.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_VGGTSDVEVNEK_	0.001	0.0001	30.8	2	_VGGTSDVEVNEK_.2	41000
Control	bad_ctrl_02.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_AAVEEGVVPGGGVALIR_	0.002	0.0002	29.6	2	_AAVEEGVVPGGGVALIR_.2	39000
Control	bad_ctrl_02.raw	Escherichia coli	P0A6Y8_ECOLI	FALSE	_AIDAGFVGSVLCK_	0.001	0.0001	30.3	2	_AIDAGFVGSVLCK_.2	43000
Control	bad_ctrl_02.raw	Escherichia coli	P0AG67_ECOLI	FALSE	_GKPQIAVIGGGFIGLELIAR_	0.002	0.0002	29.5	3	_GKPQIAVIGGGFIGLELIAR_.3	37000
Treatment	bad_treat_01.raw	Homo sapiens	P04406_HUMAN	FALSE	_VGVNGFGR_	0.001	0.0001	30.4	2	_VGVNGFGR_.2	95000
Treatment	bad_treat_01.raw	Homo sapiens	P04406_HUMAN	FALSE	_LISWYDNEFGYSNR_	0.002	0.0002	29.1	2	_LISWYDNEFGYSNR_.2	92000
Treatment	bad_treat_01.raw	Homo sapiens	P68871_HUMAN	FALSE	_VVAGVANALAHK_	0.001	0.0001	31.0	2	_VVAGVANALAHK_.2	60000
Treatment	bad_treat_01.raw	Homo sapiens	P68871_HUMAN	FALSE	_FLASVSTVLTSK_	0.003	0.0003	28.7	2	_FLASVSTVLTSK_.2	56000
Treatment	bad_treat_01.raw	Homo sapiens	P00558_HUMAN	FALSE	_GCITIIGGGDTATCCAK_	0.002	0.0002	29.4	3	_GCITIIGGGDTATCCAK_.3	88000
Treatment	bad_treat_01.raw	Homo sapiens	P06733_HUMAN	FALSE	_AAVPSGASTGIYEALELR_	0.001	0.0001	30.7	2	_AAVPSGASTGIYEALELR_.2	108000
Treatment	bad_treat_01.raw	Homo sapiens	P07195_HUMAN	FALSE	_VIGSGCNLDSAR_	0.002	0.0002	29.8	2	_VIGSGCNLDSAR_.2	90000
Treatment	bad_treat_01.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_VNQIGTLSESIK_	0.001	0.0001	30.7	2	_VNQIGTLSESIK_.2	30000
Treatment	bad_treat_01.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_TAGIQIVADDLTVTNPK_	0.002	0.0002	29.0	2	_TAGIQIVADDLTVTNPK_.2	29000
Treatment	bad_treat_01.raw	Saccharomyces cerevisiae	P00925_YEAST	FALSE	_SIVPSGASTGVHEALEMR_	0.001	0.0001	30.1	3	_SIVPSGASTGVHEALEMR_.3	15000
Treatment	bad_treat_01.raw	Saccharomyces cerevisiae	P00560_YEAST	FALSE	_TGQAPGFTYTDANK_	0.002	0.0002	29.3	2	_TGQAPGFTYTDANK_.2	70000
Treatment	bad_treat_01.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_VGGTSDVEVNEK_	0.001	0.0001	30.9	2	_VGGTSDVEVNEK_.2	60000
Treatment	bad_treat_01.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_AAVEEGVVPGGGVALIR_	0.002	0.0002	29.7	2	_AAVEEGVVPGGGVALIR_.2	57000
Treatment	bad_treat_01.raw	Escherichia coli	P0A6Y8_ECOLI	FALSE	_AIDAGFVGSVLCK_	0.001	0.0001	30.4	2	_AIDAGFVGSVLCK_.2	25000
Treatment	bad_treat_01.raw	Escherichia coli	P0AG67_ECOLI	FALSE	_GKPQIAVIGGGFIGLELIAR_	0.002	0.0002	29.6	3	_GKPQIAVIGGGFIGLELIAR_.3	80000
Treatment	bad_treat_02.raw	Homo sapiens	P04406_HUMAN	FALSE	_VGVNGFGR_	0.001	0.0001	30.2	2	_VGVNGFGR_.2	100000
Treatment	bad_treat_02.raw	Homo sapiens	P04406_HUMAN	FALSE	_LISWYDNEFGYSNR_	0.002	0.0002	28.9	2	_LISWYDNEFGYSNR_.2	98000
Treatment	bad_treat_02.raw	Homo sapiens	P68871_HUMAN	FALSE	_VVAGVANALAHK_	0.001	0.0001	30.8	2	_VVAGVANALAHK_.2	64000
Treatment	bad_treat_02.raw	Homo sapiens	P68871_HUMAN	FALSE	_FLASVSTVLTSK_	0.003	0.0003	28.4	2	_FLASVSTVLTSK_.2	60000
Treatment	bad_treat_02.raw	Homo sapiens	P00558_HUMAN	FALSE	_GCITIIGGGDTATCCAK_	0.002	0.0002	29.2	3	_GCITIIGGGDTATCCAK_.3	92000
Treatment	bad_treat_02.raw	Homo sapiens	P06733_HUMAN	FALSE	_AAVPSGASTGIYEALELR_	0.001	0.0001	30.5	2	_AAVPSGASTGIYEALELR_.2	112000
Treatment	bad_treat_02.raw	Homo sapiens	P07195_HUMAN	FALSE	_VIGSGCNLDSAR_	0.002	0.0002	29.6	2	_VIGSGCNLDSAR_.2	94000
Treatment	bad_treat_02.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_VNQIGTLSESIK_	0.001	0.0001	30.5	2	_VNQIGTLSESIK_.2	32000
Treatment	bad_treat_02.raw	Saccharomyces cerevisiae	P00924_YEAST	FALSE	_TAGIQIVADDLTVTNPK_	0.002	0.0002	28.8	2	_TAGIQIVADDLTVTNPK_.2	31000
Treatment	bad_treat_02.raw	Saccharomyces cerevisiae	P00925_YEAST	FALSE	_SIVPSGASTGVHEALEMR_	0.001	0.0001	29.9	3	_SIVPSGASTGVHEALEMR_.3	14000
Treatment	bad_treat_02.raw	Saccharomyces cerevisiae	P00560_YEAST	FALSE	_TGQAPGFTYTDANK_	0.002	0.0002	29.1	2	_TGQAPGFTYTDANK_.2	72000
Treatment	bad_treat_02.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_VGGTSDVEVNEK_	0.001	0.0001	30.7	2	_VGGTSDVEVNEK_.2	62000
Treatment	bad_treat_02.raw	Escherichia coli	P0A6F5_ECOLI	FALSE	_AAVEEGVVPGGGVALIR_	0.002	0.0002	29.5	2	_AAVEEGVVPGGGVALIR_.2	59000
Treatment	bad_treat_02.raw	Escherichia coli	P0A6Y8_ECOLI	FALSE	_AIDAGFVGSVLCK_	0.001	0.0001	30.2	2	_AIDAGFVGSVLCK_.2	24000
Treatment	bad_treat_02.raw	Escherichia coli	P0AG67_ECOLI	FALSE	_GKPQIAVIGGGFIGLELIAR_	0.002	0.0002	29.4	3	_GKPQIAVIGGGFIGLELIAR_.3	82000
```

**Problems in Bad Data:**

| Species | Expected | Actual Behavior | Problem |
|---------|----------|-----------------|---------|
| Human | 0 (unchanged) | +0.9 to +1.0 log2FC | ❌ FALSE POSITIVES |
| Yeast | +1.0 | Mixed: some +0.2, some -0.8, some +1.3 | ❌ WRONG DIRECTION & MAGNITUDE |
| E.coli | -2.0 | +0.5 to +1.0 | ❌ WRONG DIRECTION |

**Expected Results with Bad Data:**
- **deFDR**: >20% (many directional errors)
- **Asymmetry**: >2.5 or <0.4 (heavily skewed)
- **Cpk**: <0.5 for all species (not capable)
- **Classification**: Many WRONG_DIRECTION and FALSE_POSITIVE

---

## Interpreting Results

### Classification Tab

Look for:
1. **High CORRECT percentage** for spiked-in species
2. **High TRUE_NEGATIVE percentage** for background species
3. **Low/zero WRONG_DIRECTION** count
4. **Low/zero FALSE_POSITIVE** count

### Metrics Tab

| Metric | Green Flag | Yellow Flag | Red Flag |
|--------|-----------|-------------|----------|
| deFDR | < 1% | 1-5% | > 5% |
| Asymmetry | 0.8-1.2 | 0.5-2.0 | < 0.5 or > 2.0 |
| Cpk (all species) | > 1.33 | 1.0-1.33 | < 1.0 |
| RMSE | < 0.3 | 0.3-0.5 | > 0.5 |

### Visualizations Tab

1. **Volcano Plot**: Species should cluster at expected FC values
2. **Distribution Plot**: Violins should be centered on expected values (red diamonds)
3. **Expected vs Observed**: Points should fall on diagonal identity line
4. **Sunburst**: Should be dominated by green (CORRECT) and gray (TRUE_NEGATIVE)

---

## Troubleshooting

### Common Issues

#### "Only one condition found"
**Cause**: Your data only has one experimental condition.
**Solution**: The app will offer to simulate a second condition for demonstration. For real analysis, ensure you have Treatment and Control conditions.

#### "No decoy column found"
**Cause**: EG.IsDecoy column missing from export.
**Impact**: Cannot filter decoys; may affect results.
**Solution**: Re-export from Spectronaut with decoy column included.

#### "Species not detected"
**Cause**: PG.Organisms column missing or non-standard naming.
**Solution**: 
- Include PG.Organisms in export
- Use standard naming (e.g., "Homo sapiens", "Saccharomyces cerevisiae")
- Or use UniProt-style protein IDs with organism suffix (e.g., "_HUMAN")

#### High deFDR (>5%)
**Possible causes**:
1. Ratio compression in your workflow
2. Wrong expected fold changes configured
3. Poor peptide/protein ID quality
4. Batch effects between conditions

#### All proteins classified as "NOT_SIGNIFICANT"
**Possible causes**:
1. Very low actual fold changes
2. High variability between replicates
3. Standard error estimation too high
4. Significance threshold too stringent

---

## Scientific Background

### Why Multi-Category Classification?

Traditional binary classification (significant vs. not) ignores crucial information:

```
                    Traditional        Multi-Category
                    Approach           Approach
                    ─────────          ──────────────
E.coli at -2.0      Significant ✓      CORRECT ✅
E.coli at -4.0      Significant ✓      WRONG_MAGNITUDE ⚠️
E.coli at +1.0      Significant ✓      WRONG_DIRECTION ❌
Human at +0.8       Significant ✓      FALSE_POSITIVE 🔵
Human at +0.1       Not Sig.           TRUE_NEGATIVE ⚪
```

The traditional approach gives 4/5 "correct" calls.
Multi-category reveals only 2/5 are actually correct.

### TOST vs. Traditional Difference Testing

| Approach | Question Asked | Failure to Reject Means |
|----------|---------------|------------------------|
| Difference Test | "Is it different from expected?" | Could be same OR underpowered |
| Equivalence Test | "Is it equivalent to expected?" | Definitely NOT equivalent |

TOST provides positive evidence of equivalence, not just absence of difference.

### Process Capability in Proteomics

Borrowed from Six Sigma manufacturing:
- A Cpk of 1.33 means 99.99% of measurements fall within spec
- Provides intuitive quality metric for non-statisticians
- Enables benchmarking across instruments/labs

---

## Citation

If you use this tool in your research, please cite:

```bibtex
@software{proteomics_benchmark_analyzer,
  title = {Proteomics Benchmark Analyzer: Multi-category protein classification for 3-species benchmarks},
  year = {2024},
  url = {https://github.com/your-repo/proteomics-benchmark-analyzer}
}
```

### Related Publications

1. Käll, L., et al. (2019). "Integrated identification and quantification error probabilities for shotgun proteomics." *Molecular & Cellular Proteomics*.

2. Jumel, T. & Shevchenko, A. (2024). "LFQ_bout: Multispecies Benchmark Analysis for LC-MS/MS." *Journal of Proteome Research*.

3. Schubert, O.T., et al. (2015). "Building high-quality assay libraries for targeted analysis of SWATH MS data." *Nature Protocols*.

---

## License

MIT License - See LICENSE file for details.

---

## Contributing

Contributions welcome! Please submit issues and pull requests on GitHub.

---

## Changelog

### v1.0.0 (2024)
- Initial release
- Multi-category classification with TOST
- LFQ_bout metrics (deFDR, asymmetry)
- Process capability analysis (Cpk)
- Interactive Streamlit interface
