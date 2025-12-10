#!/usr/bin/env python3
"""
Proteomics Benchmark Analysis Tool
==================================

A Streamlit application for multi-category protein classification in 
3-species benchmark experiments. Implements:

1. Spectronaut data import and preprocessing
2. Protein-level quantification and fold-change calculation
3. TOST equivalence testing against expected fold-changes
4. Multi-category classification (correct, wrong magnitude, wrong direction, false positive)
5. LFQ_bout benchmark metrics (deFDR, asymmetry factor, RMSE)
6. Process capability analysis (Cpk)
7. Interactive visualizations

Author: Generated for proteomics benchmark analysis
License: MIT
"""

import streamlit as st
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import gaussian_kde
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import re
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import io

# =============================================================================
# PAGE CONFIGURATION
# =============================================================================

st.set_page_config(
    page_title="Proteomics Benchmark Analyzer",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# =============================================================================
# CUSTOM CSS
# =============================================================================

st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 1rem;
    }
    .sub-header {
        font-size: 1.2rem;
        color: #666;
        text-align: center;
        margin-bottom: 2rem;
    }
    .metric-card {
        background-color: #f0f2f6;
        border-radius: 10px;
        padding: 1rem;
        margin: 0.5rem 0;
    }
    .category-correct {
        background-color: #d4edda;
        border-left: 4px solid #28a745;
        padding: 10px;
        margin: 5px 0;
    }
    .category-magnitude {
        background-color: #fff3cd;
        border-left: 4px solid #ffc107;
        padding: 10px;
        margin: 5px 0;
    }
    .category-direction {
        background-color: #f8d7da;
        border-left: 4px solid #dc3545;
        padding: 10px;
        margin: 5px 0;
    }
    .category-fp {
        background-color: #d1ecf1;
        border-left: 4px solid #17a2b8;
        padding: 10px;
        margin: 5px 0;
    }
</style>
""", unsafe_allow_html=True)

# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class SpeciesConfig:
    """Configuration for a species in the benchmark."""
    name: str
    expected_log2fc: float
    tolerance: float  # Equivalence tolerance
    color: str

@dataclass
class ClassificationResult:
    """Result of protein classification."""
    protein_id: str
    species: str
    observed_log2fc: float
    expected_log2fc: float
    std_error: float
    p_value: float
    category: str
    tost_p_value: float
    ci_lower: float
    ci_upper: float

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def format_peptide(sequence: str) -> str:
    """Format peptide sequence, removing modifications."""
    if pd.isna(sequence):
        return None
    seq = str(sequence).strip('_')
    seq = re.sub(r'\[.*?\]', '', seq)
    return seq

def infer_species(organism: str, protein_id: str) -> str:
    """Infer species from organism column or protein ID."""
    if pd.notna(organism):
        org_lower = str(organism).lower()
        if 'homo sapiens' in org_lower or 'human' in org_lower:
            return 'human'
        elif 'saccharomyces' in org_lower or 'yeast' in org_lower:
            return 'yeast'
        elif 'escherichia' in org_lower or 'ecoli' in org_lower or 'e. coli' in org_lower:
            return 'ecoli'
    
    # Try protein ID
    if pd.notna(protein_id):
        pid = str(protein_id).upper()
        if '_HUMAN' in pid or 'HOMO' in pid:
            return 'human'
        elif '_YEAST' in pid or 'SACC' in pid:
            return 'yeast'
        elif '_ECOLI' in pid or 'ESCH' in pid:
            return 'ecoli'
    
    return 'unknown'

def calculate_protein_quantities(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate protein-level quantities using Top3 method."""
    # Group by protein and run, take top 3 peptide intensities
    protein_quant = df.groupby(['PG.ProteinGroups', 'R.FileName', 'R.Condition']).agg({
        'FG.MS2RawQuantity': lambda x: np.mean(sorted(x, reverse=True)[:3]),
        'PG.Organisms': 'first',
        'EG.IsDecoy': 'first'
    }).reset_index()
    
    protein_quant.columns = ['protein', 'run', 'condition', 'intensity', 'organism', 'is_decoy']
    return protein_quant

def calculate_fold_changes(protein_quant: pd.DataFrame, 
                           condition1: str, 
                           condition2: str) -> pd.DataFrame:
    """Calculate log2 fold changes between conditions."""
    # Pivot to wide format
    wide = protein_quant.pivot_table(
        index=['protein', 'organism', 'is_decoy'],
        columns='condition',
        values='intensity',
        aggfunc='mean'
    ).reset_index()
    
    if condition1 not in wide.columns or condition2 not in wide.columns:
        st.error(f"Conditions {condition1} and/or {condition2} not found in data")
        return None
    
    # Calculate log2 fold change
    wide['log2fc'] = np.log2(wide[condition2] / wide[condition1])
    
    # Calculate statistics (simplified - in real scenario use proper stats)
    # Here we'll use the coefficient of variation to estimate standard error
    wide['mean_intensity'] = (wide[condition1] + wide[condition2]) / 2
    wide['std_error'] = 0.3  # Placeholder - would be calculated from replicates
    
    # Simple t-test p-value (placeholder)
    wide['p_value'] = 2 * (1 - stats.norm.cdf(np.abs(wide['log2fc']) / wide['std_error']))
    
    return wide

def tost_equivalence_test(observed: float, expected: float, 
                          std_error: float, tolerance: float,
                          n: int = 3) -> Tuple[float, bool]:
    """
    Perform TOST equivalence test.
    
    Tests if observed value is equivalent to expected within tolerance bounds.
    """
    if std_error <= 0 or np.isnan(observed) or np.isnan(expected):
        return 1.0, False
    
    df = n - 1 if n > 1 else 1
    
    # Lower bound test: is observed > expected - tolerance?
    t_lower = (observed - (expected - tolerance)) / std_error
    p_lower = 1 - stats.t.cdf(t_lower, df)
    
    # Upper bound test: is observed < expected + tolerance?
    t_upper = ((expected + tolerance) - observed) / std_error
    p_upper = 1 - stats.t.cdf(t_upper, df)
    
    # TOST p-value is the maximum of the two one-sided p-values
    p_tost = max(p_lower, p_upper)
    
    # Equivalence is established if p_tost < alpha (typically 0.05)
    is_equivalent = p_tost < 0.05
    
    return p_tost, is_equivalent

def classify_protein(observed_fc: float, expected_fc: float, 
                     std_error: float, tolerance: float,
                     p_value: float, alpha: float = 0.05) -> Tuple[str, float]:
    """
    Classify protein into one of four categories:
    1. CORRECT: Correct direction AND magnitude (passes equivalence test)
    2. WRONG_MAGNITUDE: Correct direction but wrong magnitude
    3. WRONG_DIRECTION: Significant but wrong direction
    4. TRUE_NEGATIVE: Not significant when expected to be unchanged
    5. FALSE_POSITIVE: Significant when expected to be unchanged
    """
    tost_p, is_equivalent = tost_equivalence_test(observed_fc, expected_fc, std_error, tolerance)
    
    # Case 1: Background species (expected FC = 0)
    if abs(expected_fc) < 0.1:
        if p_value < alpha and abs(observed_fc) > 0.5:
            return "FALSE_POSITIVE", tost_p
        else:
            return "TRUE_NEGATIVE", tost_p
    
    # Case 2: Expected to change
    # Check direction first
    correct_direction = (observed_fc > 0 and expected_fc > 0) or (observed_fc < 0 and expected_fc < 0)
    
    if not correct_direction and p_value < alpha:
        return "WRONG_DIRECTION", tost_p
    
    # Check magnitude via equivalence
    if is_equivalent:
        return "CORRECT", tost_p
    elif correct_direction:
        return "WRONG_MAGNITUDE", tost_p
    else:
        return "NOT_SIGNIFICANT", tost_p

def calculate_defdr(classifications: List[ClassificationResult], 
                    species_expected: Dict[str, float]) -> float:
    """
    Calculate directional error FDR.
    
    deFDR = (wrong direction + false positives) / (all significant)
    """
    significant = [c for c in classifications 
                   if c.category in ['CORRECT', 'WRONG_MAGNITUDE', 'WRONG_DIRECTION', 'FALSE_POSITIVE']]
    
    if len(significant) == 0:
        return 0.0
    
    errors = [c for c in significant if c.category in ['WRONG_DIRECTION', 'FALSE_POSITIVE']]
    
    return len(errors) / len(significant)

def calculate_asymmetry_factor(log2fcs: np.ndarray) -> float:
    """
    Calculate asymmetry factor from fold-change distribution.
    
    Measured at 10% of maximum density height.
    Values 0.5-2.0 are acceptable.
    """
    if len(log2fcs) < 10:
        return np.nan
    
    log2fcs = log2fcs[~np.isnan(log2fcs)]
    if len(log2fcs) < 10:
        return np.nan
    
    try:
        kde = gaussian_kde(log2fcs)
        x_range = np.linspace(np.min(log2fcs) - 1, np.max(log2fcs) + 1, 500)
        density = kde(x_range)
        
        peak_idx = np.argmax(density)
        peak_x = x_range[peak_idx]
        threshold = 0.1 * np.max(density)
        
        # Find crossings
        left_side = density[:peak_idx]
        right_side = density[peak_idx:]
        
        left_cross = np.where(left_side < threshold)[0]
        right_cross = np.where(right_side < threshold)[0]
        
        left_x = x_range[left_cross[-1]] if len(left_cross) > 0 else x_range[0]
        right_x = x_range[peak_idx + right_cross[0]] if len(right_cross) > 0 else x_range[-1]
        
        left_dist = abs(peak_x - left_x)
        right_dist = abs(right_x - peak_x)
        
        if right_dist == 0:
            return np.inf
        
        return left_dist / right_dist
    except:
        return np.nan

def calculate_cpk(observed: np.ndarray, target: float, tolerance: float) -> float:
    """
    Calculate process capability index Cpk.
    
    Cpk = min((USL - mean)/(3σ), (mean - LSL)/(3σ))
    """
    observed = observed[~np.isnan(observed)]
    if len(observed) < 3:
        return np.nan
    
    mean = np.mean(observed)
    sigma = np.std(observed, ddof=1)
    
    if sigma == 0:
        return np.inf if abs(mean - target) < tolerance else 0
    
    usl = target + tolerance
    lsl = target - tolerance
    
    cpu = (usl - mean) / (3 * sigma)
    cpl = (mean - lsl) / (3 * sigma)
    
    return min(cpu, cpl)

# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

def plot_volcano(df: pd.DataFrame, species_colors: Dict[str, str]) -> go.Figure:
    """Create interactive volcano plot."""
    fig = go.Figure()
    
    for species, color in species_colors.items():
        mask = df['species'] == species
        subset = df[mask]
        
        fig.add_trace(go.Scatter(
            x=subset['log2fc'],
            y=-np.log10(subset['p_value'].clip(lower=1e-50)),
            mode='markers',
            name=species.capitalize(),
            marker=dict(color=color, size=8, opacity=0.7),
            text=subset['protein'],
            hovertemplate='<b>%{text}</b><br>log2FC: %{x:.2f}<br>-log10(p): %{y:.2f}<extra></extra>'
        ))
    
    # Add significance lines
    fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="gray", 
                  annotation_text="p=0.05")
    fig.add_vline(x=-0.58, line_dash="dash", line_color="gray")
    fig.add_vline(x=0.58, line_dash="dash", line_color="gray")
    
    fig.update_layout(
        title="Volcano Plot",
        xaxis_title="log2 Fold Change",
        yaxis_title="-log10(p-value)",
        height=500,
        showlegend=True
    )
    
    return fig

def plot_fc_distribution(df: pd.DataFrame, species_expected: Dict[str, float],
                         species_colors: Dict[str, str]) -> go.Figure:
    """Plot fold-change distributions with expected values."""
    fig = go.Figure()
    
    for species, expected_fc in species_expected.items():
        mask = df['species'] == species
        subset = df[mask]['log2fc'].dropna()
        
        if len(subset) > 5:
            fig.add_trace(go.Violin(
                y=subset,
                name=species.capitalize(),
                box_visible=True,
                meanline_visible=True,
                fillcolor=species_colors.get(species, 'gray'),
                opacity=0.7,
                line_color='black'
            ))
            
            # Add expected value marker
            fig.add_trace(go.Scatter(
                x=[species.capitalize()],
                y=[expected_fc],
                mode='markers',
                marker=dict(symbol='diamond', size=15, color='red', line=dict(width=2, color='black')),
                name=f'Expected ({species})',
                showlegend=False
            ))
    
    fig.update_layout(
        title="Fold-Change Distributions by Species",
        yaxis_title="log2 Fold Change",
        height=500,
        showlegend=True
    )
    
    return fig

def plot_classification_sunburst(classifications: List[ClassificationResult]) -> go.Figure:
    """Create sunburst chart of classifications."""
    # Count by species and category
    data = {}
    for c in classifications:
        key = (c.species, c.category)
        data[key] = data.get(key, 0) + 1
    
    labels = ['All Proteins']
    parents = ['']
    values = [len(classifications)]
    colors = ['#f0f0f0']
    
    category_colors = {
        'CORRECT': '#28a745',
        'WRONG_MAGNITUDE': '#ffc107',
        'WRONG_DIRECTION': '#dc3545',
        'FALSE_POSITIVE': '#17a2b8',
        'TRUE_NEGATIVE': '#6c757d',
        'NOT_SIGNIFICANT': '#adb5bd'
    }
    
    # Add species level
    species_counts = {}
    for (species, category), count in data.items():
        species_counts[species] = species_counts.get(species, 0) + count
    
    for species, count in species_counts.items():
        labels.append(species.capitalize())
        parents.append('All Proteins')
        values.append(count)
        colors.append('#e0e0e0')
    
    # Add category level
    for (species, category), count in data.items():
        labels.append(f"{category}")
        parents.append(species.capitalize())
        values.append(count)
        colors.append(category_colors.get(category, '#999999'))
    
    fig = go.Figure(go.Sunburst(
        labels=labels,
        parents=parents,
        values=values,
        marker=dict(colors=colors),
        branchvalues='total'
    ))
    
    fig.update_layout(
        title="Protein Classification Distribution",
        height=500
    )
    
    return fig

def plot_cpk_gauge(cpk: float, species: str) -> go.Figure:
    """Create gauge chart for Cpk."""
    fig = go.Figure(go.Indicator(
        mode="gauge+number+delta",
        value=cpk if not np.isnan(cpk) else 0,
        title={'text': f"Cpk - {species.capitalize()}"},
        delta={'reference': 1.33},
        gauge={
            'axis': {'range': [0, 2.5]},
            'bar': {'color': "darkblue"},
            'steps': [
                {'range': [0, 1.0], 'color': "#ffcccc"},
                {'range': [1.0, 1.33], 'color': "#ffffcc"},
                {'range': [1.33, 1.67], 'color': "#ccffcc"},
                {'range': [1.67, 2.5], 'color': "#99ff99"}
            ],
            'threshold': {
                'line': {'color': "red", 'width': 4},
                'thickness': 0.75,
                'value': 1.33
            }
        }
    ))
    
    fig.update_layout(height=300)
    return fig

def plot_expected_vs_observed(df: pd.DataFrame, species_expected: Dict[str, float],
                              species_colors: Dict[str, str]) -> go.Figure:
    """Plot observed vs expected fold changes."""
    fig = go.Figure()
    
    for species, expected_fc in species_expected.items():
        mask = df['species'] == species
        subset = df[mask]
        
        if len(subset) > 0:
            observed_median = subset['log2fc'].median()
            observed_std = subset['log2fc'].std()
            
            fig.add_trace(go.Scatter(
                x=[expected_fc],
                y=[observed_median],
                mode='markers',
                name=species.capitalize(),
                marker=dict(
                    size=20,
                    color=species_colors.get(species, 'gray'),
                    line=dict(width=2, color='black')
                ),
                error_y=dict(type='data', array=[observed_std], visible=True),
                hovertemplate=f'<b>{species}</b><br>Expected: {expected_fc:.2f}<br>Observed: %{{y:.2f}}<extra></extra>'
            ))
    
    # Add identity line
    all_fc = list(species_expected.values())
    min_fc, max_fc = min(all_fc) - 0.5, max(all_fc) + 0.5
    fig.add_trace(go.Scatter(
        x=[min_fc, max_fc],
        y=[min_fc, max_fc],
        mode='lines',
        name='Perfect Agreement',
        line=dict(dash='dash', color='gray')
    ))
    
    fig.update_layout(
        title="Expected vs Observed Fold Changes",
        xaxis_title="Expected log2 FC",
        yaxis_title="Observed log2 FC (median)",
        height=400
    )
    
    return fig

# =============================================================================
# MAIN APPLICATION
# =============================================================================

def main():
    # Header
    st.markdown('<p class="main-header">🔬 Proteomics Benchmark Analyzer</p>', unsafe_allow_html=True)
    st.markdown('<p class="sub-header">Multi-category protein classification for 3-species benchmark experiments</p>', 
                unsafe_allow_html=True)
    
    # Sidebar
    with st.sidebar:
        st.header("⚙️ Configuration")
        
        # File upload
        st.subheader("📁 Data Upload")
        uploaded_file = st.file_uploader(
            "Upload Spectronaut Report",
            type=['txt', 'tsv', 'csv'],
            help="Upload your Spectronaut export file (TSV/CSV format)"
        )
        
        st.divider()
        
        # Species configuration
        st.subheader("🧬 Species Configuration")
        
        st.markdown("**Expected log2 Fold Changes:**")
        
        human_fc = st.number_input("Human (background)", value=0.0, step=0.1,
                                   help="Expected log2FC for human proteins")
        yeast_fc = st.number_input("Yeast", value=1.0, step=0.1,
                                   help="Expected log2FC for yeast proteins")
        ecoli_fc = st.number_input("E. coli", value=-1.0, step=0.1,
                                   help="Expected log2FC for E. coli proteins")
        
        st.markdown("**Equivalence Tolerance (log2):**")
        tolerance = st.slider("Tolerance", 0.1, 1.0, 0.3, 0.05,
                             help="Acceptable deviation from expected fold change")
        
        st.divider()
        
        # Analysis parameters
        st.subheader("📊 Analysis Parameters")
        alpha = st.slider("Significance level (α)", 0.01, 0.10, 0.05, 0.01)
        min_peptides = st.slider("Min peptides per protein", 1, 5, 2)
    
    # Species configuration
    species_config = {
        'human': SpeciesConfig('human', human_fc, tolerance, '#1f77b4'),
        'yeast': SpeciesConfig('yeast', yeast_fc, tolerance, '#ff7f0e'),
        'ecoli': SpeciesConfig('ecoli', ecoli_fc, tolerance, '#2ca02c')
    }
    
    species_expected = {s: c.expected_log2fc for s, c in species_config.items()}
    species_colors = {s: c.color for s, c in species_config.items()}
    
    # Main content
    if uploaded_file is not None:
        # Load data
        try:
            # Detect separator
            content = uploaded_file.getvalue().decode('utf-8')
            sep = '\t' if '\t' in content.split('\n')[0] else ','
            uploaded_file.seek(0)
            
            df = pd.read_csv(uploaded_file, sep=sep)
            
            st.success(f"✅ Loaded {len(df):,} rows from {uploaded_file.name}")
            
            # Data preview tab
            tab1, tab2, tab3, tab4, tab5 = st.tabs([
                "📋 Data Preview", 
                "🔄 Processing", 
                "📊 Classification",
                "📈 Metrics",
                "📉 Visualizations"
            ])
            
            with tab1:
                st.header("Data Preview")
                
                col1, col2 = st.columns(2)
                
                with col1:
                    st.subheader("Available Columns")
                    st.write(list(df.columns))
                
                with col2:
                    st.subheader("Data Summary")
                    if 'R.Condition' in df.columns:
                        st.write(f"**Conditions:** {df['R.Condition'].nunique()}")
                        st.write(df['R.Condition'].value_counts())
                    if 'PG.ProteinGroups' in df.columns:
                        st.write(f"**Unique Proteins:** {df['PG.ProteinGroups'].nunique():,}")
                    if 'PG.Organisms' in df.columns:
                        st.write("**Organisms:**")
                        st.write(df['PG.Organisms'].value_counts().head(10))
                
                st.subheader("Raw Data Sample")
                st.dataframe(df.head(100), use_container_width=True)
            
            with tab2:
                st.header("Data Processing")
                
                # Check required columns
                required_cols = ['R.Condition', 'R.FileName', 'PG.ProteinGroups', 
                               'FG.MS2RawQuantity', 'PG.Organisms']
                missing_cols = [c for c in required_cols if c not in df.columns]
                
                if missing_cols:
                    st.warning(f"⚠️ Missing columns: {missing_cols}")
                    st.info("Some features may be limited without these columns.")
                
                # Process data
                st.subheader("1️⃣ Species Annotation")
                
                # Add species column
                if 'PG.Organisms' in df.columns:
                    df['species'] = df.apply(
                        lambda row: infer_species(row.get('PG.Organisms'), row.get('PG.ProteinGroups')),
                        axis=1
                    )
                else:
                    df['species'] = df['PG.ProteinGroups'].apply(lambda x: infer_species(None, x))
                
                species_dist = df['species'].value_counts()
                st.write("**Detected species distribution:**")
                
                col1, col2, col3, col4 = st.columns(4)
                col1.metric("Human", species_dist.get('human', 0))
                col2.metric("Yeast", species_dist.get('yeast', 0))
                col3.metric("E. coli", species_dist.get('ecoli', 0))
                col4.metric("Unknown", species_dist.get('unknown', 0))
                
                # Check conditions
                st.subheader("2️⃣ Condition Check")
                
                if 'R.Condition' in df.columns:
                    conditions = df['R.Condition'].unique()
                    st.write(f"**Found {len(conditions)} condition(s):** {list(conditions)}")
                    
                    if len(conditions) >= 2:
                        col1, col2 = st.columns(2)
                        with col1:
                            cond1 = st.selectbox("Control/Reference Condition", conditions, index=0)
                        with col2:
                            cond2 = st.selectbox("Treatment Condition", 
                                               [c for c in conditions if c != cond1],
                                               index=0 if len(conditions) > 1 else None)
                        
                        # Store in session state
                        st.session_state['cond1'] = cond1
                        st.session_state['cond2'] = cond2
                        st.session_state['df'] = df
                        st.session_state['species_expected'] = species_expected
                        st.session_state['species_colors'] = species_colors
                        st.session_state['tolerance'] = tolerance
                        st.session_state['alpha'] = alpha
                        
                    else:
                        st.warning("⚠️ Only one condition found. Need at least 2 for differential analysis.")
                        st.info("📌 **Demo Mode:** Simulating a second condition for demonstration purposes.")
                        
                        # Create simulated second condition
                        df_sim = df.copy()
                        df_sim['R.Condition'] = 'Simulated_Treatment'
                        
                        # Add expected fold changes based on species
                        for species, expected_fc in species_expected.items():
                            mask = df_sim['species'] == species
                            # Simulate fold change with some noise
                            noise = np.random.normal(0, 0.2, mask.sum())
                            multiplier = 2 ** (expected_fc + noise)
                            df_sim.loc[mask, 'FG.MS2RawQuantity'] *= multiplier
                        
                        df = pd.concat([df, df_sim], ignore_index=True)
                        
                        st.session_state['cond1'] = conditions[0]
                        st.session_state['cond2'] = 'Simulated_Treatment'
                        st.session_state['df'] = df
                        st.session_state['species_expected'] = species_expected
                        st.session_state['species_colors'] = species_colors
                        st.session_state['tolerance'] = tolerance
                        st.session_state['alpha'] = alpha
                        
                        st.success("✅ Simulated second condition created for demonstration.")
                else:
                    st.error("❌ R.Condition column not found!")
                
                # Filter decoys
                st.subheader("3️⃣ Decoy Filtering")
                if 'EG.IsDecoy' in df.columns:
                    decoy_count = df['EG.IsDecoy'].sum() if df['EG.IsDecoy'].dtype == bool else (df['EG.IsDecoy'] == True).sum()
                    st.write(f"**Decoy PSMs:** {decoy_count:,}")
                    df = df[df['EG.IsDecoy'] != True]
                    st.write(f"**After filtering:** {len(df):,} target PSMs")
                    st.session_state['df'] = df
                
            with tab3:
                st.header("Multi-Category Classification")
                
                if 'df' not in st.session_state:
                    st.warning("Please complete the Processing step first.")
                else:
                    df = st.session_state['df']
                    cond1 = st.session_state['cond1']
                    cond2 = st.session_state['cond2']
                    species_expected = st.session_state['species_expected']
                    tolerance = st.session_state['tolerance']
                    alpha = st.session_state['alpha']
                    
                    # Calculate protein quantities
                    with st.spinner("Calculating protein quantities..."):
                        protein_quant = calculate_protein_quantities(df)
                        fc_df = calculate_fold_changes(protein_quant, cond1, cond2)
                    
                    if fc_df is not None:
                        # Add species
                        fc_df['species'] = fc_df['organism'].apply(lambda x: infer_species(x, None))
                        
                        # Perform classification
                        st.subheader("Classification Results")
                        
                        classifications = []
                        
                        progress = st.progress(0)
                        for idx, row in fc_df.iterrows():
                            species = row['species']
                            expected_fc = species_expected.get(species, 0)
                            
                            category, tost_p = classify_protein(
                                row['log2fc'], 
                                expected_fc,
                                row['std_error'],
                                tolerance,
                                row['p_value'],
                                alpha
                            )
                            
                            classifications.append(ClassificationResult(
                                protein_id=row['protein'],
                                species=species,
                                observed_log2fc=row['log2fc'],
                                expected_log2fc=expected_fc,
                                std_error=row['std_error'],
                                p_value=row['p_value'],
                                category=category,
                                tost_p_value=tost_p,
                                ci_lower=row['log2fc'] - 1.96 * row['std_error'],
                                ci_upper=row['log2fc'] + 1.96 * row['std_error']
                            ))
                            
                            progress.progress((idx + 1) / len(fc_df))
                        
                        st.session_state['classifications'] = classifications
                        st.session_state['fc_df'] = fc_df
                        
                        # Display summary
                        st.subheader("📊 Classification Summary")
                        
                        # Count categories
                        category_counts = {}
                        for c in classifications:
                            category_counts[c.category] = category_counts.get(c.category, 0) + 1
                        
                        col1, col2, col3 = st.columns(3)
                        
                        with col1:
                            st.markdown('<div class="category-correct">', unsafe_allow_html=True)
                            st.metric("✅ Correct", category_counts.get('CORRECT', 0))
                            st.caption("Correct direction AND magnitude")
                            st.markdown('</div>', unsafe_allow_html=True)
                            
                            st.markdown('<div class="category-magnitude">', unsafe_allow_html=True)
                            st.metric("⚠️ Wrong Magnitude", category_counts.get('WRONG_MAGNITUDE', 0))
                            st.caption("Correct direction, wrong magnitude")
                            st.markdown('</div>', unsafe_allow_html=True)
                        
                        with col2:
                            st.markdown('<div class="category-direction">', unsafe_allow_html=True)
                            st.metric("❌ Wrong Direction", category_counts.get('WRONG_DIRECTION', 0))
                            st.caption("Significant but opposite direction")
                            st.markdown('</div>', unsafe_allow_html=True)
                            
                            st.markdown('<div class="category-fp">', unsafe_allow_html=True)
                            st.metric("🔵 False Positive", category_counts.get('FALSE_POSITIVE', 0))
                            st.caption("Background species showing change")
                            st.markdown('</div>', unsafe_allow_html=True)
                        
                        with col3:
                            st.metric("⚪ True Negative", category_counts.get('TRUE_NEGATIVE', 0))
                            st.metric("🔘 Not Significant", category_counts.get('NOT_SIGNIFICANT', 0))
                        
                        # Detailed results table
                        st.subheader("📋 Detailed Results")
                        
                        results_df = pd.DataFrame([{
                            'Protein': c.protein_id,
                            'Species': c.species,
                            'Observed FC': f"{c.observed_log2fc:.2f}",
                            'Expected FC': f"{c.expected_log2fc:.2f}",
                            'p-value': f"{c.p_value:.4f}",
                            'TOST p-value': f"{c.tost_p_value:.4f}",
                            'Category': c.category
                        } for c in classifications])
                        
                        # Filter by category
                        selected_categories = st.multiselect(
                            "Filter by category",
                            options=list(category_counts.keys()),
                            default=list(category_counts.keys())
                        )
                        
                        filtered_df = results_df[results_df['Category'].isin(selected_categories)]
                        st.dataframe(filtered_df, use_container_width=True)
                        
                        # Download button
                        csv = filtered_df.to_csv(index=False)
                        st.download_button(
                            "📥 Download Classification Results",
                            csv,
                            "classification_results.csv",
                            "text/csv"
                        )
            
            with tab4:
                st.header("Benchmark Metrics")
                
                if 'classifications' not in st.session_state:
                    st.warning("Please complete the Classification step first.")
                else:
                    classifications = st.session_state['classifications']
                    fc_df = st.session_state['fc_df']
                    species_expected = st.session_state['species_expected']
                    
                    # LFQ_bout metrics
                    st.subheader("📏 LFQ_bout Metrics")
                    
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        defdr = calculate_defdr(classifications, species_expected)
                        st.metric(
                            "Directional Error FDR",
                            f"{defdr:.2%}",
                            delta=f"{'✅ Good' if defdr < 0.01 else '⚠️ High'}"
                        )
                        st.caption("Target: < 1%")
                    
                    with col2:
                        asym = calculate_asymmetry_factor(fc_df['log2fc'].values)
                        st.metric(
                            "Asymmetry Factor",
                            f"{asym:.2f}" if not np.isnan(asym) else "N/A",
                            delta=f"{'✅' if 0.5 <= asym <= 2.0 else '⚠️'}" if not np.isnan(asym) else None
                        )
                        st.caption("Target: 0.5 - 2.0")
                    
                    with col3:
                        # Overall RMSE
                        rmse_vals = []
                        for c in classifications:
                            if not np.isnan(c.observed_log2fc):
                                rmse_vals.append((c.observed_log2fc - c.expected_log2fc) ** 2)
                        overall_rmse = np.sqrt(np.mean(rmse_vals)) if rmse_vals else np.nan
                        st.metric(
                            "Overall RMSE",
                            f"{overall_rmse:.3f}" if not np.isnan(overall_rmse) else "N/A"
                        )
                    
                    # Species-specific metrics
                    st.subheader("📊 Species-Specific Metrics")
                    
                    species_metrics = {}
                    for species in ['human', 'yeast', 'ecoli']:
                        species_data = fc_df[fc_df['species'] == species]['log2fc'].dropna().values
                        expected_fc = species_expected.get(species, 0)
                        
                        if len(species_data) > 2:
                            rmse = np.sqrt(np.mean((species_data - expected_fc) ** 2))
                            bias = np.mean(species_data) - expected_fc
                            cpk = calculate_cpk(species_data, expected_fc, tolerance)
                            
                            species_metrics[species] = {
                                'n': len(species_data),
                                'median': np.median(species_data),
                                'std': np.std(species_data),
                                'rmse': rmse,
                                'bias': bias,
                                'cpk': cpk
                            }
                    
                    # Display as table
                    metrics_table = pd.DataFrame(species_metrics).T
                    metrics_table.columns = ['N', 'Median FC', 'Std Dev', 'RMSE', 'Bias', 'Cpk']
                    st.dataframe(metrics_table.round(3), use_container_width=True)
                    
                    # Process Capability Gauges
                    st.subheader("🎯 Process Capability (Cpk)")
                    
                    cols = st.columns(len(species_metrics))
                    for i, (species, metrics) in enumerate(species_metrics.items()):
                        with cols[i]:
                            fig = plot_cpk_gauge(metrics['cpk'], species)
                            st.plotly_chart(fig, use_container_width=True)
                    
                    # Interpretation
                    st.subheader("📝 Interpretation Guide")
                    
                    st.markdown("""
                    | Metric | Good | Acceptable | Poor |
                    |--------|------|------------|------|
                    | **deFDR** | < 1% | 1-5% | > 5% |
                    | **Asymmetry** | 0.8-1.2 | 0.5-2.0 | < 0.5 or > 2.0 |
                    | **Cpk** | > 1.67 | 1.0-1.67 | < 1.0 |
                    | **RMSE** | < 0.3 | 0.3-0.5 | > 0.5 |
                    """)
            
            with tab5:
                st.header("Visualizations")
                
                if 'fc_df' not in st.session_state:
                    st.warning("Please complete the Classification step first.")
                else:
                    fc_df = st.session_state['fc_df']
                    classifications = st.session_state['classifications']
                    species_expected = st.session_state['species_expected']
                    species_colors = st.session_state['species_colors']
                    
                    # Volcano plot
                    st.subheader("🌋 Volcano Plot")
                    fig_volcano = plot_volcano(fc_df, species_colors)
                    st.plotly_chart(fig_volcano, use_container_width=True)
                    
                    # Distribution plot
                    st.subheader("📊 Fold-Change Distributions")
                    fig_dist = plot_fc_distribution(fc_df, species_expected, species_colors)
                    st.plotly_chart(fig_dist, use_container_width=True)
                    
                    # Expected vs Observed
                    st.subheader("🎯 Expected vs Observed")
                    fig_exp = plot_expected_vs_observed(fc_df, species_expected, species_colors)
                    st.plotly_chart(fig_exp, use_container_width=True)
                    
                    # Classification sunburst
                    st.subheader("🍩 Classification Distribution")
                    fig_sun = plot_classification_sunburst(classifications)
                    st.plotly_chart(fig_sun, use_container_width=True)
                    
                    # Add expected FC markers to volcano
                    st.subheader("📍 Fold-Change by Protein")
                    
                    # Interactive scatter with categories
                    class_df = pd.DataFrame([{
                        'protein': c.protein_id,
                        'species': c.species,
                        'observed_fc': c.observed_log2fc,
                        'expected_fc': c.expected_log2fc,
                        'category': c.category,
                        'deviation': c.observed_log2fc - c.expected_log2fc
                    } for c in classifications])
                    
                    fig_scatter = px.scatter(
                        class_df,
                        x='expected_fc',
                        y='observed_fc',
                        color='category',
                        symbol='species',
                        hover_data=['protein', 'deviation'],
                        color_discrete_map={
                            'CORRECT': '#28a745',
                            'WRONG_MAGNITUDE': '#ffc107',
                            'WRONG_DIRECTION': '#dc3545',
                            'FALSE_POSITIVE': '#17a2b8',
                            'TRUE_NEGATIVE': '#6c757d',
                            'NOT_SIGNIFICANT': '#adb5bd'
                        }
                    )
                    
                    # Add identity line
                    fig_scatter.add_trace(go.Scatter(
                        x=[-3, 3], y=[-3, 3],
                        mode='lines',
                        line=dict(dash='dash', color='gray'),
                        name='Perfect Agreement',
                        showlegend=True
                    ))
                    
                    # Add tolerance bands
                    fig_scatter.add_trace(go.Scatter(
                        x=[-3, 3], y=[-3 + tolerance, 3 + tolerance],
                        mode='lines',
                        line=dict(dash='dot', color='lightgray'),
                        name='Upper Tolerance',
                        showlegend=False
                    ))
                    fig_scatter.add_trace(go.Scatter(
                        x=[-3, 3], y=[-3 - tolerance, 3 - tolerance],
                        mode='lines',
                        line=dict(dash='dot', color='lightgray'),
                        name='Lower Tolerance',
                        showlegend=False
                    ))
                    
                    fig_scatter.update_layout(
                        title="Observed vs Expected Fold Changes by Category",
                        xaxis_title="Expected log2 FC",
                        yaxis_title="Observed log2 FC",
                        height=600
                    )
                    
                    st.plotly_chart(fig_scatter, use_container_width=True)
                    
        except Exception as e:
            st.error(f"Error loading file: {str(e)}")
            st.exception(e)
    
    else:
        # No file uploaded - show instructions
        st.info("👆 Upload a Spectronaut report file to begin analysis.")
        
        st.markdown("""
        ### 📚 How to Use This Tool
        
        1. **Upload Data**: Upload your Spectronaut export file (TSV/CSV format)
        
        2. **Configure Species**: Set expected fold changes for each species in your 3-species mix
        
        3. **Process**: The tool will:
           - Detect species annotations
           - Calculate protein-level quantities
           - Compute fold changes between conditions
        
        4. **Classify**: Each protein is classified into one of these categories:
           - ✅ **Correct**: Right direction AND magnitude (passes TOST equivalence test)
           - ⚠️ **Wrong Magnitude**: Right direction, wrong magnitude
           - ❌ **Wrong Direction**: Significant change in opposite direction
           - 🔵 **False Positive**: Background species showing unexpected change
        
        5. **Review Metrics**: Evaluate your workflow quality using:
           - **deFDR**: Directional error false discovery rate
           - **Asymmetry Factor**: Ratio compression/expansion
           - **Cpk**: Process capability index
           - **RMSE**: Root mean square error
        
        ### 📋 Required Columns
        
        Your Spectronaut export should include:
        - `R.Condition` - Experimental condition
        - `R.FileName` - Run identifier
        - `PG.ProteinGroups` - Protein group IDs
        - `PG.Organisms` - Species annotation
        - `FG.MS2RawQuantity` - Raw intensity values
        - `EG.IsDecoy` - Decoy indicator (optional but recommended)
        """)
        
        # Show sample data format
        with st.expander("📄 Sample Data Format"):
            sample_data = pd.DataFrame({
                'R.Condition': ['Control', 'Control', 'Treatment', 'Treatment'],
                'R.FileName': ['sample1.raw', 'sample2.raw', 'sample3.raw', 'sample4.raw'],
                'PG.Organisms': ['Homo sapiens', 'Homo sapiens', 'Homo sapiens', 'Homo sapiens'],
                'PG.ProteinGroups': ['P12345', 'P12345', 'P12345', 'P12345'],
                'EG.ModifiedSequence': ['_PEPTIDER_', '_PEPTIDER_', '_PEPTIDER_', '_PEPTIDER_'],
                'FG.Charge': [2, 2, 2, 2],
                'FG.MS2RawQuantity': [1000, 1100, 2000, 2200],
                'EG.IsDecoy': [False, False, False, False]
            })
            st.dataframe(sample_data)

if __name__ == "__main__":
    main()
