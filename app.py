import streamlit as st
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import plotly.express as px
import yaml, json, io, zipfile
from pathlib import Path

# ---- Utility: File Handlers and Config ----
@st.cache_data
def read_matrix(file):
    return pd.read_csv(file, sep='\t', low_memory=False)

def find_matrix_files(folder):
    files = list(Path(folder).glob("*.tsv"))
    pg = next((f for f in files if "pg.matrix" in f.name or "protein" in f.name.lower()), None)
    pr = next((f for f in files if "pr.matrix" in f.name or "precursor" in f.name.lower()), None)
    return pg, pr

def save_config(config, file_type='yaml'):
    if file_type == 'yaml':
        return yaml.dump(config).encode()
    else:
        return json.dumps(config, indent=2).encode()

def load_config(uploaded, file_type):
    if file_type == 'yaml':
        return yaml.safe_load(uploaded)
    else:
        return json.load(uploaded)

def zip_results(dict_of_files):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for fname, fdata in dict_of_files.items():
            zf.writestr(fname, fdata)
    buf.seek(0)
    return buf

# ---- Synthetic Data Demo Mode ----
def generate_synthetic_data(species_params, n_proteins=2000, n_reps=3, noise=0.13):
    np.random.seed(42)
    species_ids = {"HUMAN":0,"YEAST":1,"ECOLI":2,"CELEGANS":3}
    base = []
    for sp, fc in species_params.items():
        for i in range(n_proteins//len(species_params)):
            row = {
                "Protein.Group": f"P{sp}{i+1}",
                "Protein.Names": f"{sp}_dummy{i+1}",
            }
            a_vals = np.random.normal(20, noise, n_reps)
            b_vals = a_vals + fc
            for j in range(n_reps):
                row[f"A{j+1:02d}"] = 2**a_vals[j]
                row[f"B{j+1:02d}"] = 2**b_vals[j]
            base.append(row)
    return pd.DataFrame(base)

# ---- Preprocessing / Filtering ----
def assign_species(row, species_suffixes):
    matches = [sp for sp,suf in species_suffixes.items() if suf in row['Protein.Names']]
    return matches[0] if matches else 'MIXED'

def filter_and_annotate(df, condA, condB, species_map, thresholds):
    df = df.copy()
    df['Species'] = df.apply(assign_species, axis=1, species_suffixes=species_map)
    df = df[~df['Species'].eq('MIXED')]
    quant_mask = (
        df[condA].notna().sum(axis=1) >= thresholds['min_valid'] and
        df[condB].notna().sum(axis=1) >= thresholds['min_valid']
    )
    df = df[quant_mask]
    df['CV_A'] = df[condA].std(axis=1) / df[condA].mean(axis=1)
    df['CV_B'] = df[condB].std(axis=1) / df[condB].mean(axis=1)
    cv_mask = (df['CV_A'] <= thresholds['cv_cutoff']) & (df['CV_B'] <= thresholds['cv_cutoff'])
    return df[cv_mask]

def add_log2(df, condA, condB):
    df = df.copy()
    for col in condA+condB:
        if not col.startswith("log2_"):
            if (df[col] > 0).all():
                df[f'log2_{col}'] = np.log2(df[col])
    return df

def diffexp(df, condA, condB, fc_thresh, alpha):
    l2a = df[[f'log2_{c}' for c in condA]].values
    l2b = df[[f'log2_{c}' for c in condB]].values
    t_pvals = [ttest_ind(a, b, equal_var=False, nan_policy='omit')[1] for a,b in zip(l2a, l2b)]
    l2fc = l2a.mean(axis=1) - l2b.mean(axis=1)
    df['log2FC'] = l2fc
    df['pval'] = t_pvals
    df['signif'] = (np.abs(df['log2FC']) >= fc_thresh) & (df['pval'] <= alpha)
    df['regulation'] = np.where(df['log2FC'] >= fc_thresh, 'UP',
                        np.where(df['log2FC'] <= -fc_thresh, 'DOWN', 'NS'))
    return df

# ---- Summary/QC/Benchmark Metrics ----
def compute_metrics(df, species_params):
    stats = {}
    deFDR = (df['Species'] != df['regulation']).mean() * 100 if "regulation" in df else np.nan
    mean_cv = np.nanmean(np.r_[df['CV_A'], df['CV_B']]) * 100
    asym = np.abs(df['log2FC']).median() / (np.abs(df['log2FC']).mean() + 1e-6)
    stats['deFDR'] = deFDR
    stats['mean_CV'] = mean_cv
    stats['asymmetry'] = asym
    # TP / FP / FN, sensitivity, specificity – demo only
    expected_up = [sp for sp,fc in species_params.items() if fc>0]
    tp = sum((df['regulation']=='UP') & df['Species'].isin(expected_up))
    fp = sum((df['regulation']=='UP') & ~df['Species'].isin(expected_up))
    fn = sum((df['regulation']=='DOWN') & df['Species'].isin(expected_up))
    sens = tp/(tp+fn+1e-6)
    spec = tp/(tp+fp+1e-6)
    stats.update(dict(TP=tp, FP=fp, FN=fn, Sensitivity=sens, Specificity=spec))
    return stats

def pass_fail(stats):
    return {
        "deFDR": stats["deFDR"] <= 1.0,
        "mean CV": stats["mean_CV"] <= 5,
        "asymmetry": 0.5 <= stats["asymmetry"] <= 2.0,
    }

# ---- Main Streamlit App ----
st.set_page_config('Proteomics Benchmark', layout='wide')
st.title("Multi-Species Proteomics Workflow Benchmarking App")
st.markdown(
"""**End-to-end benchmarking and QC of DIA-NN bottom-up proteomics workflows with multi-species validation.**
*Implements the protocol described by Jumel, Tobias, Shevchenko (2024, JPR) and related works.*
"""
)
with st.expander("Show CITATION"):
    st.markdown(
        "Jumel, T. & Shevchenko, A. Multispecies Benchmark Analysis for LC-MSMS Validation and Performance Evaluation in Bottom-Up Proteomics. J. Proteome Res. 2024; see README for other references."
    )

tab1, tab2, tab3 = st.tabs(["Data", "Parameters", "Results"])

with tab1:
    st.header("Step 1: Data Ingestion")
    demo = st.checkbox("Demo Mode: Use synthetic test data")
    folder = None
    if not demo:
        folder = st.text_input("Path to DIA-NN results folder")
        if folder:
            pg, pr = find_matrix_files(folder)
            if pg:
                st.success(f"Protein Group matrix found: {pg.name}")
            if pr:
                st.info(f"Precursor matrix found: {pr.name}")
    else:
        st.info("Synthetic example data enabled.")
    if demo or (folder and pg):
        if demo:
            species_params = st.session_state.get("species_params", {"HUMAN":0,"YEAST":1,"ECOLI":-2})
            df = generate_synthetic_data(species_params)
        else:
            df = read_matrix(pg)
        st.dataframe(df.head(20))

with tab2:
    st.header("Step 2: Assign Samples & Set Parameters")
    all_cols = df.columns.tolist()
    condA = st.multiselect("Group A sample columns", [c for c in all_cols if c.startswith("A")])
    condB = st.multiselect("Group B sample columns", [c for c in all_cols if c.startswith("B")])
    min_valid = st.slider("Min valid values/group", 1, len(condA), 2)
    cv_cut = st.slider("Max CV %", 1, 50, 20)
    fc_thresh = st.slider("Log2 FC threshold", 0, 4, 1)
    alpha = st.number_input("Significance level (alpha)", 0.001, 0.2, 0.01)
    species_params = st.text_area(
        "Expected log2FC by species (JSON)", 
        '{"HUMAN":0,"YEAST":1,"ECOLI":-2,"CELEGANS":-1}'
    )
    if st.button("Save Params/Config"):
        st.download_button(
            "Download config",
            save_config(dict(condA=condA, condB=condB, min_valid=min_valid, cv_cut=cv_cut, fc_thresh=fc_thresh, alpha=alpha, species_params=json.loads(species_params))),
            "config.yaml"
        )

with tab3:
    st.header("Step 3: Analyze and Benchmark")
    # Only proceed if columns are assigned etc.
    if condA and condB:
        # Preprocess
        thresholds = dict(min_valid=min_valid, cv_cutoff=cv_cut/100)
        species_map = {"HUMAN":"HUMAN","YEAST":"YEAST","ECOLI":"ECOLI","CELEGANS":"CELEGANS"}
        dff = filter_and_annotate(df, condA, condB, species_map, thresholds)
        dff = add_log2(dff, condA, condB)
        dff = diffexp(dff, condA, condB, fc_thresh, alpha)
        stats = compute_metrics(dff, json.loads(species_params))
        pf = pass_fail(stats)
        # Dashboard
        st.subheader("Summary Dashboard")
        st.markdown(f"- deFDR: {stats['deFDR']:.2f} {'✅' if pf['deFDR'] else '❌'}")
        st.markdown(f"- mean CV: {stats['mean_CV']:.2f} {'✅' if pf['mean CV'] else '❌'}")
        st.markdown(f"- Asymmetry: {stats['asymmetry']:.3f} {'✅' if pf['asymmetry'] else '❌'}")
        st.markdown(f"- TP: {stats['TP']} | FP: {stats['FP']} | FN: {stats['FN']} | Sensitivity: {stats['Sensitivity']:.2f} | Specificity: {stats['Specificity']:.2f}")

        # Plots
        fig1 = px.scatter(dff, x='log2FC', y='pval', color='Species', title="log2FC vs p-value (Volcano)")
        fig2 = px.histogram(dff, x='CV_A', nbins=30, title="CV Distribution A")
        fig3 = px.density_contour(dff, x='log2FC', y='pval', title="log2FC Density vs p-value")
        st.plotly_chart(fig1, use_container_width=True)
        st.plotly_chart(fig2, use_container_width=True)
        st.plotly_chart(fig3, use_container_width=True)
        st.dataframe(dff.head(40))

        export_btn = st.button("Download Full Results")
        if export_btn:
            out_bytes = zip_results({
                "quant_results.tsv": dff.to_csv(sep="\t", index=False),
                "dashboard_stats.json": json.dumps(stats, indent=2),
                "config.yaml": save_config(dict(condA=condA, condB=condB, min_valid=min_valid, cv_cut=cv_cut, fc_thresh=fc_thresh, alpha=alpha, species_params=species_params)),
            })
            st.download_button("Results.zip", out_bytes, "BenchmarkResults.zip")
