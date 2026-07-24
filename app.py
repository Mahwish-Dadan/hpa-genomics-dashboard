

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

# config
st.set_page_config(
    page_title="HPA Gene Expression Analytics",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom Academic Styling
st.markdown('''
<style>
    .main-header { font-size: 2.0rem; font-weight: 700; color: #1E3A8A; margin-bottom: 0.2rem; }
    .sub-header { font-size: 0.95rem; color: #4B5563; margin-bottom: 1.2rem; }
    .metric-card { background-color: #F8FAFC; border-radius: 8px; padding: 12px; border-left: 4px solid #2563EB; }
</style>
''', unsafe_allow_html=True)

# load clean data
@st.cache_data
def load_data():
    normal_df = pd.read_parquet("clean_normal_expression.parquet")
    cell_df = pd.read_parquet("clean_cellline_expression_.parquet")
    return normal_df, cell_df

normal_df, cell_df = load_data()


# Sidebar
st.sidebar.title("🧬Controls & Filters")

# Disease & Chromosome Filter

# Chromosome Filter 
st.sidebar.markdown("### Genomic Location")
all_chroms = ["All"] + sorted([str(c) for c in cell_df["Chromosome"].dropna().unique()])
selected_chrom = st.sidebar.selectbox("Filter by Chromosome", all_chroms, index=0)

# Disease Involvement Filter
st.sidebar.markdown("### Clinical Context")
raw_diseases = cell_df["Disease involvement"].dropna().unique()
disease_options = ["All Categories"] + sorted([str(d) for d in raw_diseases if str(d).strip() != ""])
selected_disease = st.sidebar.selectbox("Filter by Disease Involvement", disease_options, index=0)

# Apply Combined Filters
filtered_lookup_df = cell_df.copy()

if selected_chrom != "All":
    filtered_lookup_df = filtered_lookup_df[filtered_lookup_df["Chromosome"].astype(str) == selected_chrom]

if selected_disease != "All Categories":
    filtered_lookup_df = filtered_lookup_df[filtered_lookup_df["Disease involvement"] == selected_disease]

# Build dropdown options
filtered_lookup_df["Lookup_Label"] = filtered_lookup_df["Gene"] + " (" + filtered_lookup_df["Gene description"].fillna("N/A") + ")"
lookup_options = sorted(filtered_lookup_df["Lookup_Label"].unique())

st.sidebar.markdown("---")
st.sidebar.markdown("### Target Search")

if not lookup_options:
    st.sidebar.warning("No targets match the selected Chromosome / Disease filter.")
    selected_label = None
else:
    selected_label = st.sidebar.selectbox(
        "Select Target Gene / Protein Description", 
        options=lookup_options,
        index=0,
        help="Type or search by Gene Symbol, Ensembl ID, or Protein Description."
    )

#filtered_gene_df = cell_df.copy()
#if selected_chrom != "All":
#    filtered_gene_df = filtered_gene_df[filtered_gene_df["Chromosome"].astype(str) == selected_chrom]

gene_list = sorted(cell_df["Gene_Unique"].unique())
selected_gene = st.sidebar.selectbox("Select Target Gene / Ensembl ID", gene_list, index=0)

use_log = st.sidebar.toggle("Transform Scale: Log2(nTPM + 1)", value=True)
val_col = "log2_nTPM" if use_log else "nTPM"
unit_label = "log2(nTPM + 1)" if use_log else "nTPM"


# expression cutoff filter
st.sidebar.markdown("### Data Filtering")

expression_cutoff = st.sidebar.slider(
    "Detection Limit Cutoff (nTPM)",
    min_value=0.0,
    max_value=5.0,
    value=0.0,  # Default = 0.0 (Shows ALL raw HPA data without removing 0.2 nTPM)
    step=0.1,
    help="HPA defines < 1.0 nTPM as 'Not Detected'. Increase slider to filter out background transcripts."
)

# Filter Datasets
# filtered_normal = normal_df[(normal_df["Gene_Unique"] == selected_gene) & (normal_df["nTPM"] >= expression_cutoff)]
# filtered_cell = cell_df[(cell_df["Gene_Unique"] == selected_gene) & (cell_df["nTPM"] >= expression_cutoff)]

#  gene_info = cell_df[cell_df["Gene_Unique"] == selected_gene].iloc[0]

if selected_label:
    selected_gene = filtered_lookup_df[filtered_lookup_df["Lookup_Label"] == selected_label]["Gene_Unique"].values[0]

    # Filter Datasets
    filtered_normal = normal_df[(normal_df["Gene_Unique"] == selected_gene) & (normal_df["nTPM"] >= expression_cutoff)].copy()
    filtered_cell = cell_df[(cell_df["Gene_Unique"] == selected_gene) & (cell_df["nTPM"] >= expression_cutoff)].copy()

    # Dynamic log calculation safety check
    if use_log:
        if "log2_nTPM" not in filtered_normal.columns and "nTPM" in filtered_normal.columns:
            filtered_normal["log2_nTPM"] = np.log2(filtered_normal["nTPM"] + 1)
        if "log2_nTPM" not in filtered_cell.columns and "nTPM" in filtered_cell.columns:
            filtered_cell["log2_nTPM"] = np.log2(filtered_cell["nTPM"] + 1)
    gene_info = cell_df[cell_df["Gene_Unique"] == selected_gene].iloc[0]

    # Clean display values
    disease_val = gene_info.get("Disease involvement", "N/A")
    subcell_loc = gene_info.get("Subcellular main location", gene_info.get("Subcellular location", "Cytoplasm"))
    
#*While HPA data is normalized, values below 1.0 nTPM fall below standard physiological detection limits
#('Not Detected' category in HPA).
#The interactive Detection Limit Cutoff allows biologists to optionally mask trace background transcripts,
#preventing low-expression mathematical artifacts when calculating Fold Change.
    
    
#KPI

    # Header
    st.markdown(f'<div class="main-header">Target Analytics: {gene_info["Gene"]}</div>', unsafe_allow_html=True)

    st.markdown(
        f'<div class="sub-header">Protein Product: <b>{gene_info.get("Gene description", gene_info["Protein class"])}</b> | '
        f'Ensembl: <b>{gene_info.get("Ensembl", "N/A")}</b> | '
        f'Chromosome: <b>Chr {gene_info.get("Chromosome", "N/A")}</b></div>',
        unsafe_allow_html=True
    )

    #st.markdown(
    #   f'<div class="sub-header">Description: <b>{gene_info.get("Gene description", "N/A")}</b> | '
      #  f'Ensembl: <b>{gene_info.get("Ensembl", "N/A")}</b> | '
      # f'Chromosome: <b>{gene_info.get("Chromosome", "N/A")}</b></div>',
      # unsafe_allow_html=True
    #)

    # KPI Calculation
    normal_breast_val = gene_info.get("Normal_Breast_nTPM", 0.0)
    max_cell_val = filtered_cell["nTPM"].max() if not filtered_cell.empty else 0.0
    top_cell_line = filtered_cell.loc[filtered_cell["nTPM"].idxmax()]["Cell Line"] if max_cell_val > 0 else "N/A"
    avg_log2_fc = filtered_cell["log2_FC_vs_Breast"].mean() if not filtered_cell.empty else 0.0

    m1, m2, m3, m4 = st.columns(4)
    with m1:
        st.metric("Normal Breast Baseline", f"{normal_breast_val:.2f} nTPM")
    with m2:
        st.metric("Peak Cancer Line", f"{top_cell_line}", f"{max_cell_val:.2f} nTPM")
    with m3:
        st.metric("Avg Log2 Fold Change", f"{avg_log2_fc:+.2f}", delta=f"{gene_info.get('Dysregulation_Status', 'N/A')}")
    with m4:
        st.metric("Tissue Specificity", f"{gene_info.get('RNA tissue specificity', 'N/A')}")

    st.divider()



    # MULTI-TAB DASHBOARD
    tab1, tab2, tab3, tab4 = st.tabs([
        "Primary Target Profiler",
        "Organ Similarity Matcher",
        "Cell Line Model Selector",
        "Multi-Gene Co-Expression"
    ])

    # TAB 1: PRIMARY TARGET PROFILER
    with tab1:

        # Metadata Overview Panel inside Tab 1
        with st.expander("View Comprehensive Biological & Protein Metadata ℹ️ ", expanded=True):
            col_a, col_b, col_c = st.columns(3)
            with col_a:
                st.write(f"• **Gene Symbol:** {gene_info['Gene']}")
                st.write(f"• **Ensembl ID:** {gene_info.get('Ensembl', 'N/A')}")
                st.write(f"• **Protein ID:** {gene_info.get('Uniprot', 'N/A')}") 
                st.write(f"• **Chromosomal Locus:** Chr {gene_info.get('Chromosome', 'N/A')}")
            with col_b:
                st.write(f"• **Protein Name:** {gene_info.get('Gene description', 'N/A')}")
                st.write(f"• **Subcellular Location:** {gene_info.get('Subcellular main location', 'N/A')}")
                st.write(f"• **HPA Specificity Class:** {gene_info.get('RNA tissue specificity', 'N/A')}")
            with col_c:
                st.write(f"• **Dysregulation Status:** {gene_info.get('Dysregulation_Status', 'N/A')}")
                st.write(f"• **Normal Breast Level:** {normal_breast_val:.2f} nTPM")
                st.write(f"• **Cancer Max Level:** {max_cell_val:.2f} nTPM ({top_cell_line})")
                st.write(f"• **Disease Involvement:** {{gene_info.get('Disease involvement', 'N/A')}") 

            
            st.write(f"• **Protein Class:** {gene_info.get('Protein class', 'N/A')}")
            
            st.write(f"• **Evidence:** {gene_info.get('Evidence', 'N/A')}")


        st.markdown("---")

        st.subheader("Healthy Breast Baseline vs. Cancer Cell Line Distributions")

        c1, c2 = st.columns([1, 2])

        with c1:
            fig_norm = px.bar(
                filtered_normal, x="Tissue", y=val_col, color="Tissue",
                title=f"Healthy Normal Tissue Profiles ({unit_label})",
                labels={val_col: unit_label}
            )
            fig_norm.update_layout(showlegend=False, height=450)
            st.plotly_chart(fig_norm, use_container_width=True)

        with c2:
            fig_cell = px.box(
                filtered_cell, y=val_col, points="all", hover_data=["Cell Line", "nTPM", "log2_FC_vs_Breast"],
                title=f"Distribution Across All Cancer Cell Lines ({unit_label})",
                labels={val_col: unit_label}
            )
            fig_cell.update_layout(height=450)
            st.plotly_chart(fig_cell, use_container_width=True)


    # TAB 2: ORGAN MIMICRY & ECTOPIC MATCHER
    with tab2:
        st.subheader("Organ Similarity Matcher")
        st.caption("Identify whether cancer cell line expression deviates from breast tissue and matches another healthy organ's profile.")

        if not filtered_cell.empty and not filtered_normal.empty:
            selected_line = st.selectbox("Select Cell Line to Analyze Organ Proximity", sorted(filtered_cell["Cell Line"].unique()))
            cell_val = filtered_cell[filtered_cell["Cell Line"] == selected_line]["nTPM"].values[0]

            # Calculate Absolute Expression Distance across Normal Organs
            organ_match = filtered_normal.copy()
            organ_match["Abs_Difference"] = (organ_match["nTPM"] - cell_val).abs()
            organ_match = organ_match.sort_values(by="Abs_Difference", ascending=True)

            closest_organ = organ_match.iloc[0]["Tissue"]
            closest_diff = organ_match.iloc[0]["Abs_Difference"]

            st.info(f"**Lineage Match Result:** Cell line **{selected_line}** (nTPM = {cell_val:.2f}) most closely matches **Healthy {closest_organ} Tissue** (nTPM = {organ_match.iloc[0]['nTPM']:.2f}, Δ = {closest_diff:.2f}).")

            fig_match = px.bar(
                organ_match, x="Tissue", y="nTPM", color="Abs_Difference",
                color_continuous_scale="Viridis_r",
                title=f"Organ Expression Distance for {selected_line} (Lower Delta = Closer Match)",
                labels={"nTPM": "Normal Tissue nTPM", "Abs_Difference": "|Δ nTPM|"}
            )
            # Add Horizontal Line for Cell Line Value
            fig_match.add_hline(y=cell_val, line_dash="dash", line_color="red", annotation_text=f"{selected_line} Level ({cell_val:.1f} nTPM)")
            fig_match.update_layout(height=450)
            st.plotly_chart(fig_match, use_container_width=True)


    # TAB 3: CELL LINE MODEL SELECTOR
    with tab3:
        st.subheader("Custom In Vitro Model Selector")

        total_available_lines = len(cell_df["Cell Line"].unique())

        ctrl_col1, ctrl_col2 = st.columns(2)
        with ctrl_col1:
            sort_order = st.radio("Display Ranking", ["Top Expressors (Highest nTPM)", "Bottom Expressors (Lowest/Negative Controls)"], horizontal=True)
        with ctrl_col2:
            num_display = st.slider("Number of Cell Lines to Display", min_value=5, max_value=total_available_lines, value=min(20, total_available_lines), step=1)

        is_ascending = True if "Bottom" in sort_order else False
        ranked_df = filtered_cell.sort_values(by=val_col, ascending=is_ascending).head(num_display)

        fig_rank = px.bar(
            ranked_df, x="Cell Line", y=val_col, color="log2_FC_vs_Breast",
            color_continuous_scale="RdBu_r",
            title=f"{sort_order.split(' ')[0]} {num_display} Cell Lines for {gene_info['Gene']}",
            labels={val_col: unit_label, "log2_FC_vs_Breast": "Log2 FC vs Breast"}
        )
        fig_rank.update_layout(height=480, xaxis_tickangle=-45)
        st.plotly_chart(fig_rank, use_container_width=True)


    # TAB 4: MULTI-GENE CO-EXPRESSION
    with tab4:
        st.subheader("Multi-Gene Expression Heatmap")

        col_g, col_cl = st.columns(2)
        with col_g:
            multi_selected_genes = st.multiselect("Select Genes to Compare", gene_list, default=gene_list[:min(5, len(gene_list))])
        with col_cl:
            all_lines = sorted(cell_df["Cell Line"].unique())
            default_subset = all_lines[:10]
            select_all = st.checkbox("Select All 62 Cell Lines", value=False)

            if select_all:
                selected_cell_lines = st.multiselect("2. Select Cell Lines", all_lines, default=all_lines)
            else:
                selected_cell_lines = st.multiselect("2. Select Cell Lines (Default: Top 10)", all_lines, default=default_subset)

        if len(multi_selected_genes) > 1 and len(selected_cell_lines) > 0:
            sub_df = cell_df[(cell_df["Gene_Unique"].isin(multi_selected_genes)) & (cell_df["Cell Line"].isin(selected_cell_lines))]

            pivot_df = sub_df.pivot_table(
                index="Gene", columns="Cell Line", values=val_col, aggfunc="mean"
            ).fillna(0)

            st.markdown("### Expression Intensity Heatmap")
            fig_heat = px.imshow(
                pivot_df, labels=dict(x="Cell Line", y="Gene", color=unit_label),
                title=f"Co-Expression Matrix Across {len(selected_cell_lines)} Cell Lines",
                color_continuous_scale="Magma", aspect="auto"
            )
            fig_heat.update_layout(height=480)
            st.plotly_chart(fig_heat, use_container_width=True)

            # Pearson Correlation Matrix Section
            st.markdown("---")
            show_corr = st.checkbox("Enable Gene-Gene Pearson Correlation Analysis", value=False)

            if show_corr:
                if len(multi_selected_genes) >= 3:
                    # st.markdown("#### Gene-Gene Expression Correlation Matrix (Pearson r)")

                    corr_matrix = pivot_df.T.corr()

                    c_fig, c_exp = st.columns([1.2, 1])

                    with c_fig:

                        fig_corr = px.imshow(
                            corr_matrix, text_auto=".2f",
                            color_continuous_scale="RdBu_r", zmin=-1, zmax=1,
                            title="Gene Co-Regulation Correlation (Pearson r)"
                        )
                        fig_corr.update_layout(height=400)
                        st.plotly_chart(fig_corr, use_container_width=True)

                    with c_exp:
                            st.markdown("### Correlation Interpretation Guide 💡 ")
                            st.info("""
                            **Pearson Coefficient ($r$) Ranges:**
                            * **$r \ge +0.70$ (Strong Co-Regulation):** Genes tend to be activated or suppressed together across cell lines, suggesting shared signaling pathways or transcriptional controls.
                            * **$-0.30 < r < +0.30$ (Independent):** Expression profiles operate independently.
                            * **$r \le -0.70$ (Inverse Co-Regulation):** High expression of one gene coincides with suppression of the other (e.g., mutually exclusive oncogenic drivers).
                            """)

            else:
                st.warning("Please select at least 2 genes and 1 cell line to display the matrix.")


    # export option

    st.sidebar.markdown("---")
    st.sidebar.download_button(
        "Download Filtered Data (CSV)",
        data=filtered_cell.to_csv(index=False).encode('utf-8'),
        file_name=f"{gene_info['Gene']}_expression_data.csv",
        mime="text/csv"
    )

    
