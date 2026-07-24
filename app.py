import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import gc

# Config
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

@st.cache_data
def load_data():
    cell_df = pd.read_parquet("clean_cellline_expression_.parquet")

    normal_df = pd.read_parquet("clean_normal_expression.parquet")

    # Pre-calculate log2_nTPM ONCE in cache
    if "log2_nTPM" not in normal_df.columns and "nTPM" in normal_df.columns:
        normal_df["log2_nTPM"] = np.log2(normal_df["nTPM"] + 1)
        
    if "log2_nTPM" not in cell_df.columns and "nTPM" in cell_df.columns:
        cell_df["log2_nTPM"] = np.log2(cell_df["nTPM"] + 1)
        
    return normal_df, cell_df

# Callback function to wipe dead variables from RAM when a toggle changes
def free_ram_callback():
    gc.collect()

normal_df, cell_df = load_data()
try:
    
    # Sidebar
    st.sidebar.title("🧬 Controls & Filters")
    
    # Chromosome Filter 
    st.sidebar.markdown("### Genomic Location")
    all_chroms = ["All"] + sorted([str(c) for c in cell_df["Chromosome"].dropna().unique()])
    selected_chrom = st.sidebar.selectbox("Filter by Chromosome", all_chroms, index=0, on_change=free_ram_callback)
    
    # Disease Involvement Filter
    st.sidebar.markdown("### Clinical Context")
    raw_diseases = cell_df["Disease involvement"].dropna().unique()
    disease_options = ["All Categories"] + sorted([str(d) for d in raw_diseases if str(d).strip() != ""])
    selected_disease = st.sidebar.selectbox("Filter by Disease Involvement", disease_options, index=0, on_change=free_ram_callback)
    
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
    
    gene_list = sorted(cell_df["Gene_Unique"].unique())
    
    use_log = st.sidebar.toggle("Transform Scale: Log2(nTPM + 1)", value=True, on_change=free_ram_callback)
    val_col = "log2_nTPM" if use_log else "nTPM"
    unit_label = "Log₂ (nTPM + 1)" if use_log else "nTPM"
    
    # Expression cutoff filter
    st.sidebar.markdown("### Data Filtering")
    expression_cutoff = st.sidebar.slider(
        "Detection Limit Cutoff (nTPM)",
        min_value=0.0,
        max_value=5.0,
        value=0.0,
        step=0.1,
        help="HPA defines < 1.0 nTPM as 'Not Detected'. Increase slider to filter out background transcripts."
    )
    
    if selected_label:
        selected_gene = filtered_lookup_df[filtered_lookup_df["Lookup_Label"] == selected_label]["Gene_Unique"].values[0]
    
        # Filter Datasets
        filtered_normal = normal_df[(normal_df["Gene_Unique"] == selected_gene) & (normal_df["nTPM"] >= expression_cutoff)].copy()
        filtered_cell = cell_df[(cell_df["Gene_Unique"] == selected_gene) & (cell_df["nTPM"] >= expression_cutoff)].copy()
            
        gene_info = cell_df[cell_df["Gene_Unique"] == selected_gene].iloc[0]
    
        # Header
        st.markdown(f'<div class="main-header">Target Analytics: {gene_info["Gene"]}</div>', unsafe_allow_html=True)
    
        st.markdown(
            f'<div class="sub-header">Protein Product: <b>{gene_info.get("Gene description", "N/A")}</b> | '
            f'Protein ID: <b>{gene_info.get("Uniprot", "N/A")}</b> | '
            f'Ensembl: <b>{gene_info.get("Ensembl", "N/A")}</b> | ' 
            f'Chromosome: <b>Chr {gene_info.get("Chromosome", "N/A")}</b></div>',
            unsafe_allow_html=True
        )
    
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
        #tab1, tab2, tab3, tab4 = st.tabs([
        #   "Primary Target Profiler",
        #    "Organ Similarity Matcher",
        #    "Cell Line Model Selector",
        #    "Multi-Gene Co-Expression"
        #])

        selected_tab = st.radio(
            "Navigation", 
            ["Primary Target Profiler", "Organ Similarity Matcher", "Cell Line Model Selector", "Multi-Gene Co-Expression"], 
            horizontal=True,
            label_visibility="collapsed"
        )
    
        # TAB 1: PRIMARY TARGET PROFILER
        # with tab1:
        if selected_tab == "Primary Target Profiler":
            with st.expander("View Metadata ℹ️", expanded=True):
                col_a, col_b, col_c = st.columns(3)
                with col_a:
                    st.write(f"• **Gene Symbol:** {gene_info['Gene']}")
                    st.write(f"• **Ensembl ID:** {gene_info.get('Ensembl', 'N/A')}")
                    st.write(f"• **Chromosomal Locus:** Chr {gene_info.get('Chromosome', 'N/A')}")
                    st.write(f"• **Subcellular Location:** {gene_info.get('Subcellular main location', 'N/A')}")
                with col_b:
                    st.write(f"• **Protein ID:** {gene_info.get('Uniprot', 'N/A')}") 
                    st.write(f"• **Protein Name:** {gene_info.get('Gene description', 'N/A')}")
                    st.write(f"• **HPA Specificity Class:** {gene_info.get('RNA tissue specificity', 'N/A')}")
                    st.write(f"• **Disease Involvement:** {gene_info.get('Disease involvement', 'N/A')}") 
                with col_c:
                    st.write(f"• **Normal Breast Level:** {normal_breast_val:.2f} nTPM")
                    st.write(f"• **Cancer Max Level:** {max_cell_val:.2f} nTPM ({top_cell_line})")
                    st.write(f"• **Dysregulation Status:** {gene_info.get('Dysregulation_Status', 'N/A')}")
    
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
    
        # TAB 2: ORGAN SIMILARITY MATCHER
        #with tab2:
        elif selected_tab == "Organ Similarity Matcher":
            st.subheader("Organ Similarity Matcher")
            st.caption("Identify whether cancer cell line expression deviates from breast tissue and matches another healthy organ's profile.")
    
            if not filtered_cell.empty and not filtered_normal.empty:
                selected_line = st.selectbox("Select Cell Line to Analyze Organ Proximity", sorted(filtered_cell["Cell Line"].unique()))
                cell_val = filtered_cell[filtered_cell["Cell Line"] == selected_line]["nTPM"].values[0]
    
                organ_match = filtered_normal.copy()
                organ_match["Abs_Difference"] = (organ_match["nTPM"] - cell_val).abs()
                organ_match = organ_match.sort_values(by="Abs_Difference", ascending=True)
    
                closest_organ = organ_match.iloc[0]["Tissue"]
                closest_diff = organ_match.iloc[0]["Abs_Difference"]

                delta_val = closest_diff
                
                # Categorize confidence level based on biological thresholds
                if delta_val <= 0.1:
                    confidence_label = "Strong Match / High Similarity"
                    explanation = "The expression level is nearly identical to this normal tissue baseline."
                elif delta_val <= 0.5:
                    confidence_label = "Moderate Match / Plausible Variant"
                    explanation = "The expression is relatively close, but represents a minor variance."
                else:
                    confidence_label = "Weak / Relative Match Only"
                    explanation = "This is simply the closest tissue available in the panel, but the expression difference is substantial."
                
                # Display in Streamlit
                st.metric(
                    label=f"Closest Matching Tissue: {closest_organ}",
                    value=f"{organ_match.iloc[0]['nTPM']:.2f} nTPM",
                    delta=f"Δ {delta_val:.2f} difference",
                    delta_color="inverse"
                )
                
                st.info(f"**Similarity Level:** {confidence_label}\n\n*{explanation}*")
    
                #st.info(f"**Lineage Match Result:** Cell line **{selected_line}** (nTPM = {cell_val:.2f}) most closely matches **Healthy {closest_organ} Tissue** (nTPM = {organ_match.iloc[0]['nTPM']:.2f}, Δ = {closest_diff:.2f}).")
    
                fig_match = px.bar(
                    organ_match, x="Tissue", y="nTPM", color="Abs_Difference",
                    color_continuous_scale="Viridis_r",
                    title=f"Organ Expression Distance for {selected_line} (Lower Delta = Closer Match)",
                    labels={"nTPM": "Normal Tissue nTPM", "Abs_Difference": "|Δ nTPM|"}
                )
                fig_match.add_hline(y=cell_val, line_dash="dash", line_color="red", annotation_text=f"{selected_line} Level ({cell_val:.1f} nTPM)")
                fig_match.update_layout(height=450)
                st.plotly_chart(fig_match, use_container_width=True)
    
        # TAB 3: CELL LINE MODEL SELECTOR
        # with tab3:
        elif selected_tab == "Cell Line Model Selector":
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
        # with tab4:
        elif selected_tab == "Multi-Gene Co-Expression":
            st.subheader("Multi-Gene Expression Heatmap")
    
            col_g, col_cl = st.columns(2)
            with col_g:
                multi_selected_genes = st.multiselect("Select Genes to Compare", gene_list, default=gene_list[:min(5, len(gene_list))])
            with col_cl:
                all_lines = sorted(cell_df["Cell Line"].unique())
                default_subset = all_lines[:min(10, len(all_lines))]
                select_all = st.checkbox("Select All Cell Lines", value=False, on_change=free_ram_callback)
    
                if select_all:
                    selected_cell_lines = st.multiselect("Select Cell Lines", all_lines, default=all_lines)
                else:
                    selected_cell_lines = st.multiselect("Select Cell Lines", all_lines, default=default_subset)
    
            if len(multi_selected_genes) >= 1 and len(selected_cell_lines) >= 1:
                sub_df = cell_df[
                    (cell_df["Gene_Unique"].isin(multi_selected_genes)) & 
                    (cell_df["Cell Line"].isin(selected_cell_lines))
                ].copy()
    
                if not sub_df.empty:
                    # Pivot using Gene_Unique to prevent duplicate row index issues
                    pivot_df = sub_df.pivot_table(
                        index="Gene_Unique", columns="Cell Line", values=val_col, aggfunc="mean"
                    ).fillna(0)

                    # Map Gene_Unique back to clean Gene Symbol for display
                    gene_map = cell_df.set_index("Gene_Unique")["Gene"].to_dict()
                    pivot_df.index = [gene_map.get(idx, idx) for idx in pivot_df.index]
    
                    st.markdown("### Expression Intensity Heatmap")
                    fig_heat = px.imshow(
                        pivot_df, labels=dict(x="Cell Line", y="Gene Symbol", color=unit_label),
                        title=f"Co-Expression Matrix Across {len(selected_cell_lines)} Cell Lines",
                        color_continuous_scale="Magma", aspect="auto"
                    )
                    fig_heat.update_layout(height=480)
                    st.plotly_chart(fig_heat, use_container_width=True)
    
                    # Pearson Correlation Matrix Section
                    st.markdown("---")
                    show_corr = st.checkbox("Enable Gene-Gene Pearson Correlation Analysis", value=False)
    
                    if show_corr:
                        if len(multi_selected_genes) >= 2:
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
                                st.markdown("### Correlation Interpretation Guide 💡")
                                st.info("""
                                **Pearson Coefficient ($r$) Ranges:**
                                * **$r \\ge +0.70$ (Strong Co-Regulation):** Genes tend to be activated or suppressed together across cell lines.
                                * **$-0.30 < r < +0.30$ (Independent):** Expression profiles operate independently.
                                * **$r \\le -0.70$ (Inverse Co-Regulation):** High expression of one gene coincides with suppression of the other.
                                """)
                        else:
                            st.warning("Please select at least 2 genes to display the correlation matrix.")
            else:
                st.warning("Please select at least 1 gene and 1 cell line to display the matrix.")
    
        # Export option
        st.sidebar.markdown("---")
        st.sidebar.download_button(
            "Download Filtered Data (CSV)",
            data=filtered_cell.to_csv(index=False).encode('utf-8'),
            file_name=f"{gene_info['Gene']}_expression_data.csv",
            mime="text/csv"
        )

        

except Exception as e:
    st.error("⚠️ An unexpected error occurred while running the analytics dashboard.")
    with st.expander("Show error details (for developers)"):
        st.exception(e)
        
    if st.button("Reset & Reload App", type="primary"):
        st.cache_data.clear()
        st.rerun()


# Place at the absolute end of your script
gc.collect()
