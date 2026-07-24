# Human Protein Atlas (HPA) Gene Expression Analytics 🧬 

### **Live Web Application**
> **Click here to explore the live dashboard:**  
> **(https://hpa-genomics-dashboard.streamlit.app/)**

---

## Overview

The **HPA Gene Expression Analytics Dashboard** is an interactive, bioinformatic tool designed to explore and visualize transcriptomic profiles across healthy normal tissues and cancer cell lines derived from the **Human Protein Atlas (HPA)** data models.

Built with performance, memory optimization, and clinical relevance in mind, this application enables researchers to interrogate gene expression levels, evaluate tumor vs. healthy tissue fold changes, and identify top *in vitro* cancer cell line models for experimental targeting.

---

## Key Features

* **Targeted Search & Genomic Filters:** Search genes by **Gene Symbol**, **Ensembl ID**, or **Protein Name/Description**. Filter dynamically by **Chromosome** and **Disease Involvement**.
* **Scalable Log Transformation:** Dynamic toggling between raw normalized Transcripts Per Million (**nTPM**) and logarithmic transformation $\log_2(\text{nTPM} for high dynamic range expression datasets.
* **Primary Target Profiler:** Compare healthy tissue baselines against cancer cell line distributions with interactive boxplots and hover tooltips.
* **Cell Line Model Selector:** Dynamic ranking of top or bottom expressors (useful for identifying positive/negative controls for *in vitro* validation assays).
* **Organ Similarity Matcher:** Quantify expression delta ($\Delta\text{nTPM}$) to identify which healthy tissue profile a selected cancer line most closely aligns with.
* **Multi-Gene Co-Expression & Heatmaps:** Compare multiple targets across cell line panels with integrated **Pearson Correlation Analysis** ($r$) for gene co-regulation discovery.
* **Data Export:** Instant CSV export of custom filtered sub-datasets.

---

## Performance Optimization & Stability

Large multi-tab bioinformatic apps running on cloud-hosted environments (like Streamlit Cloud) often suffer from Out-Of-Memory (OOM) `SIGKILL` server crashes. This dashboard implements memory-conscious architecture:

1. **Callback Garbage Collection:** Uses custom `on_change` event callbacks to trigger `gc.collect()`, immediately purging unused UI objects and Plotly figure buffers upon user input.
2. **Modular Tab Architecture:** Segregates heavy data pivoting and matrix rendering to minimize peak RAM consumption.
3. **Parquet Storage Engine:** Data is stored in compressed `.parquet` format for fast zero-copy reads and low disk I/O overhead.

