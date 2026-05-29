# 16p11.2-iPSC-scRNAseq
# 16p11.2 Microdeletion scRNA-seq Analysis Pipeline

This repository contains the complete computational pipeline and custom R scripts utilized for the data analysis in our manuscript:

**"The 16p11.2 microdeletion enhances gene expression variability between human IPSC derived forebrain interneuron progenitor cells in culture."**

**Preprint Link:** [medRxiv](https://www.medrxiv.org/content/10.64898/2026.05.21.26353723v1)

## Pipeline Overview

To ensure full transparency and reproducibility, our analysis is organized into a sequential, six-step pipeline. These scripts cover the workflow from raw data processing to the final integration of differential expression and variability metrics.

* **`Step_0_data_processing.R`**
  Raw 10X Genomics data ingestion, stringent quality control filtering, and initial Seurat object construction.
* **`Step_1_dot_plot_heatmap_for_Figure_1.R`**
  Cell type annotation, canonical marker visualization (DotPlots and Heatmaps), and transcriptomic correlation against in vivo references.
* **`Step_2_GSVA.R`**
  Gene Set Variation Analysis to identify robust pathway enrichment differences between genotypes.
* **`Step_3_SCENIC_analysis.R`**
  Gene Regulatory Network (GRN) inference and regulon activity scoring using the SCENIC framework.
* **`Step_4_EVA.R`**
  Expression Variation Analysis (EVA). This script evaluates multivariate dispersion to quantify the transcriptional heterogeneity driven by the microdeletion.
* **`Step_5_plot_BETA_vs_logFC.R`**
  Integration of pseudobulk differential expression (logFC) with differential variability (BETA) for final figure generation.

## System Requirements
* R (v4.0 or higher)
* Core Packages: `Seurat`, `GSVA`, `vegan` (for EVA), `ggplot2`, and `pheatmap`.
* SCENIC analysis relies on the Python-based `pySCENIC` framework alongside the R-based downstream visualizations.
