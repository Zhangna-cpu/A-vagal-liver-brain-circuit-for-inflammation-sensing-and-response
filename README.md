# A vagal liver–brain circuit for inflammation sensing and response

This repository contains R scripts/R Markdown notebooks used to generate figures for the manuscript **“A vagal liver–brain circuit for inflammation sensing and response”**.

The code is organized under `code/analysis/` (analysis workflows) and `code/figures/` (figure rendering / visualization).

> **Note:** Many scripts expect local paths (e.g., `setwd(...)` and absolute file paths). Please update paths before running.

---

## Repository structure

- `code/analysis/`  
  Analysis pipelines
- `code/figures/`  
  Figure-specific plotting scripts

---

## Figure-to-code mapping (Manuscript)

### Extended Data Fig. 10b–c
Generated with:
- `code/figures/SingleCell_CalciumHeatmap_DeltaFoverF_ResponseRanking.R`

### Fig. 1j–k; Extended Data Fig. 2a–s; Extended Data Fig. 8a–s
Generated with:
- `code/analysis/NG_Liver_scRNAseq_UPB_Labeling_Clustering_DEG_GSEA_Power_CellChat.R`

### Fig. 1N
Generated with:
- `code/figures/Enrichr analysis_revisualization.Rmd`

### Fig. 2D–G; Fig. 5D–G
Generated with:
- `code/analysis/Calcium_imaging_data_preposcessing.Rmd`
- `code/analysis/Calcium_imaging_signal_sorting.Rmd`
