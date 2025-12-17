# A vagal liver–brain circuit for inflammation sensing and response

This repository contains R scripts/R Markdown notebooks used to generate figures for the manuscript **“A vagal liver–brain circuit for inflammation sensing and response”**.

The code is organized under `code/analysis/` (analysis workflows) and `code/figures/` (figure rendering / visualization).

---

## Repository structure

- `code/analysis/`  
  Analysis pipelines
- `code/figures/`  
  Figure-specific plotting scripts

---

## Figure-to-code mapping (Manuscript)

### Extended Data Fig. 10b–c; Extended Data Fig. 12o-q, s-u
Generated with:
- `code/visulization/Calcium imaging Heatmap_DeltaFoverF_ResponseRanking.R`

### Fig. 1j–k; Extended Data Fig. 2; Extended Data Fig. 8; Extended Data Fig. 12f, i, l; Extended Data Fig. 13o-z
Generated with:
- `code/analysis/NG_Liver_scRNAseq_DEG_GSEA_Power_CellChat.R`

### Fig. 1n
Generated with:
- `code/visulization/Enrichr analysis_visualization.Rmd`

### Fig. 2d–g; Fig. 5d–g; Extended Data Fig. 10a
Generated with:
- `code/analysis/Calcium_imaging_signal_sorting.Rmd`
