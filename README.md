# A vagal liver–brain circuit for inflammation sensing and response

This repository contains R code and small demo datasets used to reproduce **selected analyses / figures** from the manuscript:

**“A vagal liver–brain circuit for inflammation sensing and response”**

The repository is organized as research scripts (not an R package). For Nature-style code review, the repo includes:
- **Version-pinned environment** via `renv.lock`
- **Small demo inputs** under `data/demo/`
- **Reference (expected) outputs** under `expected_output/demo/`

---

## System requirements

### Tested software / OS
- **R**: 4.3.1
- **Platform tested**: x86_64-w64-mingw32/x64 (Windows 64-bit)

> The full list of R package dependencies (with versions) is recorded in `renv.lock`.

### Hardware
- No non-standard hardware required.

### Typical install time
- `renv::restore()` may take **~10–60 minutes** on first run depending on network speed and whether packages are installed from binaries vs source.

---

## Installation / Demo

### 1) Download / clone this repository
Use GitHub “Download ZIP” or `git clone`, then open an R session at the repository root.

### 2) Restore the exact package versions
In R (from the repo root):

```r
install.packages("renv", repos = "https://cloud.r-project.org")
renv::restore()
```

### 3) Demo inputs, scripts, and expected outputs
#### (A) Calcium imaging demo
- **Input demo data**: `data/demo/calcium/All the responsve neurons 2.xlsx`
- **Script**: `code/analysis/Calcium_imaging_signal_sorting.R`
- **What it does**: Processes calcium imaging ΔF/F traces, performs response typing/classification, and generates line plots summarizing response patterns.
- **Expected outputs**: See `expected_output/demo/calcium/`.

#### (B) Enrichr visualization demo
- **Input demo data**: `data/demo/enrichr/HMDB.xlsx`
- **Script**: `code/visulization/Enrichr_analysis_visualization.R`
- **What it does**: Generates Enrichr-related visualization(s) (e.g., lollipop-style plots) based on the HMDB demo table.
- **Expected outputs**: See `expected_output/demo/enrichr/`.

#### Expected runtime
> The expected run time for the demo on a “normal” desktop computer is around **10–15 minutes** (depending on package installation state and system performance).
---

## Repository structure

- `code/analysis/`  
  - Analysis pipelines
  - `Calcium_imaging_signal_sorting.R`
  - `NG_Liver_scRNAseq_DEG_GSEA_Power_CellChat.R`
- `code/visualization/`  
  - Figure-specific plotting scripts
  - `Enrichr_analysis_visualization.R`
  - `Calcium imaging Heatmap_DeltaFoverF_ResponseRanking.R`
- `data/demo/`
  - Small demo input datasets
  - `calcium`
  - `enrichr`
- `expected_output/demo/`
  - Reference outputs for the demo datasets
  - `calcium`
  - `enrichr`
- `renv.lock`  
  - Full dependency list with pinned versions for reproducibility

> Note: some “full analysis” workflows may require larger raw datasets not included here. The demo is intended to allow reviewers to run representative steps end-to-end.

---

## Figure-to-code mapping (Manuscript)

**Extended Data Fig. 10b–c; Extended Data Fig. 12o-q, s-u**  
Generated with:
- `code/visulization/Calcium_imaging_Heatmap_DeltaFoverF_ResponseRanking.R`

**Fig. 1j–k; Extended Data Fig. 2; Extended Data Fig. 8; Extended Data Fig. 12f, i, l; Extended Data Fig. 13o-z**  
Generated with:
- `code/analysis/NG_Liver_scRNAseq_DEG_GSEA_Power_CellChat.R`

**Fig. 1n**  
Generated with:
- `code/visulization/Enrichr_analysis_visualization.Rmd`

**Fig. 2d–g; Fig. 5d–g; Extended Data Fig. 10a**  
Generated with:
- `code/analysis/Calcium_imaging_signal_sorting.Rmd`
