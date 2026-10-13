# GSE305595 Figure Reproduction

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20301128.svg)](https://doi.org/10.5281/zenodo.20301128)
[![GEO](https://img.shields.io/badge/GEO-GSE305595-blue.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305595)

Code to reproduce Figures 1–6 for the manuscript: **A Single-Cell Transcriptomic Atlas Identifies Hierarchical Skeletal Progenitors and Zonal Patterning in Mouse Mandibular Cartilage** (Cell Reports).

## Abstract
This repository contains the full data-processing and visualization pipeline used to generate the scRNA-seq figure panels for the manuscript. To maintain a lightweight repository size, large intermediate data objects (such as Seurat RDS files and scanpy H5AD files) are not included. Users must download the raw objects from GEO/Zenodo, place them in the correct directory structure, and execute the final baking script before running the plotting scripts.

---

## 1. Setup & Working Directory

> **IMPORTANT:**
> All scripts in this repository use relative paths. **Before running any scripts, you MUST set your working directory to the root of this repository.**

In R:
```r
setwd("/path/to/GSE305595")
```

---

## 2. Data Acquisition & Directory Structure

Download the processed datasets from the public repositories:
*   **GEO:** [GSE305595](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305595)
*   **Zenodo:** [10.5281/zenodo.20301128](https://doi.org/10.5281/zenodo.20301128)

Place the downloaded files into the `data/raw/` directory EXACTLY as shown below:

```text
GSE305595/
├── scripts/
├── output/                     # Generated figures will be saved here
└── data/
    ├── Final_Annotated_RDS/    # Will be populated by Step 1
    └── raw/
        ├── TMJ_INT_ALL.rds                 # <- Download from Zenodo
        ├── TMJ_13W_Mes.rds                 # <- Download from Zenodo
        ├── TMJ_P0_Mes.rds                  # <- Download from Zenodo
        ├── TMJ_E16_Mes.rds                 # <- Download from Zenodo
        ├── TMJ_E16_Chondro.rds             # <- Download from Zenodo
        ├── TMJ_E16_Chondro_scVelo/
        │   └── scVelo.h5ad                 # <- Download from Zenodo
        ├── TMJ_E16_Chondro_STREAM/
        │   └── STREAM_out.pkl              # <- Download from Zenodo
        └── [Mapping CSVs included in repo]
```

---

## 3. Dependencies

### R Environment
*   **R** ≥ 4.3.0
*   `Seurat`
*   `SeuratWrappers`
*   `SeuratDisk`
*   `dplyr`, `ggplot2`, `RColorBrewer`, `scales`, `viridis`, `ggrepel`, `gridExtra`, `stringr`, `openxlsx`
*   `ComplexHeatmap`, `circlize`
*   `clusterProfiler`, `org.Mm.eg.db`
*   `CellChat`

### Python Environments
*   **scVelo/Scanpy**: `scvelo`, `scanpy`, `pandas`, `numpy`, `matplotlib`
*   **STREAM**: We recommend a dedicated conda environment with `stream` installed.

---

## 4. Quick Start Execution Order

Run the scripts in the following order.

### Step 1: Bake Annotated Datasets
```r
source("scripts/99_bake_final_rds.R")
```
*(This applies biological cluster mappings and factor levels, saving the final objects to `data/Final_Annotated_RDS/`)*

### Step 2: Generate Figures
```r
source("scripts/01_fig1_atlas.R")       # Fig 1A-C
source("scripts/01_fig1_FEA.R")         # Fig 1D-F
source("scripts/02_fig2_13w.R")         # Fig 2
source("scripts/03_fig3_P0.R")          # Fig 3
source("scripts/04_fig4_E16.R")         # Fig 4
source("scripts/05_fig5_trajectory.R")  # Fig 5A-F (Calls Python scripts for 5C, 5G, 5H)
source("scripts/06_fig6e_crabp1.R")     # Fig 6E
```

*(Optional: Run `07_cellchat.R` to regenerate CellChat analysis. This is computationally intensive.)*

---

## 5. Figure Index

| Figure | Script | Output File |
| :--- | :--- | :--- |
| **Fig 1A** | `01_fig1_atlas.R` | `Fig1/Fig1A_UMAP_Stage.png` |
| **Fig 1B** | `01_fig1_atlas.R` | `Fig1/Fig1B_UMAP_Cluster.png` |
| **Fig 1C** | `01_fig1_atlas.R` | `Fig1/Fig1C_FeaturePlot_*.png` |
| **Fig 1D** | `01_fig1_FEA.R` | `Fig1/Fig1D_E16.5_GO_DotPlot.png` |
| **Fig 1E** | `01_fig1_FEA.R` | `Fig1/Fig1E_P0_GO_DotPlot.png` |
| **Fig 1F** | `01_fig1_FEA.R` | `Fig1/Fig1F_13w_GO_DotPlot.png` |
| **Fig 2** | `02_fig2_13w.R` | `Fig2/*.png` |
| **Fig 3** | `03_fig3_P0.R` | `Fig3/*.png` |
| **Fig 4** | `04_fig4_E16.R` | `Fig4/*.png` |
| **Fig 5A** | `05_fig5_trajectory.R` | `Fig5/Fig5A_UMAP.png` |
| **Fig 5B** | `05_fig5_trajectory.R` | `Fig5/Fig5B_FeaturePlot_*.png` |
| **Fig 5C** | `05_fig5_scvelo.py` | `Fig5/Fig5C_scVelo_Stream.png` |
| **Fig 5D** | `05_fig5_trajectory.R` | `Fig5/Fig5D_UMAP_Pseudotime.png` |
| **Fig 5E** | `05_fig5_trajectory.R` | `Fig5/Fig5E_FeaturePlot_*.png` |
| **Fig 5F** | `05_fig5_trajectory.R` | `Fig5/Fig5F_DotPlot.png` |
| **Fig 5G** | `05_fig5_stream.py` | `Fig5/Fig5G_STREAM_Branches.png` |
| **Fig 5H** | `05_fig5_stream.py` | `Fig5/Fig5H_STREAM_SubwayMap.png`, `Fig5H_STREAM_*.png` |
| **Fig 6E** | `06_fig6e_crabp1.R` | `Fig6/Fig6E_Crabp1_DotPlot.png` |

## Checks and deploys

This repository uses no GitHub Actions (LabOnoM decision, 2026-10-13): it holds figure-reproduction code with no automated checks, so re-run the affected scripts by hand before merging and record the result in the pull request; nothing runs automatically.

---
## License
This code is released under the **CC-BY-4.0** License.
