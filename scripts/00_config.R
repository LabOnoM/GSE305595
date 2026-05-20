# 00_config.R
# Shared configuration and utility functions for TosaCellReports Figure Reproduction

suppressMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(SeuratDisk)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
  library(viridis)
  library(ggrepel)
  library(gridExtra)
  library(ComplexHeatmap)
  library(circlize)
  library(stringr)
  library(openxlsx)
})

# Define Paths
# --- PATHS ---
PROJECT_ROOT <- getwd()
DATA_DIR       <- file.path(PROJECT_ROOT, "data", "raw")
BAKED_DATA_DIR <- file.path(PROJECT_ROOT, "data", "Final_Annotated_RDS")
OUT_DIR        <- file.path(PROJECT_ROOT, "output")

# Final Baked Datasets (Annotations and Factors Pre-Applied)
DATA_P0_MES      <- file.path(BAKED_DATA_DIR, "P0_Mesenchymal_Annotated.rds")
DATA_E16_MES     <- file.path(BAKED_DATA_DIR, "E16.5_Mesenchymal_Annotated.rds")
DATA_E16_CHONDRO <- file.path(BAKED_DATA_DIR, "E16.5_Chondro_Annotated.rds")
DATA_13W_MES     <- file.path(BAKED_DATA_DIR, "13w_Mesenchymal_Annotated.rds")

# Unbaked/Original Atlas Datasets
DATA_INT_ALL <- file.path(DATA_DIR, "DataInt.TMJ_all.rds")
DATA_INT_MES <- file.path(DATA_DIR, "TMJ_all_IntData_Mes.rds")

# Ensure output directories exist
fig_dirs <- paste0("Fig", 1:6)
sapply(file.path(OUT_DIR, fig_dirs), dir.create, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "CellChat"), showWarnings = FALSE, recursive = TRUE)

# Standard Aesthetics
# Defines standard ggplot2 theme for reproducible, clean figures
theme_publication <- function(base_size=14, base_family="sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

# Standard Color Palette for clusters (consistent with distinct color palettes)
cluster_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
  "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", 
  "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"
)

cluster_colors_13w <- c(
  "SZ Cell" = "#00C19F",
  "SZ-associated Osteogenic cell" = "#00BA38",
  "Thy1+ PZ Cell" = "#619CFF",
  "Acta2+ PZ Cell" = "#D39200",
  "Proliferating PZ cell" = "#00B9E3",
  "Early Chondrocyte" = "#DB72FB",
  "Chondrocyte" = "#F8766D",
  "Hypertrophic Chondrocyte" = "#93AA00",
  "Chondrocyte-associated Osteogenic cell" = "#FF61C3"
)

# Feature plot manuscript aesthetic
scale_color_manuscript <- function() {
  scale_colour_gradientn(colours = c("#C0C0C080", "#00000080","#0000FF80","#FF000080","#FFA50080"))
}

# --- BIOLOGICAL NAME MAPPINGS ---
# Deduced from marker gene profile analysis and Combined Manuscript Nov 22.pdf

annotation_map_E16_Chondro <- c(
  "0" = "PZ Cell-1",
  "1" = "Early Chondrocyte",
  "2" = "PZ Cell-3",
  "3" = "Proliferating PZ Cell-2",
  "4" = "Chondrocyte",
  "5" = "PZ Cell-2",
  "6" = "Proliferating PZ Cell-1",
  "7" = "Proliferating PZ Cell-3"
)

annotation_map_E16 <- c(
  "0" = "Acta2+ PZ Cell",
  "1" = "Proliferating PZ Cell",
  "2" = "Mesenchymal Cells",
  "3" = "Early Chondrocyte",
  "4" = "Chondrocyte",
  "5" = "Gli1+ PZ Cell",
  "6" = "Pthlh+ SZ/Disc Progenitor",
  "7" = "Mesenchymal Cells",
  "8" = "Proliferating SZ/Disc Progenitor",
  "9" = "Lgr5+ SZ/Disc Progenitor",
  "10" = "Mesenchymal Cells"
)

annotation_map_P0 <- c(
  "0" = "Early Chondrocyte",
  "1" = "Acta2+ PZ Cell",
  "2" = "PZ Cell",
  "3" = "Chondrocyte",
  "4" = "PZ Cell",
  "5" = "Proliferating PZ Cell",
  "6" = "SZ Cell",
  "7" = "Hypertrophic Chondrocyte",
  "8" = "Chondrocyte-associated Osteogenic cell",
  "9" = "Proliferating PZ Cell"
)

annotation_map_13w_biological <- c(
  "Lgr5+" = "SZ Cell",
  "Dlx5+" = "SZ-associated Osteogenic cell",
  "Prg4+" = "Thy1+ PZ Cell",
  "Acta2+" = "Acta2+ PZ Cell",
  "Mki67+" = "Proliferating PZ cell",
  "Sox9+" = "Early Chondrocyte",
  "Acan+" = "Chondrocyte",
  "Col10a1+" = "Hypertrophic Chondrocyte",
  "Sp7+" = "Chondrocyte-associated Osteogenic cell"
)

# Utility Function: Apply Mappings
apply_annotation_map <- function(seurat_obj, map_dict) {
  # Default to seurat_clusters if available
  if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    Idents(seurat_obj) <- seurat_obj$seurat_clusters
  }
  seurat_obj <- RenameIdents(seurat_obj, map_dict)
  seurat_obj$BiologicalName <- Idents(seurat_obj)
  return(seurat_obj)
}

# Utility Function: Apply CSV Annotation Mappings (For Fig 2 & 3)
apply_csv_annotation <- function(seurat_obj, csv_path) {
  if (file.exists(csv_path)) {
    mapping <- read.csv(csv_path, stringsAsFactors=FALSE)
    # Remove duplicate barcodes caused by multiple HTML tooltip entries
    mapping <- mapping[!duplicated(mapping$Cell.ID), ]
    
    # The cell IDs in the CSV match the barcodes in the Seurat object
    rownames(mapping) <- mapping$Cell.ID
    common_cells <- intersect(rownames(seurat_obj@meta.data), mapping$Cell.ID)
    
    seurat_obj$BiologicalName <- NA
    seurat_obj@meta.data[common_cells, "BiologicalName"] <- mapping[common_cells, "NewName"]
    
    # Filter out NA cells (cells without a new name annotation)
    seurat_obj <- subset(seurat_obj, subset = !is.na(BiologicalName))
    
    Idents(seurat_obj) <- seurat_obj$BiologicalName
  } else {
    warning(paste("CSV not found:", csv_path))
  }
  return(seurat_obj)
}

# Utility Function: Save Plot
save_plot <- function(plot, filename, width, height, dpi = 300) {
  ggsave(filename, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  message("Saved: ", filename)
}

message("00_config.R loaded successfully. Paths, annotations, and standard aesthetics defined.")
