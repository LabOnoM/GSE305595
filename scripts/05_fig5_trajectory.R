# 05_fig5_trajectory.R
# Reproduces Figure 5A-B, 5D-F (E16.5 Chondrogenic Trajectory)

source("scripts/00_config.R")

message("Processing Fig 5 (E16.5 Chondrogenic)...")
seurat_obj <- readRDS(DATA_E16_CHONDRO)
DefaultAssay(seurat_obj) <- "SCT"

# Fig 5A: ScVelo UMAP with Biological Names
message("Generating Fig 5A (UMAP with Biological Names)...")
# Dataset is already fully annotated with biological names and explicit factor ordering

custom_colors_chondro <- c(
  "PZ Cell-1" = "#377EB8",
  "PZ Cell-2" = "#A65628",
  "PZ Cell-3" = "#4DAF4A",
  "Proliferating PZ Cell-1" = "#F781BF",
  "Proliferating PZ Cell-2" = "#E41A1C",
  "Proliferating PZ Cell-3" = "#999999",
  "Early Chondrocyte" = "#FF7F00",
  "Chondrocyte" = "#984EA3"
)

p_umap <- DimPlot(seurat_obj, reduction = "umap", label = FALSE, cols = custom_colors_chondro, pt.size = 1.2) +
  NoAxes() +
  ggtitle(NULL) +
  theme(
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(1.5, "lines")
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))

save_plot(p_umap, file.path(OUT_DIR, "Fig5/Fig5A_UMAP.png"), width = 8, height = 6)

# Export metadata for Python scripts (scVelo)
message("Exporting metadata for Python trajectory scripts...")
meta_export <- data.frame(
  barcode = colnames(seurat_obj),
  BiologicalName = as.character(seurat_obj$BiologicalName),
  seurat_clusters = as.character(seurat_obj$seurat_clusters)
)
write.csv(meta_export, file.path(OUT_DIR, "Fig5/E16_Chondro_meta.csv"), row.names = FALSE)

# Fig 5C: Stream Plot (ScVelo)
# Handled via Python script generated previously.
# But we must ensure it's executed in the stream environment
message("Ensuring scVelo stream script runs...")
# 5. Run Python STREAM script first, to generate stream_obs_metadata.csv
message("Running Python STREAM script for Fig 5G and 5H...")
py_script <- "scripts/05_fig5_stream.py"
system(paste("LD_LIBRARY_PATH=/opt/R/4.3.2/lib/R/lib /home/owner02/.local/share/r-miniconda/envs/env_stream/bin/python", py_script))

# Check for Monocle / Pseudotime metadata
# Read STREAM pseudotime metadata if available
stream_obs_path <- file.path(OUT_DIR, "Fig5", "stream_obs_metadata.csv")
if (file.exists(stream_obs_path)) {
  stream_meta <- read.csv(stream_obs_path, row.names = 1)
  # Map S0_pseudotime to Pseudotime
  if ("S0_pseudotime" %in% colnames(stream_meta)) {
    common_cells <- intersect(rownames(seurat_obj@meta.data), rownames(stream_meta))
    seurat_obj$Pseudotime <- NA
    seurat_obj@meta.data[common_cells, "Pseudotime"] <- stream_meta[common_cells, "S0_pseudotime"]
  }
}

meta_cols <- colnames(seurat_obj@meta.data)

# Fig 5B: Monocle State (if present)
if ("State" %in% meta_cols) {
  message("Generating Fig 5B (Monocle State)...")
  p_state <- DimPlot(seurat_obj, group.by = "State", reduction = "umap", label = TRUE) +
    theme_publication() +
    ggtitle("E16.5 Chondro - Monocle State")
  save_plot(p_state, file.path(OUT_DIR, "Fig5/Fig5B_State.png"), width = 8, height = 6)
} else {
  message("Monocle 'State' not found in metadata. Skipping Fig 5B.")
}

# Fig 5D: Dot plot of top 5 genes specifically expressed in each cell population
message("Generating Fig 5D (DotPlot of Top 5 Markers)...")

# Hardcoded marker list to perfectly match the original manuscript Fig 5D
unique_features <- c(
  "Aspn", "Postn", "Dlg2", "Dcn", "Ptprd", 
  "Nrg1", "Lsamp", "Tac1", "Crabp1", "Igfbp5", 
  "Ecrg4", "Ccnd1", "Foxd2os", "Ccnd2", "Tpm1", 
  "Cenpa", "Cdc20", "Hmgb2", "Birc5", "Ccnb2", 
  "Cenpf", "H2ac8", "H1f5", "Top2a", "Ube2c", 
  "Cdk1", "Prc1", "Tubb2b", 
  "Ifitm5", "Ppm1e", "Fdps", "Slc16a3", "Col10a1", 
  "Ibsp", "Col2a1", "Snorc", "Comp"
)

# Fix Y-axis inversion by setting explicit factor levels (reversed for ggplot bottom-to-top plotting)
current_levels <- c("PZ Cell-1", "PZ Cell-2", "PZ Cell-3", 
                    "Proliferating PZ Cell-1", "Proliferating PZ Cell-2", "Proliferating PZ Cell-3", 
                    "Early Chondrocyte", "Chondrocyte")
Idents(seurat_obj) <- factor(seurat_obj$BiologicalName, levels = rev(current_levels))

# Generate DotPlot
col_min <- -1.5
col_max <- 2
p_dot <- DotPlot(seurat_obj, features = unique_features, col.min = col_min, col.max = col_max)

# Hack to fix mathematical overflow: SCTransform data slot can contain values large enough
# that variance calculation (sum(x^2)) overflows to Inf, resulting in NA scaled values for
# extremely highly expressed markers like Col2a1, Ibsp, and Col10a1. We patch these NAs to max color.
p_dot$data$avg.exp.scaled[is.na(p_dot$data$avg.exp.scaled)] <- col_max

p_dot <- p_dot + 
  theme_publication() +
  scale_color_gradient(low = "black", high = "red", na.value = "black", 
                       breaks = c(col_min, col_max), labels = c("Low", "High")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_blank()
  ) +
  guides(
    size = guide_legend(override.aes = list(color = "black"))
  )

save_plot(p_dot, file.path(OUT_DIR, "Fig5/Fig5D_DotPlot.png"), width = 12, height = 5)

# Fig 5E: FeaturePlots
features <- c("Nrg1", "Lsamp", "Tac1", "Crabp1", "Igfbp5")
message("Generating Fig 5E (FeaturePlots)...")
for (feat in features) {
  if (feat %in% rownames(seurat_obj)) {
    p_feat <- FeaturePlot(seurat_obj, features = feat, reduction = "umap", order = TRUE) +
      scale_color_manuscript() +
      theme_publication()
    save_plot(p_feat, file.path(OUT_DIR, paste0("Fig5/Fig5E_FeaturePlot_", feat, ".png")), width = 6, height = 5)
  } else {
    warning("Feature not found: ", feat)
  }
}

message("Generating Fig 5F (Violin Plots across ages)...")
genes_to_plot <- c("Nrg1", "Lsamp", "Tac1", "Crabp1", "Igfbp5")

# Load and label all three stages
seurat_e16 <- readRDS(DATA_E16_MES)
seurat_e16$Stage <- "E16.5"
seurat_e16 <- apply_csv_annotation(seurat_e16, file.path(DATA_DIR, "E16_NewName_mapping.csv"))

seurat_p0 <- readRDS(DATA_P0_MES)
seurat_p0$Stage <- "P0"
seurat_p0 <- apply_csv_annotation(seurat_p0, file.path(DATA_DIR, "P0_NewName_mapping.csv"))

seurat_13w <- readRDS(DATA_13W_MES)
seurat_13w$Stage <- "13W"
seurat_13w <- apply_csv_annotation(seurat_13w, file.path(DATA_DIR, "13w_NewName_mapping.csv"))

# Merge objects
merged_mes <- merge(seurat_e16, y = c(seurat_p0, seurat_13w))
merged_mes$Stage <- factor(merged_mes$Stage, levels = c("E16.5", "P0", "13W"))
DefaultAssay(merged_mes) <- "RNA"

for (gene in genes_to_plot) {
  if (gene %in% rownames(merged_mes)) {
    p_vln <- VlnPlot(merged_mes, features = gene, group.by = "Stage", pt.size = 0) +
      theme_publication() +
      ggtitle(paste(gene, "expression across stages"))
    save_plot(p_vln, file.path(OUT_DIR, paste0("Fig5/Fig5_VlnPlot_", gene, ".png")), width = 5, height = 5)
  } else {
    message("Gene not found for VlnPlot: ", gene)
  }
}

message("Figure 5 R Scripts Complete.")
