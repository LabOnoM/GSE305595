# 03_fig3_P0.R
# Reproduces Figure 3 panels for P0 Mesenchymal cells

source("scripts/00_config.R")

message("Processing Fig 3 (P0 Mesenchymal)...")
seurat_obj <- readRDS(DATA_P0_MES)
DefaultAssay(seurat_obj) <- "RNA"

# Dataset is already fully annotated with biological names and explicit factor ordering

custom_colors <- c(
  "SZ Cell" = "#00BFC4",
  "PZ Cell" = "#FF61CC",
  "Acta2+ PZ Cell" = "#CD9600",
  "Proliferating PZ Cell" = "#00A9FF",
  "Early Chondrocyte" = "#C77CFF",
  "Chondrocyte" = "#F8766D",
  "Hypertrophic Chondrocyte" = "#00BE67",
  "Chondrocyte-associated Osteogenic cell" = "#7CAE00"
)

# Fig 3A: UMAP
message("Generating Fig 3A (UMAP)...")

# Rotate UMAP 90 degrees counter-clockwise to match manuscript Figure 2A
emb <- Embeddings(seurat_obj, "umap")
new_x <- -emb[, 2]
new_y <- emb[, 1]
emb[, 1] <- new_x
emb[, 2] <- new_y
seurat_obj[["umap"]]@cell.embeddings <- emb

p_umap <- DimPlot(seurat_obj, reduction = "umap", label = FALSE, cols = custom_colors, pt.size = 1.2) +
  theme_publication() +
  NoAxes() +
  ggtitle(NULL) +
  theme(
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(1.5, "lines")
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))

save_plot(p_umap, file.path(OUT_DIR, "Fig3/Fig3A_UMAP.png"), width = 8, height = 6)

# Fig 3B-F: FeaturePlots
features <- c("Lgr5", "Pthlh", "Prg4", "Prrx1", "Runx2", "Col1a1", "Thy1", 
              "Acta2", "Mki67", "Sox9", "Acan", "Col2a1", "Col10a1", "Ihh")

message("Generating Fig 3B-F (FeaturePlots)...")
for (feat in features) {
  if (feat %in% rownames(seurat_obj)) {
    p_feat <- FeaturePlot(seurat_obj, features = feat, reduction = "umap", order = TRUE) +
      scale_color_manuscript() +
      theme_publication()
    save_plot(p_feat, file.path(OUT_DIR, paste0("Fig3/Fig3_FeaturePlot_", feat, ".png")), width = 6, height = 5)
  } else {
    warning("Feature not found: ", feat)
  }
}

message("Figure 3 Scripts Complete.")
