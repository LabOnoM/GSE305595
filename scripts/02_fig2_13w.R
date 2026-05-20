# 02_fig2_13w.R
# Reproduces Figure 2 panels for 13w Mesenchymal cells

source("scripts/00_config.R")

message("Processing Fig 2 (13w Mesenchymal)...")
seurat_obj <- readRDS(DATA_13W_MES)
DefaultAssay(seurat_obj) <- "RNA"

# Dataset is already fully annotated with biological names and explicit factor ordering

# Invert UMAP Y-axis to match the orientation of Figure 1A in the manuscript
emb <- Embeddings(seurat_obj, "umap")
emb[, 2] <- -emb[, 2] # Flip UMAP_2
seurat_obj[["umap"]]@cell.embeddings <- emb

# Fig 2A: UMAP
message("Generating Fig 2A (UMAP)...")
p_umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, cols = cluster_colors_13w) +
  theme_publication() +
  NoAxes() +
  ggtitle("13w Mesenchymal - UMAP")
save_plot(p_umap, file.path(OUT_DIR, "Fig2/Fig2A_UMAP.png"), width = 8, height = 6)

# Fig 2B-F: FeaturePlots
features <- c("Lgr5", "Pthlh", "Prg4", "Prrx1", "Runx2", "Col1a1", "Thy1", 
              "Acta2", "Mki67", "Sox9", "Acan", "Col2a1", "Col10a1", "Ihh")

message("Generating Fig 2B-F (FeaturePlots)...")
for (feat in features) {
  if (feat %in% rownames(seurat_obj)) {
    p_feat <- FeaturePlot(seurat_obj, features = feat, reduction = "umap", order = TRUE) +
      scale_color_manuscript() +
      theme_publication()
    save_plot(p_feat, file.path(OUT_DIR, paste0("Fig2/Fig2_FeaturePlot_", feat, ".png")), width = 6, height = 5)
  } else {
    warning("Feature not found: ", feat)
  }
}

message("Figure 2 Scripts Complete.")
