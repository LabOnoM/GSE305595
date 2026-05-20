# 01_fig1_atlas.R
# Reproduces Figure 1A-C panels for the all-stage integrated atlas

source("scripts/00_config.R")

message("Processing Fig 1 (Atlas - All Stages)...")
seurat_obj <- readRDS(DATA_INT_ALL)
DefaultAssay(seurat_obj) <- "RNA"

# Fig 1A: UMAP by stage
message("Generating Fig 1A (UMAP by Stage)...")
p_umap_stage <- DimPlot(seurat_obj, group.by = "orig.ident", reduction = "umap", label = FALSE) +
  theme_publication() +
  ggtitle("TMJ Atlas - By Stage")
save_plot(p_umap_stage, file.path(OUT_DIR, "Fig1/Fig1A_UMAP_Stage.png"), width = 8, height = 6)

# Fig 1B: UMAP by cluster
message("Generating Fig 1B (UMAP by Cluster)...")
p_umap_cluster <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE) +
  theme_publication() +
  ggtitle("TMJ Atlas - By Cluster")
save_plot(p_umap_cluster, file.path(OUT_DIR, "Fig1/Fig1B_UMAP_Cluster.png"), width = 8, height = 6)

# Fig 1C: FeaturePlots
features <- c("Prrx1", "Sox9")

message("Generating Fig 1C (FeaturePlots)...")
for (feat in features) {
  if (feat %in% rownames(seurat_obj)) {
    p_feat <- FeaturePlot(seurat_obj, features = feat, reduction = "umap", order = TRUE) +
      scale_color_viridis(option = "C") +
      theme_publication()
    save_plot(p_feat, file.path(OUT_DIR, paste0("Fig1/Fig1C_FeaturePlot_", feat, ".png")), width = 6, height = 5)
  } else {
    warning("Feature not found: ", feat)
  }
}

message("Figure 1 Atlas Scripts Complete.")
