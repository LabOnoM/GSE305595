# 04_fig4_E16.R
# Reproduces Figure 4 panels for E16.5 Mesenchymal cells

source("scripts/00_config.R")

message("Processing Fig 4 (E16.5 Mesenchymal)...")
seurat_obj <- readRDS(DATA_E16_MES)
DefaultAssay(seurat_obj) <- "RNA"

# Dataset is already fully annotated with biological names and explicit factor ordering

custom_colors <- c(
  "Mesenchymal Cells" = "#93AA00",
  "Proliferating SZ/Disc Progenitor" = "#00B9E3",
  "Lgr5+ SZ/Disc Progenitor" = "#00BA38",
  "Pthlh+ SZ/Disc Progenitor" = "#619CFF",
  "Gli1+ PZ Cell" = "#FF61C3",
  "Acta2+ PZ Cell" = "#D39200",
  "Proliferating PZ Cell" = "#00C19F",
  "Early Chondrocyte" = "#DB72FB",
  "Chondrocyte" = "#F8766D"
)

legend_labels <- c(
  "Mesenchymal Cells",
  "Proliferating SZ/Disc Progenitor",
  expression(italic("Lgr5")^"+"~"SZ/Disc Progenitor"),
  expression(italic("Pthlh")^"+"~"SZ/Disc Progenitor"),
  expression(italic("Gli1")^"+"~"PZ Cell"),
  expression(italic("Acta2")^"+"~"PZ Cell"),
  "Proliferating PZ Cell",
  "Early Chondrocyte",
  "Chondrocyte"
)

# Fig 4A: UMAP
message("Generating Fig 4A (UMAP)...")
p_umap <- DimPlot(seurat_obj, reduction = "umap", label = FALSE, cols = custom_colors, pt.size = 1.2) +
  NoAxes() +
  ggtitle(NULL) +
  scale_color_manual(values = custom_colors, labels = legend_labels) +
  theme(
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(1.5, "lines")
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))
save_plot(p_umap, file.path(OUT_DIR, "Fig4/Fig4A_UMAP.png"), width = 8, height = 6)

# Fig 4B-H: FeaturePlots
features <- c("Prrx1", "Sox9", "Runx2", "Lgr5", "Pthlh", "Acan", "Col2a1", 
              "Dlx5", "Ptch1", "Gli1", "Col10a1", "Ihh", "Fgfr3", "Mki67", 
              "Acta2", "Spp1", "Pth1r")

message("Generating Fig 4B-H (FeaturePlots)...")
for (feat in features) {
  if (feat %in% rownames(seurat_obj)) {
    p_feat <- FeaturePlot(seurat_obj, features = feat, reduction = "umap", order = TRUE) +
      scale_color_manuscript() +
      theme_publication()
    save_plot(p_feat, file.path(OUT_DIR, paste0("Fig4/Fig4_FeaturePlot_", feat, ".png")), width = 6, height = 5)
  } else {
    warning("Feature not found: ", feat)
  }
}

message("Figure 4 Scripts Complete.")
