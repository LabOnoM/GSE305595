# 06_fig6e_crabp1.R
# Reproduces Figure 6E DotPlot (Crabp1+ vs Crabp1- comparison)

source("scripts/00_config.R")

message("Processing Fig 6E (Crabp1 Validation)...")
seurat_obj <- readRDS(DATA_E16_CHONDRO)
DefaultAssay(seurat_obj) <- "SCT"

# 1. Define Crabp1 expression groups
message("  Defining Crabp1-positive vs Crabp1-negative groups...")
crabp1_exp <- GetAssayData(seurat_obj, slot = "data")["Crabp1", ]
seurat_obj$Crabp1_status <- ifelse(crabp1_exp > 0, "Crabp1-positive", "Crabp1-negative")

# Enforce factor ordering so Crabp1-positive is on top
seurat_obj$Crabp1_status <- factor(seurat_obj$Crabp1_status, levels = c("Crabp1-negative", "Crabp1-positive"))

# 2. Generate DotPlot
message("  Generating Fig 6E (DotPlot)...")
features <- c("Prrx1", "Runx2", "Acta2", "Gli1", "Ptch1", "Sox9", "Acan", "Col2a1", "Col10a1")

# Ensure features are present
features_present <- features[features %in% rownames(seurat_obj)]

p_dotplot <- DotPlot(seurat_obj, features = features_present, group.by = "Crabp1_status") +
  theme_publication() +
  scale_color_gradient(low = "black", high = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Figure 6E: Crabp1 Validation")

save_plot(p_dotplot, file.path(OUT_DIR, "Fig6/Fig6E_Crabp1_DotPlot.png"), width = 7, height = 5)

message("Figure 6E Script Complete.")
