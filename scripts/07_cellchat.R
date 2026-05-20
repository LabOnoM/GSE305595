# 07_cellchat.R
# Reproduces CellChat panels for Figures 2, 3, and 4

source("scripts/00_config.R")
library(CellChat)

options(stringsAsFactors = FALSE)
# Use mouse database
CellChatDB <- CellChatDB.mouse

run_cellchat_pipeline <- function(rds_path, output_name, title, pathways.show) {
  message("Processing CellChat for: ", output_name)
  
  seurat_obj <- readRDS(rds_path)
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Ensure Idents are correctly set
  if (is.null(Idents(seurat_obj))) {
    Idents(seurat_obj) <- seurat_obj$NewName
  }
  
  message("  Creating CellChat object...")
  cellchat <- createCellChat(object = seurat_obj, group.by = "ident")
  
  # Set the database
  cellchat@DB <- CellChatDB
  
  # Subset data
  message("  Subsetting signaling data...")
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # Compute communication probability
  message("  Computing communication probability...")
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
  # Filter out communications with few cells
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # Infer cellular communication network
  message("  Inferring communication network...")
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # Visualize aggregate network (Circle plot)
  message("  Generating aggregate network plot...")
  out_file_circle <- file.path(OUT_DIR, "CellChat", paste0(output_name, "_Aggregate_Circle.png"))
  png(out_file_circle, width = 800, height = 800, res = 150)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                   label.edge= F, title.name = paste(title, "- Number of interactions"))
  dev.off()
  
  # Generate pathway specific plots
  message("  Generating pathway specific plots...")
  available_pathways <- cellchat@netP$pathways
  pathways.show <- intersect(pathways.show, available_pathways)
  
  for (pathway in pathways.show) {
    out_file_pathway <- file.path(OUT_DIR, "CellChat", paste0(output_name, "_Pathway_", pathway, ".png"))
    png(out_file_pathway, width = 800, height = 800, res = 150)
    netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
    dev.off()
    
    # Dot plot or contribution
    out_file_contrib <- file.path(OUT_DIR, "CellChat", paste0(output_name, "_Contrib_", pathway, ".png"))
    p_contrib <- netAnalysis_contribution(cellchat, signaling = pathway) + theme_publication()
    save_plot(p_contrib, out_file_contrib, width = 6, height = 4)
  }
  
  message("  CellChat complete for ", output_name)
  return(cellchat)
}

# Define pathways to show based on manuscript
pathways_13w <- c("WNT", "BMP", "TGFb", "FGF", "IGF")
pathways_P0 <- c("PTH", "IHH", "WNT", "BMP", "TGFb", "FGF", "IGF")
pathways_E16 <- c("PTH", "IHH", "WNT", "BMP", "TGFb", "FGF", "IGF", "NOTCH", "VEGF")

# Ensure CellChat uses groupSize from data
groupSize <- NULL # will be calculated internally by CellChat functions when needed, but netVisual_circle uses vertex.weight directly. 
# We'll let CellChat handle the defaults if vertex.weight = as.numeric(table(cellchat@idents))
# To fix the scope issue for vertex.weight in the function:
# We should redefine run_cellchat_pipeline to calculate groupSize
run_cellchat_pipeline <- function(rds_path, output_name, title, pathways.show, annot_map) {
  message("Processing CellChat for: ", output_name)
  seurat_obj <- readRDS(rds_path)
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Apply Biological Annotation Map
  seurat_obj <- apply_annotation_map(seurat_obj, annot_map)
  seurat_obj$ident <- Idents(seurat_obj)
  
  cellchat <- createCellChat(object = seurat_obj, group.by = "ident")
  cellchat@DB <- CellChatDB.mouse
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  groupSize <- as.numeric(table(cellchat@idents))
  
  out_file_circle <- file.path(OUT_DIR, "CellChat", paste0(output_name, "_Aggregate_Circle.png"))
  png(out_file_circle, width = 800, height = 800, res = 150)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                   label.edge= F, title.name = paste(title, "- Interactions"))
  dev.off()
  
  available_pathways <- cellchat@netP$pathways
  pathways.show <- intersect(pathways.show, available_pathways)
  
  for (pathway in pathways.show) {
    out_file_pathway <- file.path(OUT_DIR, "CellChat", paste0(output_name, "_Pathway_", pathway, ".png"))
    png(out_file_pathway, width = 800, height = 800, res = 150)
    netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
    dev.off()
    
    out_file_contrib <- file.path(OUT_DIR, "CellChat", paste0(output_name, "_Contrib_", pathway, ".png"))
    p_contrib <- netAnalysis_contribution(cellchat, signaling = pathway) + theme_publication()
    ggsave(out_file_contrib, plot = p_contrib, width = 6, height = 4, dpi = 150, bg = "white")
  }
  return(cellchat)
}

run_cellchat_pipeline(DATA_13W_MES, "Fig2G_13w", "13w Mesenchymal", pathways_13w, annotation_map_13w_biological)
run_cellchat_pipeline(DATA_P0_MES, "Fig3G_P0", "P0 Mesenchymal", pathways_P0, annotation_map_P0)
run_cellchat_pipeline(DATA_E16_MES, "Fig4I_E16.5", "E16.5 Mesenchymal", pathways_E16, annotation_map_E16)

message("CellChat Scripts Complete.")
