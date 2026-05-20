# 01_fig1_FEA.R
# Reproduces Figure 1D (E16.5), 1E (P0), 1F (13w) GO Enrichment DotPlots

source("scripts/00_config.R")
library(clusterProfiler)
library(org.Mm.eg.db)

message("Loading integrated Atlas object for cross-stage comparison...")
seurat_obj <- readRDS(DATA_INT_ALL)
DefaultAssay(seurat_obj) <- "RNA"

# The manuscript specifies: 
# "To evaluate developmental stage-specific transcriptional programs... we performed DEG and GO enrichment analyses on Sox9+ cartilage progenitor and chondrocyte populations from each stage."
# Progenitors (Prrx1+Sox9+): Clusters 2, 3, 7, 8, 9
# Chondrocytes (Prrx1-Sox9+): Clusters 0, 1, 4, 10, 11, 12, 17
target_clusters <- c(0, 1, 2, 3, 4, 7, 8, 9, 10, 11, 12, 17)
seurat_sub <- subset(seurat_obj, idents = target_clusters)

# Set identities to the stage to compare E16.5 vs P0 vs 13w
Idents(seurat_sub) <- seurat_sub$orig.ident

message("Finding markers between ages within the Cartilage/Chondrocyte lineage...")
markers <- FindAllMarkers(seurat_sub, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers_sig <- markers %>% filter(p_val_adj < 0.05)

# Map to Entrez
gene_symbols <- unique(markers_sig$gene)
entrez_map <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
markers_sig <- markers_sig %>% left_join(entrez_map, by = c("gene" = "SYMBOL")) %>% filter(!is.na(ENTREZID))

# Define titles and output names
meta_info <- list(
  "E16-5WT" = list(title = "GO BP Enrichment - E16.5", out = "Fig1D_E16.5"),
  "0dWT"    = list(title = "GO BP Enrichment - P0",    out = "Fig1E_P0"),
  "13wWT"   = list(title = "GO BP Enrichment - 13w",   out = "Fig1F_13w")
)

# Define Exact Terms based on Manuscript
terms_e16 <- c(
  "regulation of epithelial cell proliferation",
  "chromatin remodeling",
  "RNA splicing, via transesterification reactions",
  "positive regulation of cell cycle phase transition",
  "nucleosome assembly"
)

terms_p0 <- c(
  "skeletal system morphogenesis",
  "extracellular matrix organization",
  "cartilage development",
  "endochondral bone morphogenesis",
  "genomic imprinting"
)

terms_13w <- c(
  "extracellular matrix organization",
  "wound healing",
  "ossification",
  "collagen metabolic process",
  "collagen fibril organization"
)

create_dotplot <- function(fea_data, target_terms, title) {
  # Exact string matching
  plot_data <- fea_data[fea_data$Description %in% target_terms, ]
  
  # Set order to match the target_terms precisely
  plot_data$Description <- factor(plot_data$Description, levels = rev(target_terms))
  
  # Convert GeneRatio to numeric
  plot_data$GeneRatio_Num <- sapply(plot_data$GeneRatio, function(x) {
    parts <- as.numeric(strsplit(x, "/")[[1]])
    parts[1] / parts[2]
  })
  
  # Calculate count
  plot_data$Count <- sapply(plot_data$BgRatio, function(x) {
    as.numeric(strsplit(x, "/")[[1]][1])
  })
  
  p <- ggplot(plot_data, aes(x = GeneRatio_Num, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue") +
    theme_publication() +
    labs(title = title, x = "GeneRatio", y = "")
  
  return(p)
}

# Define stage to terms mapping
terms_map <- list(
  "E16-5WT" = terms_e16,
  "0dWT"    = terms_p0,
  "13wWT"   = terms_13w
)

for (stage in names(meta_info)) {
  message("Processing: ", stage)
  stage_genes <- markers_sig %>% filter(cluster == stage) %>% pull(ENTREZID)
  
  ego <- enrichGO(gene          = stage_genes,
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)
  
  if (!is.null(ego)) {
    target_terms <- terms_map[[stage]]
    
    # Check what is found
    found_terms <- ego@result$Description[ego@result$Description %in% target_terms]
    if (length(found_terms) > 0) {
      p <- create_dotplot(ego@result, target_terms, meta_info[[stage]]$title)
      
      out_file <- file.path(OUT_DIR, "Fig1", paste0(meta_info[[stage]]$out, "_GO_DotPlot.png"))
      save_plot(p, out_file, width = 10, height = 6)
    } else {
      message("Warning: None of the target GO terms were found for ", stage)
    }
  }
}

message("Figure 1 FEA Complete.")
