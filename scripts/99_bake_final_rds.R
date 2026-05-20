# 99_bake_final_rds.R
# Script to build a fully standalone, final annotated dataset for the reproducibility bundle.
# This prevents scripts from needing to re-apply annotations and factor levels at runtime.

source("scripts/00_config.R")
library(Seurat)

# Create output directory for baked RDS files
BAKED_DIR <- "/mnt/Public/00.Database/TosaCellReports/FigureReproduction/data/Final_Annotated_RDS"
if (!dir.exists(BAKED_DIR)) dir.create(BAKED_DIR, recursive = TRUE)

message("--- Baking P0 Dataset ---")
p0 <- readRDS(file.path(DATA_DIR, "TMJ_P0_Mes.rds"))
DefaultAssay(p0) <- "RNA"
p0 <- apply_annotation_map(p0, annotation_map_P0)
# Define exactly the desired order from Fig 3
desired_order_P0 <- c(
  "SZ Cell",
  "PZ Cell",
  "Acta2+ PZ Cell",
  "Proliferating PZ Cell",
  "Early Chondrocyte",
  "Chondrocyte",
  "Hypertrophic Chondrocyte",
  "Chondrocyte-associated Osteogenic cell"
)
Idents(p0) <- factor(Idents(p0), levels = desired_order_P0)
p0$BiologicalName <- factor(p0$BiologicalName, levels = desired_order_P0)
saveRDS(p0, file.path(BAKED_DIR, "P0_Mesenchymal_Annotated.rds"))
rm(p0); gc()

message("--- Baking E16.5 Mesenchymal Dataset ---")
e16 <- readRDS(file.path(DATA_DIR, "TMJ_E16_Mes.rds"))
DefaultAssay(e16) <- "RNA"
e16 <- apply_annotation_map(e16, annotation_map_E16)
desired_order_E16 <- c(
  "Mesenchymal Cells",
  "Proliferating SZ/Disc Progenitor",
  "Lgr5+ SZ/Disc Progenitor",
  "Pthlh+ SZ/Disc Progenitor",
  "Gli1+ PZ Cell",
  "Acta2+ PZ Cell",
  "Proliferating PZ Cell",
  "Early Chondrocyte",
  "Chondrocyte"
)
Idents(e16) <- factor(Idents(e16), levels = desired_order_E16)
e16$BiologicalName <- factor(e16$BiologicalName, levels = desired_order_E16)
saveRDS(e16, file.path(BAKED_DIR, "E16.5_Mesenchymal_Annotated.rds"))
rm(e16); gc()

message("--- Baking E16.5 Chondrogenic Dataset ---")
e16_chondro <- readRDS(file.path(DATA_DIR, "TMJ_E16_Chondro.rds"))
DefaultAssay(e16_chondro) <- "RNA"
e16_chondro <- apply_annotation_map(e16_chondro, annotation_map_E16_Chondro)
desired_order_chondro <- c(
  "PZ Cell-1",
  "PZ Cell-2",
  "PZ Cell-3",
  "Proliferating PZ Cell-1",
  "Proliferating PZ Cell-2",
  "Proliferating PZ Cell-3",
  "Early Chondrocyte",
  "Chondrocyte"
)
Idents(e16_chondro) <- factor(Idents(e16_chondro), levels = desired_order_chondro)
e16_chondro$BiologicalName <- factor(e16_chondro$BiologicalName, levels = desired_order_chondro)


saveRDS(e16_chondro, file.path(BAKED_DIR, "E16.5_Chondro_Annotated.rds"))
rm(e16_chondro); gc()

message("--- Baking 13w Dataset ---")
# Need raw data since baked data doesn't exist correctly yet
w13 <- readRDS(file.path(DATA_DIR, "TMJ_13w_Mes.rds"))
DefaultAssay(w13) <- "RNA"
# First get gene markers from CSV
w13 <- apply_csv_annotation(w13, file.path(DATA_DIR, "13w_NewName_mapping.csv"))
# apply_csv_annotation sets BiologicalName to the gene markers (NewName in CSV)
Idents(w13) <- w13$BiologicalName
# Then map gene markers to biological names
w13 <- RenameIdents(w13, annotation_map_13w_biological)
w13$BiologicalName <- Idents(w13)

desired_order_13w <- c(
  "SZ Cell",
  "SZ-associated Osteogenic cell",
  "Thy1+ PZ Cell",
  "Acta2+ PZ Cell",
  "Proliferating PZ cell",
  "Early Chondrocyte",
  "Chondrocyte",
  "Hypertrophic Chondrocyte",
  "Chondrocyte-associated Osteogenic cell"
)
Idents(w13) <- factor(Idents(w13), levels = desired_order_13w)
w13$BiologicalName <- factor(w13$BiologicalName, levels = desired_order_13w)

saveRDS(w13, file.path(BAKED_DATA_DIR, "13w_Mesenchymal_Annotated.rds"))
rm(w13); gc()

message("All datasets successfully baked and saved to: ", BAKED_DIR)
