# 05_fig5_scvelo.py
# Reproduces Figure 5C (scVelo Velocity Stream Plot)
# Uses the pre-computed dynamical model from the original analysis h5ad.
# Injects biological cluster names from Seurat metadata.
#
# Run with: /home/owner02/.local/share/r-miniconda/envs/scRNAseq/bin/python 05_fig5_scvelo.py

import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

print("Processing Fig 5C (scVelo Velocity Stream)...")

BASE_DIR = os.getcwd()
H5AD_FILE = os.path.join(BASE_DIR, "data/raw/TMJ_E16_Chondro_scVelo/scVelo.h5ad")
META_FILE = os.path.join(BASE_DIR, "output/Fig5/E16_Chondro_meta.csv")
OUT_DIR = os.path.join(BASE_DIR, "output/Fig5")
os.makedirs(OUT_DIR, exist_ok=True)

# ── 1. Load pre-computed scVelo h5ad (dynamical model already fitted) ──
print("  Loading pre-computed scVelo h5ad...")
adata = scv.read(H5AD_FILE)
print("  Shape:", adata.shape)
print("  Layers:", list(adata.layers.keys()))

# ── 2. Inject biological names from Seurat metadata ──
print("  Injecting biological cluster names from Seurat...")
meta = pd.read_csv(META_FILE)
meta["bc_stripped"] = meta["barcode"].str.replace(r"-1$", "", regex=True)

# Build barcode -> BiologicalName mapping
bc_to_name = dict(zip(meta["bc_stripped"], meta["BiologicalName"]))

# Map onto adata (the h5ad barcodes lack the -1 suffix)
adata.obs["BiologicalName"] = adata.obs.index.map(bc_to_name)
# Ensure categorical
adata.obs["BiologicalName"] = pd.Categorical(adata.obs["BiologicalName"])
print("  Cluster distribution:")
print(adata.obs["BiologicalName"].value_counts().to_string())

# ── 3. Define color palette matching manuscript ──
# Order: PZ Cell-1(0), Early Chondrocyte(1), PZ Cell-3(2), Proliferating PZ Cell-2(3),
#         Chondrocyte(4), PZ Cell-2(5), Proliferating PZ Cell-1(6), Proliferating PZ Cell-3(7)
cluster_colors = {
    "PZ Cell-1": "#1f77b4",               # blue
    "Early Chondrocyte": "#ff7f0e",        # orange
    "PZ Cell-3": "#2ca02c",                # green
    "Proliferating PZ Cell-2": "#d62728",  # red
    "Chondrocyte": "#9467bd",              # purple
    "PZ Cell-2": "#8c564b",               # brown
    "Proliferating PZ Cell-1": "#e377c2",  # pink
    "Proliferating PZ Cell-3": "#7f7f7f",  # grey
}

# Set palette in adata.uns so scvelo uses it
categories = adata.obs["BiologicalName"].cat.categories.tolist()
adata.uns["BiologicalName_colors"] = [cluster_colors.get(c, "#333333") for c in categories]

# ── 4. Recompute velocity_umap for the stream plot ──
# The h5ad has velocity_umap pre-computed, but let's ensure it's properly projected
# for the current UMAP embedding
print("  Projecting velocity onto UMAP...")
scv.tl.velocity_graph(adata)

# ── 5. Generate velocity stream plot ──
print("  Generating Fig 5C (Velocity Stream)...")
fig, ax = plt.subplots(figsize=(10, 8))
scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="BiologicalName",
    legend_loc="right margin",
    density=3.0,
    arrow_size=2.0,
    linewidth=1.0,
    save=False,
    ax=ax,
    show=False,
    title=""
)

plt.tight_layout()
out_path = os.path.join(OUT_DIR, "Fig5C_scVelo_Stream.png")
plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print("  Saved:", out_path)
print("Figure 5C Script Complete.")
