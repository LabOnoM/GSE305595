# 05_fig5_stream.py
# Reproduces Figure 5G-H (STREAM Trajectory)

import stream as st
import matplotlib.pyplot as plt
import os

print("Processing Fig 5G-H (STREAM)...")
BASE_DIR = os.getcwd()
DATA_FILE = os.path.join(BASE_DIR, "data/raw/TMJ_E16_Chondro_STREAM/STREAM_out.pkl")
OUT_DIR = os.path.join(BASE_DIR, "output/Fig5")

os.makedirs(OUT_DIR, exist_ok=True)

# 1. Load STREAM data
print("  Loading STREAM dataset...")
adata = st.read(DATA_FILE)

import pandas as pd
label_map = {
    '0': 'PZ Cell-1',
    '1': 'Early Chondrocyte',
    '2': 'PZ Cell-3',
    '3': 'Proliferating PZ Cell-2',
    '4': 'Chondrocyte',
    '5': 'PZ Cell-2',
    '6': 'Proliferating PZ Cell-1',
    '7': 'Proliferating PZ Cell-3'
}
label_map.update({int(k): v for k, v in label_map.items()})
adata.obs['label_str'] = adata.obs['label'].map(label_map).astype(str)

# Map numeric labels to the strings in adata.obs['label']
# STREAM requires categorical variables to be string arrays for plot_flat_tree
adata.obs['label'] = adata.obs['label_str']

custom_colors = {
    "PZ Cell-1": "#377EB8",
    "PZ Cell-2": "#A65628",
    "PZ Cell-3": "#4DAF4A",
    "Proliferating PZ Cell-1": "#F781BF",
    "Proliferating PZ Cell-2": "#E41A1C",
    "Proliferating PZ Cell-3": "#999999",
    "Early Chondrocyte": "#FF7F00",
    "Chondrocyte": "#984EA3"
}

# STREAM expects a dictionary for colors, NOT a list (unlike scanpy)
adata.uns['label_color'] = custom_colors
# Let's also enforce legend ordering by hacking STREAM's st.plot_flat_tree kwarg or we just accept seaborn's default order, which is fine since colors are correctly mapped.

# 2. Fig 5G: Flat tree plot
print("  Generating Fig 5G (plot_flat_tree)...")
st.plot_flat_tree(adata, color=['label', 'S0_pseudotime'], save_fig=True, fig_path=OUT_DIR, fig_name='Fig5G_STREAM_Branches.png')

import shutil
# flat tree saves exactly what we name it if fig_name is passed, but it might append suffix or similar depending on STREAM version.
# Let's check how it saves it:
fig5g_path_label = os.path.join(OUT_DIR, "flat_tree_label.png")
fig5g_path_pseudo = os.path.join(OUT_DIR, "flat_tree_S0_pseudotime.png")
if os.path.exists(fig5g_path_label):
    shutil.move(fig5g_path_label, os.path.join(OUT_DIR, "Fig5G_STREAM_Branches.png"))
if os.path.exists(fig5g_path_pseudo):
    shutil.move(fig5g_path_pseudo, os.path.join(OUT_DIR, "Fig5G_STREAM_Pseudotime.png"))


# 3. Fig 5H: Marker gene expression and Biological Labels on Subway Map
print("  Generating Fig 5H (plot_stream_sc / Subway Map)...")
# First plot the biological labels
st.plot_stream_sc(adata, root='S0', color=['label'], save_fig=True, fig_path=OUT_DIR, fig_format='png')

fig5h_src = os.path.join(OUT_DIR, "S0", "stream_sc_label.png")
if os.path.exists(fig5h_src):
    shutil.move(fig5h_src, os.path.join(OUT_DIR, "Fig5H_STREAM_SubwayMap.png"))

# Next plot the key markers
markers = ['Prrx1', 'Sox9', 'Crabp1', 'Runx2', 'Col2a1']
available_markers = [m for m in markers if m in adata.var_names]
if available_markers:
    st.plot_stream_sc(adata, root='S0', color=available_markers, save_fig=True, 
                      fig_path=OUT_DIR, fig_format='png')
    
    fig5h_dir = os.path.join(OUT_DIR, "S0")
    for m in available_markers:
        src = os.path.join(fig5h_dir, f"stream_sc_{m.lower()}.png")
        if os.path.exists(src):
            shutil.copy(src, os.path.join(OUT_DIR, f"Fig5H_STREAM_{m}.png"))

# 4. Export obs for Seurat pseudotime UMAP mapping
print("  Exporting Pseudotime for Fig 5D...")
adata.obs.to_csv(os.path.join(OUT_DIR, "stream_obs_metadata.csv"))

print("Figure 5G-H Script Complete.")
