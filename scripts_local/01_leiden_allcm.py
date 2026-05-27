# Created: 2026-04-20 15:10
# Updated: 2026-04-23 10:25
"""
scanpy leiden clustering on ArchR Harmony embedding (CM subset).

Input : outputs/dar_archr_inputs_0420/harmony_CM.tsv  (cells x 30 Harmony dims)
Output: outputs/dar_archr_inputs_0420/leiden_CM.tsv            (cell, leiden_res04, UMAP1, UMAP2)
        outputs/dar_archr_inputs_0420/leiden_sweep_CM.tsv      (cell, sweep of resolutions)
        outputs/dar_archr_inputs_0420/leiden_CM_umap.png       (diagnostic UMAP)
"""
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
from pathlib import Path

sc.settings.verbosity = 3
sc.settings.seed = 0

PROJECT = Path(__file__).resolve().parents[1]
INPUT   = PROJECT / "outputs/dar_archr_inputs_0420/harmony_CM.tsv"
OUTDIR  = PROJECT / "outputs/dar_archr_inputs_0420"
OUTDIR.mkdir(parents=True, exist_ok=True)

print(f"[read] {INPUT}")
H = pd.read_csv(INPUT, sep="\t", index_col=0)
print(f"  shape: {H.shape}  (cells x dims)")

adata = ad.AnnData(
    X=np.zeros((H.shape[0], 1), dtype="float32"),
    obs=pd.DataFrame(index=H.index),
)
adata.obsm["X_harmony"] = H.values.astype("float32")

print("[neighbors] n_neighbors=30  on X_harmony")
sc.pp.neighbors(adata, n_neighbors=30, use_rep="X_harmony", random_state=0)

print("[umap]")
sc.tl.umap(adata, random_state=0)

print("[leiden] resolution=0.4  key=leiden_res04  (primary for DAR)")
sc.tl.leiden(adata, resolution=0.4, random_state=0, key_added="leiden_res04")

sweep = {}
for r in [0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]:
    key = f"leiden_r{str(r).replace('.', '')}"
    sc.tl.leiden(adata, resolution=r, random_state=0, key_added=key)
    sweep[key] = adata.obs[key].copy()
    print(f"  res={r}: {adata.obs[key].nunique()} clusters")

primary = pd.DataFrame({
    "cell": adata.obs.index,
    "leiden_res04": adata.obs["leiden_res04"].astype(int).values,
    "UMAP1": adata.obsm["X_umap"][:, 0],
    "UMAP2": adata.obsm["X_umap"][:, 1],
})
primary_path = OUTDIR / "leiden_CM.tsv"
primary.to_csv(primary_path, sep="\t", index=False)
print(f"[write] {primary_path}  ({len(primary)} cells)")

sweep_df = pd.DataFrame({k: v.astype(int).values for k, v in sweep.items()}, index=adata.obs.index)
sweep_df.index.name = "cell"
sweep_path = OUTDIR / "leiden_sweep_CM.tsv"
sweep_df.to_csv(sweep_path, sep="\t")
print(f"[write] {sweep_path}")

fig, ax = plt.subplots(figsize=(6, 5))
sc.pl.umap(adata, color="leiden_res04", legend_loc="on data", ax=ax, show=False, frameon=False, title="Leiden res=0.4")
fig.tight_layout()
png_path = OUTDIR / "leiden_CM_umap.png"
fig.savefig(png_path, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"[write] {png_path}")

counts = adata.obs["leiden_res04"].value_counts().sort_index()
print("\n[summary] leiden_res04 cluster sizes:")
print(counts.to_string())
print(f"\nTotal cells: {adata.n_obs}  |  primary clusters: {counts.shape[0]}")
print("done.")
