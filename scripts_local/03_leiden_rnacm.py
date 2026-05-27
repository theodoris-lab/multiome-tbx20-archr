# Created: 2026-04-23 13:25
# Updated: 2026-04-23 13:25
"""
Phase 4 of DAR_analysis_plan_rnaCM_0423.

scanpy leiden clustering on ArchR Harmony embedding (rnaCM subset).
Mirrors scripts_local/04_leiden_CM_local.py for the original run.

Input : outputs/dar_rnaCM_inputs_0423/harmony_rnaCM.tsv   (scp'd from Wynton)
Output: outputs/dar_rnaCM_inputs_0423/leiden_rnaCM.tsv    (cellName, leiden_res04, UMAP1/2)
        outputs/dar_rnaCM_inputs_0423/leiden_sweep_rnaCM.tsv
        outputs/dar_rnaCM_inputs_0423/leiden_rnaCM_umap.png
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3
sc.settings.seed = 0

PROJECT = Path("/Users/bkim/vscode/V2_multiome_2026-04-14")
INPUT = PROJECT / "outputs/dar_rnaCM_inputs_0423/harmony_rnaCM.tsv"
OUTDIR = PROJECT / "outputs/dar_rnaCM_inputs_0423"
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
    sweep[key] = adata.obs[key].astype(str)
    n = adata.obs[key].nunique()
    print(f"  res={r}: {n} clusters")

adata.obs["UMAP1"] = adata.obsm["X_umap"][:, 0]
adata.obs["UMAP2"] = adata.obsm["X_umap"][:, 1]

primary = adata.obs[["leiden_res04", "UMAP1", "UMAP2"]].copy()
primary.index.name = "cellName"
primary.to_csv(OUTDIR / "leiden_rnaCM.tsv", sep="\t")

sweep_df = pd.DataFrame(sweep, index=adata.obs.index)
sweep_df.index.name = "cellName"
sweep_df.to_csv(OUTDIR / "leiden_sweep_rnaCM.tsv", sep="\t")

fig, ax = plt.subplots(figsize=(6, 5))
clusters = sorted(adata.obs["leiden_res04"].unique(), key=lambda x: int(x))
for c in clusters:
    m = adata.obs["leiden_res04"] == c
    ax.scatter(adata.obs.loc[m, "UMAP1"], adata.obs.loc[m, "UMAP2"],
               s=2, alpha=0.6, label=f"C{int(c)+1}", rasterized=True)
ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
ax.set_title(f"rnaCM Leiden res=0.4 ({len(clusters)} clusters, {adata.n_obs:,} cells)")
ax.legend(markerscale=3, fontsize=7, loc="best")
plt.tight_layout()
plt.savefig(OUTDIR / "leiden_rnaCM_umap.png", dpi=150)
plt.close()

print(f"\n[wrote]")
print(f"  {OUTDIR / 'leiden_rnaCM.tsv'}")
print(f"  {OUTDIR / 'leiden_sweep_rnaCM.tsv'}")
print(f"  {OUTDIR / 'leiden_rnaCM_umap.png'}")
