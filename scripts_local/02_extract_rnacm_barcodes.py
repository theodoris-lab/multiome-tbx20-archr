# Created: 2026-04-23 13:20
# Updated: 2026-04-23 13:20
"""
Compute scanpy leiden_res0.1 on rna.h5ad, select CM clusters, and export
barcodes + metadata for ArchR subsetting (scripts_hive/09_rnacm_subset.R).

CM clusters: leiden_res0.1 clusters {"0", "7"} — update if re-clustering changes IDs.

Outputs
-------
outputs/dar_rnaCM_inputs_0423/rnaCM_barcodes.tsv        (cellName)
outputs/dar_rnaCM_inputs_0423/rnaCM_cell_metadata.tsv   (rich metadata)
outputs/dar_rnaCM_inputs_0423/rnaCM_leiden_res01.tsv    (all cells, for audit)
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3
sc.settings.seed = 0

PROJECT = Path(__file__).resolve().parents[1]
H5AD = PROJECT / "data/anndata/rna.h5ad"
OUTDIR = PROJECT / "outputs/dar_rnaCM_inputs_0423"
OUTDIR.mkdir(parents=True, exist_ok=True)

CM_CLUSTERS = {"0", "7"}
EXPECTED_CM = 17408

print(f"[read] {H5AD}")
adata = ad.read_h5ad(H5AD)
print(f"  n_obs={adata.n_obs:,}  obsm={list(adata.obsm.keys())}")

assert "X_pca" in adata.obsm, "X_pca missing from rna.h5ad"

print("[neighbors] X_pca, n_neighbors=15, n_pcs=30")
sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15, n_pcs=30, random_state=0)

print("[leiden] resolution=0.1  key=leiden_res0.1")
sc.tl.leiden(adata, resolution=0.1, key_added="leiden_res0.1", random_state=0)

counts = adata.obs["leiden_res0.1"].value_counts().sort_index()
print("  cluster sizes:")
for k, v in counts.items():
    flag = "  <- CM" if str(k) in CM_CLUSTERS else ""
    print(f"    {k}: {v:,}{flag}")

cm_mask = adata.obs["leiden_res0.1"].astype(str).isin(CM_CLUSTERS)
n_cm = int(cm_mask.sum())
print(f"\n[CM] {n_cm:,} cells  (expected ~{EXPECTED_CM:,})")
if abs(n_cm - EXPECTED_CM) > 500:
    print(f"  WARNING: CM count drift > 500 from expectation ({EXPECTED_CM:,}).")

cm = adata.obs.loc[cm_mask, [
    "arc_library", "gex_barcode", "multiseq_group", "multiseq_sample",
    "seurat_clusters", "leiden_res0.1",
]].copy()
cm.index.name = "cellName"
cm.reset_index(inplace=True)

sample_prefix = cm["cellName"].str.split("#").str[0]
assert sample_prefix.isin({"G1", "G2", "G3"}).all(), \
    f"cellName prefix not in {{G1,G2,G3}}: unique={sample_prefix.unique()}"
print(f"[format ok] cellName prefix ∈ {{G1,G2,G3}}, n_samples={sample_prefix.value_counts().to_dict()}")

bc_path = OUTDIR / "rnaCM_barcodes.tsv"
meta_path = OUTDIR / "rnaCM_cell_metadata.tsv"
all_path = OUTDIR / "rnaCM_leiden_res01.tsv"

cm[["cellName"]].to_csv(bc_path, sep="\t", index=False)
cm.to_csv(meta_path, sep="\t", index=False)

full = adata.obs[["leiden_res0.1"]].copy()
full.index.name = "cellName"
full.to_csv(all_path, sep="\t")

print("\n[wrote]")
print(f"  {bc_path}  ({len(cm):,} lines)")
print(f"  {meta_path}")
print(f"  {all_path}  (all {adata.n_obs:,} cells)")

print("\n[genotype breakdown of CM]")
print(cm["multiseq_group"].value_counts().to_string())
print("\n[timepoint breakdown of CM]")
print(cm["multiseq_sample"].value_counts().sort_index().to_string())
