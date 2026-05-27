# Created: 2026-04-20 14:45
# =============================================================================
# Phase 2 Part 5.2 — scanpy leiden on ArchR Harmony embedding (CM subset)
# -----------------------------------------------------------------------------
# Prereq: 03_CM_subset_harmony.R 완료 → ArchRSubset_CM/harmony_CM.tsv 생성됨
# Output: ArchRSubset_CM/leiden_CM.tsv  (cellName, leiden_res02..res10)
# Next:   05_inject_leiden.R 로 ArchRProject에 주입
# -----------------------------------------------------------------------------
# 실행 환경: scanpy + leidenalg 설치된 conda env (컨테이너 밖, Wynton 호스트)
#   module load CBI miniconda3   # 또는 기본 설정된 경로
#   conda activate <scanpy_env>
#   python 04_leiden_CM.py
# =============================================================================

import os
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

WORK_DIR = "/gladstone/theodoris/lab/bkim/multi_multi/archr_dar/ArchRSubset_CM"
IN_TSV   = os.path.join(WORK_DIR, "harmony_CM.tsv")
OUT_TSV  = os.path.join(WORK_DIR, "leiden_CM.tsv")

print(f"Loading Harmony embedding: {IN_TSV}")
X = pd.read_csv(IN_TSV, sep="\t", index_col=0)
print(f"  shape: {X.shape}")

adata = ad.AnnData(X=X.values)
adata.obs_names = X.index
adata.obsm["X_harmony"] = X.values

# kNN graph on Harmony dims
sc.pp.neighbors(adata, use_rep="X_harmony", n_neighbors=30, random_state=0)

# 논문 (Kathiriya-Rao 2025, ATAC Clusters_res0.4)에 맞춰 res=0.4를 주 클러스터로.
# 이웃 resolution도 같이 계산해 sanity check / 대안 선택 가능하게.
resolutions = [0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
for r in resolutions:
    key = f"leiden_res{int(r*10):02d}"
    sc.tl.leiden(adata, resolution=r, random_state=0, key_added=key)
    n = adata.obs[key].nunique()
    print(f"  {key}: {n} clusters")

# 주 클러스터를 'leiden'으로 복사 (ArchR 주입 편의)
adata.obs["leiden"] = adata.obs["leiden_res04"]

out = adata.obs.copy()
out.index.name = "cellName"
out.to_csv(OUT_TSV, sep="\t")
print(f"\nSaved: {OUT_TSV}")
print("\nleiden_res04 cluster size:")
print(adata.obs["leiden_res04"].value_counts().sort_index().to_string())
