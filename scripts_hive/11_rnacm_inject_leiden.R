# Created: 2026-04-23 13:25
# Updated: 2026-04-23 16:55
# ==============================================================================
# Phase 5 — Inject scanpy leiden labels into ArchRSubset_rnaCM, set up
#           Clusters + ClusterByGenotype. Mirrors 05_inject_leiden.R.
# Prereq: local scanpy leiden 결과 (leiden_rnaCM.tsv) Wynton 업로드 완료.
# Input:  ArchRSubset_rnaCM/leiden_rnaCM.tsv  (row index = cellName, cols = leiden_res04, ...)
# Output: ArchRSubset_rnaCM/ 에 Clusters, ClusterByGenotype 주입, UMAP plots
#         outputs/tables/Leiden_vs_Genotype_rnaCM.tsv
# Next:   04_rnaCM_dar_peak_motif.R
# ==============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(data.table)
})

.libPaths(.libPaths()[!grepl("/wynton/home/.*/R/x86_64", .libPaths())])

WORK_DIR <- "/gladstone/theodoris/lab/bkim/multi_multi/archr_dar_rnaCM"
setwd(WORK_DIR)
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)

addArchRThreads(threads = 8)
addArchRGenome("hg38")

projSub <- loadArchRProject(file.path(WORK_DIR, "ArchRSubset_rnaCM"), showLogo = FALSE)

lc_path <- file.path(WORK_DIR, "ArchRSubset_rnaCM", "leiden_rnaCM.tsv")
stopifnot(file.exists(lc_path))
lc <- fread(lc_path, data.table = FALSE)
rownames(lc) <- lc$cellName

ord <- match(projSub$cellNames, rownames(lc))
stopifnot(!any(is.na(ord)))

leiden_int <- as.integer(as.character(lc[ord, "leiden_res04"]))
projSub$Clusters          <- paste0("C", leiden_int + 1L)
projSub$ClusterByGenotype <- paste0(projSub$Clusters, "_x_", projSub$Genotype)

# Note: 보조 leiden_res 컬럼은 추가하지 않음 — ArchRProject 가 [[<- 미지원.
# 필요 시 addCellColData() 로 주입 가능.

cat("===== Clusters × Genotype (leiden res=0.4, rnaCM) =====\n")
tab <- table(projSub$Clusters, projSub$Genotype)
print(tab)
fwrite(as.data.frame.matrix(tab),
       file.path(WORK_DIR, "outputs/tables/Leiden_vs_Genotype_rnaCM.tsv"),
       sep = "\t", row.names = TRUE)

cat("\n===== ClusterByGenotype n cells =====\n")
cbg_tab <- sort(table(projSub$ClusterByGenotype), decreasing = TRUE)
print(cbg_tab)

pl <- plotEmbedding(projSub, colorBy = "cellColData", name = "Clusters",
                    embedding = "UMAP", plotAs = "points", rastr = FALSE)
pg <- plotEmbedding(projSub, colorBy = "cellColData", name = "Genotype",
                    embedding = "UMAP", plotAs = "points", rastr = FALSE,
                    pal = c(Wild = "blue", Het = "orange", Hom = "red"))
pcbg <- plotEmbedding(projSub, colorBy = "cellColData", name = "ClusterByGenotype",
                      embedding = "UMAP", plotAs = "points", rastr = FALSE)
plotPDF(pl, pg, pcbg, name = "rnaCM_UMAP_leiden.pdf",
        ArchRProj = projSub, addDOC = FALSE, width = 6, height = 6)

saveArchRProject(projSub)

cat("\n완료. Leiden_vs_Genotype_rnaCM.tsv 검토 후 04_rnaCM_dar_peak_motif.R 내 DAR 그룹 수정.\n")
cat("  - WT 많은 cluster(들) → bgdGroups\n")
cat("  - Het 쏠린 cluster   → HET_FG\n")
cat("  - Hom 쏠린 cluster   → HOM_FG\n")
cat("기대: 원본 run과 유사한 분포 (C1/C2 WT+Het, C3 Hom). cluster IDs는 다를 수 있음.\n")
