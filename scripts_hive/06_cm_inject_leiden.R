# Created: 2026-04-20 14:45
# ==============================================================================
# Phase 2 Part 5.3 — Inject scanpy leiden labels into ArchRSubset_CM
# ------------------------------------------------------------------------------
# Prereq: 04_leiden_CM.py 완료 → ArchRSubset_CM/leiden_CM.tsv
# Output: ArchRSubset_CM/ 에 Clusters + ClusterByGenotype 저장
#         outputs/plots/CM_subset_UMAP_leiden.pdf
#         outputs/tables/Leiden_vs_Genotype.tsv
# Next:   06_dar_peak_motif.R (qsub 배치)
# ==============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(data.table)
})

.libPaths(.libPaths()[!grepl("/wynton/home/.*/R/x86_64", .libPaths())])

WORK_DIR <- "/gladstone/theodoris/lab/bkim/multi_multi/archr_dar"
setwd(WORK_DIR)

addArchRThreads(threads = 8)
addArchRGenome("hg38")

projCM <- loadArchRProject(path = file.path(WORK_DIR, "ArchRSubset_CM"))

lc <- fread(file.path(WORK_DIR, "ArchRSubset_CM", "leiden_CM.tsv"),
            data.table = FALSE)
rownames(lc) <- lc$cellName

# 순서 맞춤 + ArchR 관행 라벨 (C1, C2, …)
ord <- match(projCM$cellNames, rownames(lc))
stopifnot(!any(is.na(ord)))
leiden_int <- as.integer(as.character(lc[ord, "leiden_res04"]))
projCM$Clusters          <- paste0("C", leiden_int + 1L)   # 0-based → 1-based
projCM$ClusterByGenotype <- paste0(projCM$Clusters, "_x_", projCM$Genotype)

# 참고용 다른 resolution도 같이 주입
for (r_key in grep("^leiden_res", colnames(lc), value = TRUE)) {
  ival <- as.integer(as.character(lc[ord, r_key]))
  projCM[[paste0("Leiden_", r_key)]] <- paste0("C", ival + 1L)
}

cat("===== Clusters × Genotype (leiden res=0.4) =====\n")
tab <- table(projCM$Clusters, projCM$Genotype)
print(tab)
fwrite(as.data.frame.matrix(tab),
       file.path(WORK_DIR, "outputs/tables/Leiden_vs_Genotype.tsv"),
       sep = "\t", row.names = TRUE)

cat("\n===== ClusterByGenotype n cells =====\n")
cbg_tab <- sort(table(projCM$ClusterByGenotype), decreasing = TRUE)
print(cbg_tab)

# UMAP colored by leiden + genotype
pl <- plotEmbedding(projCM, colorBy = "cellColData", name = "Clusters",
                    embedding = "UMAP", plotAs = "points", rastr = FALSE)
pg <- plotEmbedding(projCM, colorBy = "cellColData", name = "Genotype",
                    embedding = "UMAP", plotAs = "points", rastr = FALSE,
                    pal = c(Wild = "blue", Het = "orange", Hom = "red"))
pcbg <- plotEmbedding(projCM, colorBy = "cellColData", name = "ClusterByGenotype",
                      embedding = "UMAP", plotAs = "points", rastr = FALSE)
plotPDF(pl, pg, pcbg, name = "CM_subset_UMAP_leiden.pdf",
        ArchRProj = projCM, addDOC = FALSE, width = 6, height = 6)

saveArchRProject(projCM)

cat("\n완료. Leiden_vs_Genotype.tsv 검토 후 DAR 그룹 선정:\n")
cat(" - WT 많은 cluster(들) → bgdGroups\n")
cat(" - Het 쏠린 cluster  → HET_FG (useGroups)\n")
cat(" - Hom 쏠린 cluster  → HOM_FG\n")
cat("최소 cluster당 500 cells 있어야 addGroupCoverages 안정적.\n")
