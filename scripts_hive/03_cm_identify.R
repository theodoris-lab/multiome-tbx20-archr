# Created: 2026-04-20 14:45
# ==============================================================================
# Phase 2 Part 4 — CM cluster identification on all-cells ArchRProject
# ------------------------------------------------------------------------------
# Prereq: Phase 2.3 완료 (LSI + Harmony + Clusters(res=0.4) + UMAP 저장됨)
# Output: outputs/plots/CM_markers_UMAP.pdf
#         outputs/tables/Cluster_vs_SeuratClusters.tsv
#         outputs/tables/Cluster_vs_Genotype.tsv
#         outputs/tables/GeneScore_median_by_cluster.tsv
# Next:   plot 검토 후 CM_CLUSTERS 결정 → 03_CM_subset_harmony.R 에 기록
# ==============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(ggplot2)
  library(data.table)
})

# 호스트 R 라이브러리 섞임 방지
.libPaths(.libPaths()[!grepl("/wynton/home/.*/R/x86_64", .libPaths())])

WORK_DIR <- "/gladstone/theodoris/lab/bkim/multi_multi/archr_dar"
setwd(WORK_DIR)
dir.create("outputs/plots",  recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

addArchRThreads(threads = 8)
addArchRGenome("hg38")

proj <- loadArchRProject(path = file.path(WORK_DIR, "ArchRProject"))
cat("Cells:", length(proj$cellNames), "\n")
cat("Clusters:\n"); print(table(proj$Clusters))

# ---- 1) MAGIC impute weights --------------------------------------------------
if (is.null(getImputeWeights(proj))) {
  proj <- addImputeWeights(proj)
  saveArchRProject(proj)
}

# ---- 2) Gene score UMAP plots -------------------------------------------------
# CM 마커 + 비-CM (배제용) 마커
markerGenes <- c(
  "TNNT2","MYH6","MYH7","NPPA","PITX2","HAMP",       # CM
  "COL1A1","COL3A1","DCN","POSTN",                    # Fibroblast
  "EPCAM","CLDN4",                                    # Epithelial
  "TOP2A","CENPF",                                    # Dividing
  "NEFM","MSX1",                                      # Neural
  "AFP","SERPINA1","TTR"                              # Endoderm
)

p <- plotEmbedding(
  ArchRProj     = proj,
  colorBy       = "GeneScoreMatrix",
  name          = markerGenes,
  embedding     = "UMAP",
  imputeWeights = getImputeWeights(proj),
  plotAs        = "points",
  rastr         = FALSE,
  quantCut      = c(0.01, 0.95)
)
p2 <- lapply(p, function(x) {
  x + guides(color = "none", fill = "none") +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0,0,0,0), "cm"),
          axis.text = element_blank(), axis.ticks = element_blank())
})
plotPDF(plotList = p2, name = "CM_markers_UMAP.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# Clusters UMAP (라벨 확인용)
pc <- plotEmbedding(proj, colorBy = "cellColData", name = "Clusters",
                    embedding = "UMAP", plotAs = "points", rastr = FALSE)
pg <- plotEmbedding(proj, colorBy = "cellColData", name = "Genotype",
                    embedding = "UMAP", plotAs = "points", rastr = FALSE,
                    pal = c(Wild = "blue", Het = "orange", Hom = "red"))
ps <- plotEmbedding(proj, colorBy = "cellColData", name = "SeuratClusters",
                    embedding = "UMAP", plotAs = "points", rastr = FALSE)
plotPDF(pc, pg, ps, name = "UMAP_overview.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# ---- 3) Gene score 중앙값 by cluster ------------------------------------------
# GeneScoreMatrix → cluster별 median → TNNT2 등 상위 클러스터 식별
gsm <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
gsm_genes <- rowData(gsm)$name
keep <- gsm_genes %in% markerGenes
mat  <- assay(gsm)[keep, , drop = FALSE]
rownames(mat) <- gsm_genes[keep]

clust <- proj$Clusters
med_by_cluster <- sapply(sort(unique(clust)), function(cl) {
  apply(mat[, clust == cl, drop = FALSE], 1, median)
})
med_df <- as.data.frame(med_by_cluster)
med_df$gene <- rownames(med_df)
fwrite(med_df[, c("gene", sort(unique(clust)))],
       "outputs/tables/GeneScore_median_by_cluster.tsv", sep = "\t")

# ---- 4) Cluster × SeuratClusters, Cluster × Genotype 교차표 ------------------
tab_sc <- table(proj$Clusters, proj$SeuratClusters)
tab_gt <- table(proj$Clusters, proj$Genotype)
fwrite(as.data.frame.matrix(tab_sc),
       "outputs/tables/Cluster_vs_SeuratClusters.tsv",
       sep = "\t", row.names = TRUE)
fwrite(as.data.frame.matrix(tab_gt),
       "outputs/tables/Cluster_vs_Genotype.tsv",
       sep = "\t", row.names = TRUE)

cat("\n===== Cluster × SeuratClusters =====\n"); print(tab_sc)
cat("\n===== Cluster × Genotype =====\n");       print(tab_gt)
cat("\n===== GeneScore medians (marker × cluster) =====\n")
print(round(med_by_cluster, 3))

cat("\n완료. outputs/plots/ 와 outputs/tables/ 검토하여 CM_CLUSTERS 결정:\n")
cat(" - TNNT2/MYH6 중앙값 높음\n")
cat(" - COL1A1/COL3A1/EPCAM 낮음\n")
cat(" - SeuratClusters의 CM 라벨과 교차 일치\n")
