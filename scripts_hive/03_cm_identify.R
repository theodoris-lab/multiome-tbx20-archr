# Created: 2026-04-20 14:45
# ==============================================================================
# CM cluster identification on all-cells ArchRProject
# ------------------------------------------------------------------------------
# Prereq: 02_archr_dar_analysis.R Stage A complete (LSI + Harmony + Clusters + UMAP saved)
# Output: outputs/plots/CM_markers_UMAP.pdf
#         outputs/tables/Cluster_vs_SeuratClusters.tsv
#         outputs/tables/Cluster_vs_Genotype.tsv
#         outputs/tables/GeneScore_median_by_cluster.tsv
# Next:   review plots, decide CM_CLUSTERS → record in 04_cm_subset_harmony.R
# ==============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(ggplot2)
  library(data.table)
})

# Prevent mixing with host R library paths
.libPaths(.libPaths()[!grepl("/wynton/home/.*/R/x86_64", .libPaths())])

WORK_DIR <- Sys.getenv("ARCHR_WORK_DIR", "/path/to/archr_project")
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
# CM markers + non-CM markers (for exclusion)
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

# Clusters UMAP (for label verification)
pc <- plotEmbedding(proj, colorBy = "cellColData", name = "Clusters",
                    embedding = "UMAP", plotAs = "points", rastr = FALSE)
pg <- plotEmbedding(proj, colorBy = "cellColData", name = "Genotype",
                    embedding = "UMAP", plotAs = "points", rastr = FALSE,
                    pal = c(Wild = "blue", Het = "orange", Hom = "red"))
ps <- plotEmbedding(proj, colorBy = "cellColData", name = "SeuratClusters",
                    embedding = "UMAP", plotAs = "points", rastr = FALSE)
plotPDF(pc, pg, ps, name = "UMAP_overview.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# ---- 3) Gene score median by cluster ------------------------------------------
# GeneScoreMatrix → per-cluster median → identify top clusters by TNNT2 etc.
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

# ---- 4) Cluster × SeuratClusters, Cluster × Genotype cross-tabulation --------
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

cat("\nDone. Review outputs/plots/ and outputs/tables/ to decide CM_CLUSTERS:\n")
cat(" - high TNNT2/MYH6 median score\n")
cat(" - low COL1A1/COL3A1/EPCAM\n")
cat(" - consistent with CM labels in SeuratClusters cross-tab\n")
