# Created: 2026-04-23 13:25
# Updated: 2026-04-23 16:55
# ==============================================================================
# Inject scanpy leiden labels into ArchRSubset_rnaCM, set up
#           Clusters + ClusterByGenotype. Mirrors 06_cm_inject_leiden.R logic.
# Prereq: upload leiden_rnaCM.tsv (from scripts_local/03_leiden_rnacm.py) to cluster.
# Input:  ArchRSubset_rnaCM/leiden_rnaCM.tsv  (row index = cellName, cols = leiden_res04, ...)
# Output: ArchRSubset_rnaCM/ with Clusters, ClusterByGenotype injected + UMAP plots
#         outputs/tables/Leiden_vs_Genotype_rnaCM.tsv
# Next:   12_rnacm_dar_peak_motif.R
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

# Note: additional leiden_res columns are not injected — ArchRProject does not support [[<-.
# Use addCellColData() if needed.

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

cat("\nDone. Review Leiden_vs_Genotype_rnaCM.tsv and update DAR groups in 12_rnacm_dar_peak_motif.R.\n")
cat("  - WT-dominant cluster(s) → bgdGroups\n")
cat("  - Het-skewed cluster     → HET_FG\n")
cat("  - Hom-skewed cluster     → HOM_FG\n")
cat("Cluster IDs may differ from the CM-subset run; check genotype proportions.\n")
