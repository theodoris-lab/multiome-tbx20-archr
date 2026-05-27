# Created: 2026-04-20 14:45
# ==============================================================================
# Inject scanpy leiden labels into ArchRSubset_CM
# ------------------------------------------------------------------------------
# Prereq: 05_leiden_cm.py complete → ArchRSubset_CM/leiden_CM.tsv
# Output: ArchRSubset_CM/ with Clusters + ClusterByGenotype saved
#         outputs/plots/CM_subset_UMAP_leiden.pdf
#         outputs/tables/Leiden_vs_Genotype.tsv
# Next:   07_dar_peak_motif.R (submit as qsub batch)
# ==============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(data.table)
})

.libPaths(.libPaths()[!grepl("/wynton/home/.*/R/x86_64", .libPaths())])

WORK_DIR <- Sys.getenv("ARCHR_WORK_DIR", "/path/to/archr_project")
setwd(WORK_DIR)

addArchRThreads(threads = 8)
addArchRGenome("hg38")

projCM <- loadArchRProject(path = file.path(WORK_DIR, "ArchRSubset_CM"))

lc <- fread(file.path(WORK_DIR, "ArchRSubset_CM", "leiden_CM.tsv"),
            data.table = FALSE)
rownames(lc) <- lc$cellName

# Align order + apply ArchR-style labels (C1, C2, …)
ord <- match(projCM$cellNames, rownames(lc))
stopifnot(!any(is.na(ord)))
leiden_int <- as.integer(as.character(lc[ord, "leiden_res04"]))
projCM$Clusters          <- paste0("C", leiden_int + 1L)   # 0-based → 1-based
projCM$ClusterByGenotype <- paste0(projCM$Clusters, "_x_", projCM$Genotype)

# Also inject adjacent resolutions for reference
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

cat("\nDone. Review Leiden_vs_Genotype.tsv to select DAR groups:\n")
cat(" - WT-dominant cluster(s) → bgdGroups\n")
cat(" - Het-skewed cluster     → HET_FG (useGroups)\n")
cat(" - Hom-skewed cluster     → HOM_FG\n")
cat("Need >= 500 cells per cluster for stable addGroupCoverages.\n")
