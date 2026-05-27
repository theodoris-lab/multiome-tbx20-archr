# Created: 2026-04-23 13:25
# Updated: 2026-04-23 13:25
# ==============================================================================
# IterativeLSI + Harmony + UMAP within rnaCM subset.
#            Export Harmony embedding (TSV) for local scanpy Leiden.
# Prereq: 09_rnacm_subset.R complete (ArchRSubset_rnaCM/ exists)
# Output: ArchRSubset_rnaCM/ (LSI + Harmony + UMAP added)
#         ArchRSubset_rnaCM/harmony_rnaCM.tsv
#         outputs/plots/rnaCM_UMAP_before_leiden.pdf
# Next:   download harmony_rnaCM.tsv, run scripts_local/03_leiden_rnacm.py
# ==============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(data.table)
})

.libPaths(.libPaths()[!grepl("/wynton/home/.*/R/x86_64", .libPaths())])

WORK_DIR <- "/gladstone/theodoris/lab/bkim/multi_multi/archr_dar_rnaCM"
setwd(WORK_DIR)

addArchRThreads(threads = 8)
addArchRGenome("hg38")

projSub <- loadArchRProject(file.path(WORK_DIR, "ArchRSubset_rnaCM"), showLogo = FALSE)
cat("Loaded ArchRSubset_rnaCM:", length(projSub$cellNames), "cells\n")

projSub <- addIterativeLSI(
  ArchRProj  = projSub, useMatrix = "TileMatrix", name = "IterativeLSI",
  iterations = 4, varFeatures = 25000, dimsToUse = 1:30,
  clusterParams = list(resolution = c(0.2, 0.4, 0.6), sampleCells = 10000, n.start = 10),
  seed = 1, force = TRUE
)

projSub <- addHarmony(
  ArchRProj   = projSub, reducedDims = "IterativeLSI", name = "Harmony",
  groupBy     = "Sample", force = TRUE
)

projSub <- addUMAP(
  ArchRProj   = projSub, reducedDims = "Harmony", name = "UMAP",
  nNeighbors  = 30, minDist = 0.3, metric = "cosine",
  seed = 1, force = TRUE
)

harmony <- getReducedDims(projSub, reducedDims = "Harmony")
out_path <- file.path(WORK_DIR, "ArchRSubset_rnaCM", "harmony_rnaCM.tsv")
write.table(harmony, out_path, sep = "\t", quote = FALSE, col.names = NA)
cat("Harmony exported:", out_path, " (", nrow(harmony), "cells ×", ncol(harmony), "dims)\n")

pg <- plotEmbedding(projSub, colorBy = "cellColData", name = "Genotype",
                    embedding = "UMAP", plotAs = "points", rastr = FALSE,
                    pal = c(Wild = "blue", Het = "orange", Hom = "red"))
ps <- plotEmbedding(projSub, colorBy = "cellColData", name = "Sample",
                    embedding = "UMAP", plotAs = "points", rastr = FALSE)
plotPDF(pg, ps, name = "rnaCM_UMAP_before_leiden.pdf",
        ArchRProj = projSub, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(projSub)
cat("DONE. Next: download harmony_rnaCM.tsv and run scripts_local/03_leiden_rnacm.py\n")
