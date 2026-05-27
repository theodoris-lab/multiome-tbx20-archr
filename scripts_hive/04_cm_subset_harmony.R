# Created: 2026-04-20 14:45
# ==============================================================================
# CM subset + LSI + Harmony + UMAP, export Harmony embedding
# ------------------------------------------------------------------------------
# Prereq: review 03_cm_identify.R outputs → decide CM_CLUSTERS
# Output: ArchRSubset_CM/ (subsetted project)
#         ArchRSubset_CM/harmony_CM.tsv  (cells × 30 Harmony dims, for scanpy)
#         outputs/plots/CM_subset_UMAP_before_leiden.pdf
# Next:   run 05_leiden_cm.py
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

# ----------------------------------------------------------------------
# USER: fill in CM clusters after reviewing 03_cm_identify.R outputs
CM_CLUSTERS <- c("")   # e.g. c("C3","C4","C6","C7") — must be filled in
# ----------------------------------------------------------------------
stopifnot(length(CM_CLUSTERS) > 0 && all(nzchar(CM_CLUSTERS)))

proj <- loadArchRProject(path = file.path(WORK_DIR, "ArchRProject"))

cellsCM <- proj$cellNames[proj$Clusters %in% CM_CLUSTERS]
cat("CM cells:", length(cellsCM),
    "/", length(proj$cellNames),
    sprintf("(%.1f%%)", 100 * length(cellsCM) / length(proj$cellNames)), "\n")

projCM <- subsetArchRProject(
  ArchRProj       = proj,
  cells           = cellsCM,
  outputDirectory = file.path(WORK_DIR, "ArchRSubset_CM"),
  dropCells       = TRUE,
  force           = TRUE
)

projCM <- addIterativeLSI(
  ArchRProj  = projCM, useMatrix = "TileMatrix", name = "IterativeLSI",
  iterations = 4, varFeatures = 25000, dimsToUse = 1:30,
  clusterParams = list(resolution = c(0.2,0.4,0.6), sampleCells = 10000, n.start = 10),
  seed = 1, force = TRUE
)

projCM <- addHarmony(
  ArchRProj   = projCM, reducedDims = "IterativeLSI", name = "Harmony",
  groupBy     = "Sample", force = TRUE
)

projCM <- addUMAP(
  ArchRProj   = projCM, reducedDims = "Harmony", name = "UMAP",
  nNeighbors  = 30, minDist = 0.3, metric = "cosine",
  seed = 1, force = TRUE
)

# Harmony embedding + cellName export (for scanpy leiden)
harmony <- getReducedDims(projCM, reducedDims = "Harmony")
out_path <- file.path(WORK_DIR, "ArchRSubset_CM", "harmony_CM.tsv")
write.table(harmony, out_path, sep = "\t", quote = FALSE, col.names = NA)
cat("Harmony exported:", out_path, " (", nrow(harmony), "cells ×", ncol(harmony), "dims)\n")

# Overview plots before leiden
pg <- plotEmbedding(projCM, colorBy = "cellColData", name = "Genotype",
                    embedding = "UMAP", plotAs = "points", rastr = FALSE,
                    pal = c(Wild = "blue", Het = "orange", Hom = "red"))
ps <- plotEmbedding(projCM, colorBy = "cellColData", name = "Sample",
                    embedding = "UMAP", plotAs = "points", rastr = FALSE)
plotPDF(pg, ps, name = "CM_subset_UMAP_before_leiden.pdf",
        ArchRProj = projCM, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(projCM)
cat("ArchRSubset_CM saved. Next: run 05_leiden_cm.py\n")
