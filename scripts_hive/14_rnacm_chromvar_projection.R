# Created: 2026-04-24 10:30
# Updated: 2026-04-24 17:05
# ==============================================================================
# Phase 11 — Het chromVAR projection onto WT ↔ HOM motif-deviation axis.
#
# Motivation: rnaCM HET DAR = 0 peaks; at chromatin-accessibility level there is
# no detectable Het signature. But per-cell motif deviation scores (chromVAR)
# could still capture an intermediate Het position along the WT→HOM axis in TF
# activity space, particularly for TALE/GATA/NKX motifs identified in §6.
#
# Method:
#   1. Add MotifDeviations matrix to ArchRSubset_rnaCM (if not present)
#   2. Compute WT_FG centroid (mean deviation z-score across WT cells) and
#      HOM_FG centroid across the same ClusterByGenotype groups used in §5.
#   3. Axis vector = HOM_centroid - WT_centroid (normalised)
#   4. Per-cell projection score = (cell_deviation - WT_centroid) · axis_unit
#      → 0 = WT-like, 1 = HOM-like, <0 = beyond-WT, >1 = beyond-HOM
#   5. Export per-cell scores + per-motif axis loadings.
#
# Outputs (→ download to outputs/chromvar_proj_rnaCM_0424/):
#   chromvar_deviation_z_matrix.tsv.gz  — motifs × cells z-scores (optional, large)
#   chromvar_cell_projection.tsv        — per-cell projection score + metadata
#   chromvar_axis_loadings.tsv          — per-motif contribution to WT→HOM axis
#   chromvar_projection_summary.tsv     — group means (WT/HET/HOM × timepoint)
# ==============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(data.table)
  library(SummarizedExperiment)
  library(chromVAR)
  library(Matrix)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

.libPaths(.libPaths()[!grepl("/wynton/home/.*/R/x86_64", .libPaths())])

WORK_DIR <- "/gladstone/theodoris/lab/bkim/multi_multi/archr_dar_rnaCM"
setwd(WORK_DIR)
addArchRThreads(threads = 8)
addArchRGenome("hg38")

OUTDIR <- "outputs/chromvar_proj_rnaCM"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

projSub <- loadArchRProject("ArchRSubset_rnaCM", showLogo = FALSE)
cat("Loaded ArchRSubset_rnaCM:", length(projSub$cellNames), "cells\n")
cat("ClusterByGenotype table:\n")
print(sort(table(projSub$ClusterByGenotype), decreasing = TRUE))

# ---------- 1. Motif deviations ----------
has_dev <- "MotifMatrix" %in% getAvailableMatrices(projSub)
if (!has_dev) {
  cat("[deviations] MotifMatrix not found — computing with addDeviationsMatrix()…\n")
  # Motif annotations should already exist from Phase 6; add if missing
  if (!("Motif" %in% names(projSub@peakAnnotation))) {
    cat("[deviations] adding cisbp motif annotations first\n")
    projSub <- addMotifAnnotations(projSub, motifSet = "cisbp", name = "Motif",
                                   force = TRUE)
  }
  projSub <- addBgdPeaks(projSub, force = TRUE)
  projSub <- addDeviationsMatrix(projSub, peakAnnotation = "Motif",
                                 matrixName = "MotifMatrix", force = TRUE)
  saveArchRProject(projSub)
} else {
  cat("[deviations] MotifMatrix already present — reusing\n")
}

dev_se <- getMatrixFromProject(projSub, useMatrix = "MotifMatrix")
cat("[deviations] matrix: ", nrow(dev_se), " motifs × ", ncol(dev_se), " cells\n")

# use z-score assay (handles peak-count bias)
stopifnot("z" %in% assayNames(dev_se))
Z <- assays(dev_se)$z  # rows = motifs, cols = cells

# ---------- 2. Define WT / HOM centroids along Cluster × Genotype ----------
# multiseq_group/sample aren't reliably in cellColData; derive genotype from CxG
cell_meta <- data.table(
  cellName = colnames(Z),
  CxG      = projSub$ClusterByGenotype[match(colnames(Z), projSub$cellNames)],
  sample   = projSub$Sample[match(colnames(Z), projSub$cellNames)]
)
cell_meta[, raw_geno := sub("^C[0-9]+_x_", "", CxG)]
cell_meta[, genotype := fifelse(raw_geno == "wt", "WT",
                        fifelse(raw_geno == "heterozygous", "HET",
                        fifelse(raw_geno == "homozygous", "HOM", NA_character_)))]
cat("cell_meta genotype table (from CxG parse):\n")
print(table(cell_meta$genotype, useNA = "ifany"))
cat("cell_meta sample table:\n")
print(table(cell_meta$sample, useNA = "ifany"))

# use same grouping logic as Phase 6: WT = C1/C3 × wt, HOM = C4 × hom
WT_groups  <- c("C1_x_wt",           "C3_x_wt")
HET_groups <- c("C1_x_heterozygous", "C3_x_heterozygous")
HOM_groups <- c("C4_x_homozygous")

wt_cells  <- cell_meta$cellName[cell_meta$CxG %in% WT_groups]
het_cells <- cell_meta$cellName[cell_meta$CxG %in% HET_groups]
hom_cells <- cell_meta$cellName[cell_meta$CxG %in% HOM_groups]
cat(sprintf("WT n=%d | HET n=%d | HOM n=%d\n",
            length(wt_cells), length(het_cells), length(hom_cells)))

wt_centroid  <- rowMeans(as.matrix(Z[, wt_cells,  drop = FALSE]))
hom_centroid <- rowMeans(as.matrix(Z[, hom_cells, drop = FALSE]))

axis_vec <- hom_centroid - wt_centroid
axis_norm <- sqrt(sum(axis_vec^2))
axis_unit <- axis_vec / axis_norm
cat(sprintf("axis ||HOM-WT|| = %.3f\n", axis_norm))

# ---------- 3. Per-cell projection score ----------
# s_i = (z_i - wt_centroid) · axis_unit / axis_norm
#   → 0 = WT mean, 1 = HOM mean
center_Z <- sweep(as.matrix(Z), 1, wt_centroid, "-")
scores <- as.numeric(t(axis_unit / axis_norm) %*% center_Z)
names(scores) <- colnames(Z)

proj <- data.table(
  cellName   = colnames(Z),
  CxG        = cell_meta$CxG,
  genotype   = cell_meta$genotype,
  sample     = cell_meta$sample,
  projection = scores
)
fwrite(proj, file.path(OUTDIR, "chromvar_cell_projection.tsv"), sep = "\t")
cat("saved: chromvar_cell_projection.tsv  (", nrow(proj), " rows)\n")

# ---------- 4. Group-level summary ----------
summ <- proj[!is.na(genotype),
             .(n_cells = .N,
               proj_mean = mean(projection),
               proj_median = as.numeric(median(projection)),
               proj_sd = sd(projection)),
             by = .(genotype, sample)]
summ[, timepoint := paste0("D", sub("^[a-z]+", "", sample))]   # drops prefix letters (wt/ho/ht/w)
fwrite(summ, file.path(OUTDIR, "chromvar_projection_summary.tsv"), sep = "\t")
cat("saved: chromvar_projection_summary.tsv\n")

# ---------- 5. Axis loadings per motif ----------
loadings <- data.table(
  TF           = rownames(Z),
  wt_mean_z    = wt_centroid,
  hom_mean_z   = hom_centroid,
  delta        = axis_vec,
  abs_delta    = abs(axis_vec),
  axis_loading = axis_unit
)
setorder(loadings, -abs_delta)
fwrite(loadings, file.path(OUTDIR, "chromvar_axis_loadings.tsv"), sep = "\t")
cat("saved: chromvar_axis_loadings.tsv (top 10):\n")
print(head(loadings, 10))

cat("\n[done] chromVAR projection outputs in", OUTDIR, "\n")
