# Created: 2026-04-24 17:50
# Updated: 2026-04-24 18:05
# ==============================================================================
# Phase 12 — MEIS1 binding-site validation: per-TALE-TF enrichment in HOM-DOWN DAR
#
# Asks: within the TALE family, is MEIS specifically the dominant driver in
# HOM-DOWN DARs, or do TGIF/PKNOX contribute equally?
#
# Method:
#   1. Load ArchRSubset_rnaCM (already has Motif annotation)
#   2. Extract per-peak × motif binary match matrix via getMatches()
#   3. Load markerTest_HOM_vs_WT.rds for peak-level FDR/Log2FC
#   4. DOWN = FDR<=0.1 & Log2FC<=-1, BG = NS (FDR>0.1)
#   5. Per-TALE-TF Fisher enrichment (one-sided greater) in DOWN vs BG
#   6. Positive control: same for WT-side chromVAR top TFs (MEF2/TBX/MAFB)
#   7. Within-DOWN co-occurrence matrix across TALE TFs
#
# Outputs (-> outputs/meis_validation_rnaCM/):
#   tale_perTF_enrichment.tsv         per-motif Fisher result for TALE family
#   wt_side_perTF_enrichment.tsv      positive control (other WT-side TFs)
#   down_peak_tale_cooccurrence.tsv   TF × TF cooccurrence count in DOWN peaks
#   down_peak_tale_count_dist.tsv     # of TALE motifs per DOWN peak
# ==============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(data.table)
  library(SummarizedExperiment)
  library(Matrix)
})

.libPaths(.libPaths()[!grepl("/wynton/home/.*/R/x86_64", .libPaths())])

WORK_DIR <- "/gladstone/theodoris/lab/bkim/multi_multi/archr_dar_rnaCM"
setwd(WORK_DIR)
addArchRThreads(threads = 8)
addArchRGenome("hg38")

OUTDIR <- "outputs/meis_validation_rnaCM"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

projSub <- loadArchRProject("ArchRSubset_rnaCM", showLogo = FALSE)
cat("Loaded ArchRSubset_rnaCM:", length(projSub$cellNames), "cells\n")

# ---------- 1. Motif × peak match matrix ----------
matches_se <- getMatches(projSub, name = "Motif")
M <- assay(matches_se)   # sparse binary, rows = peaks, cols = motifs
cat("matches matrix:", nrow(M), "peaks x", ncol(M), "motifs\n")

# ---------- 2. Load markerTest for peak-level FDR/Log2FC ----------
mt_path <- "outputs/dar_rnaCM/markerTest_HOM_vs_WT.rds"
stopifnot(file.exists(mt_path))
mt <- readRDS(mt_path)
mt_rd <- rowData(mt)
cat("rowData columns:\n"); print(colnames(mt_rd))
mt_dt <- data.table(
  chr    = as.character(mt_rd$seqnames),
  start  = as.integer(mt_rd$start),
  end    = as.integer(mt_rd$end),
  Log2FC = as.numeric(assay(mt, "Log2FC")[, 1]),
  FDR    = as.numeric(assay(mt, "FDR")[, 1])
)
cat("markerTest peaks:", nrow(mt_dt), "\n")

# matches_se and markerTest both use the project peakset in identical order
stopifnot(nrow(M) == nrow(mt_dt))

# ---------- 3. DOWN vs background ----------
down_idx <- which(mt_dt$FDR <= 0.1 & mt_dt$Log2FC <= -1)
up_idx   <- which(mt_dt$FDR <= 0.1 & mt_dt$Log2FC >=  1)
ns_idx   <- which(mt_dt$FDR > 0.1)
cat(sprintf("DOWN n=%d | UP n=%d | NS n=%d\n",
            length(down_idx), length(up_idx), length(ns_idx)))

# ---------- 4. Per-TF Fisher enrichment ----------
fisher_per_tf <- function(motif_col, peak_idx_target, peak_idx_bg) {
  hits <- as.numeric(M[, motif_col])
  a <- sum(hits[peak_idx_target])
  b <- length(peak_idx_target) - a
  c <- sum(hits[peak_idx_bg])
  d <- length(peak_idx_bg) - c
  if (a + c == 0) {
    return(data.table(motif = motif_col, n_target_hit = a, n_target = a + b,
                      n_bg_hit = c, n_bg = c + d,
                      pct_target = NA_real_, pct_bg = NA_real_,
                      OR = NA_real_, pval = NA_real_))
  }
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
  data.table(motif = motif_col,
             n_target_hit = a, n_target = a + b,
             n_bg_hit = c, n_bg = c + d,
             pct_target = 100 * a / (a + b),
             pct_bg = 100 * c / (c + d),
             OR = unname(ft$estimate),
             pval = ft$p.value)
}

all_motifs <- colnames(M)

# TALE family motifs (cisbp naming: TF_NNN suffix)
tale_pat  <- "^(MEIS1|MEIS2|MEIS3|PKNOX1|PKNOX2|TGIF1|TGIF2|TGIF2LX|TGIF2LY|PBX1|PBX2|PBX3|PBX4|IRX1|IRX2|IRX3|IRX4|IRX5|IRX6)_"
tale_cols <- grep(tale_pat, all_motifs, value = TRUE)
cat("\nTALE motifs found in cisbp annotation:\n")
print(tale_cols)

tale_res <- rbindlist(lapply(tale_cols, fisher_per_tf,
                             peak_idx_target = down_idx,
                             peak_idx_bg = ns_idx))
tale_res[, fdr := p.adjust(pval, method = "BH")]
tale_res[, TF_clean := sub("_\\d+$", "", motif)]
setorder(tale_res, pval)
fwrite(tale_res, file.path(OUTDIR, "tale_perTF_enrichment.tsv"), sep = "\t")
cat("\n[TALE enrichment in HOM-DOWN vs NS] (sorted by pval)\n")
print(tale_res[, .(motif, pct_target, pct_bg, OR, pval, fdr)])

# Positive control: WT-side chromVAR top TFs
wt_pat  <- "^(MEF2A|MEF2B|MEF2C|MEF2D|TBX4|TBX5|MAFB|MAF|SNAI3|MGA|CIC)_"
wt_cols <- grep(wt_pat, all_motifs, value = TRUE)
wt_res <- rbindlist(lapply(wt_cols, fisher_per_tf,
                           peak_idx_target = down_idx,
                           peak_idx_bg = ns_idx))
wt_res[, fdr := p.adjust(pval, method = "BH")]
wt_res[, TF_clean := sub("_\\d+$", "", motif)]
setorder(wt_res, pval)
fwrite(wt_res, file.path(OUTDIR, "wt_side_perTF_enrichment.tsv"), sep = "\t")
cat("\n[WT-side positive control]\n")
print(wt_res[, .(motif, pct_target, pct_bg, OR, pval, fdr)])

# ---------- 5. Within-DOWN co-occurrence ----------
co_mat <- as.matrix(M[down_idx, tale_cols, drop = FALSE])
storage.mode(co_mat) <- "integer"
cooc <- crossprod(co_mat)   # TF × TF count of co-occurrence
cooc_dt <- as.data.table(cooc, keep.rownames = "TF_a")
fwrite(cooc_dt, file.path(OUTDIR, "down_peak_tale_cooccurrence.tsv"), sep = "\t")

peak_tale_count <- rowSums(co_mat)
count_dist <- as.data.table(table(n_TALE_motifs = peak_tale_count))
fwrite(count_dist, file.path(OUTDIR, "down_peak_tale_count_dist.tsv"), sep = "\t")
cat("\n[TALE motifs per DOWN peak]\n")
print(count_dist)

cat("\n[done] outputs in", OUTDIR, "\n")
