# Created: 2026-04-23 20:15
# Updated: 2026-04-23 20:35
# ==============================================================================
# Phase 10 — ArchR browser tracks with HOM-vs-WT DAR shading.
#
# Reference: downloads/kathiriya-rao-2025-main/Fig5_atrialday20_cmsubset_ArchR.R
#   L288-305 (Fig 5H/I — plotBrowserTrack + DAR shading)
#
# 6 loci (user-selected 2026-04-23):
#   TBX20, MEIS1, TGIF1, GATA2, ISL1, NKX2-5
#   - TBX20              : study locus (chr7:35.42Mb HOM-up peak)
#   - MEIS1, TGIF1       : HOM_DOWN TALE axis (maturation arrest)
#   - GATA2, ISL1        : HOM_UP progenitor re-opening (dedifferentiation)
#   - NKX2-5             : negative control (cardiac backbone should be unchanged)
#
# Groups shown (WT_BG + HET_FG + HOM_FG, set in 04_rnaCM_dar_peak_motif.R):
#   C1_x_wt, C3_x_wt,
#   C1_x_heterozygous, C3_x_heterozygous,
#   C4_x_homozygous
#
# Output: outputs/dar_rnaCM/browser_tracks/<GENE>.pdf
# ==============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(grid)
})

.libPaths(.libPaths()[!grepl("/wynton/home/.*/R/x86_64", .libPaths())])

WORK_DIR <- "/gladstone/theodoris/lab/bkim/multi_multi/archr_dar_rnaCM"
setwd(WORK_DIR)
addArchRThreads(threads = 8)
addArchRGenome("hg38")

OUTDIR <- "outputs/dar_rnaCM/browser_tracks"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---------- 1) Load ArchR project and DAR marker SE ----------
projSub <- loadArchRProject("ArchRSubset_rnaCM", showLogo = FALSE)
cat("Loaded ArchRSubset_rnaCM:", length(projSub$cellNames), "cells\n")

mt_path <- "outputs/dar_rnaCM/markerTest_HOM_vs_WT.rds"
stopifnot(file.exists(mt_path))
mt_hom <- readRDS(mt_path)

# Significant DAR GRanges (FDR ≤ 0.1, |Log2FC| ≥ 0.25) — same cutoff as DAR table
sig_gr_list <- getMarkers(
  mt_hom,
  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.25",
  returnGR = TRUE
)
# getMarkers returns a SimpleList keyed by useGroups entry; HOM_FG = C4_x_homozygous
sig_name <- names(sig_gr_list)[1]
dar_gr   <- sig_gr_list[[sig_name]]
cat(sprintf("HOM_vs_WT sig peaks (%s): %d\n", sig_name, length(dar_gr)))

# Up / down split for two-color shading
dar_up <- dar_gr[mcols(dar_gr)$Log2FC >=  0.25]
dar_dn <- dar_gr[mcols(dar_gr)$Log2FC <= -0.25]
cat(sprintf("  up: %d  down: %d  (actual 500 bp peaks)\n",
            length(dar_up), length(dar_dn)))

# Visual widen: at ±500 kb zoom a 500 bp peak is <0.05% of plot width and
# renders as sub-pixel. Widen to ±5 kb around each summit (10 kb total) so the
# shaded bars are clearly visible. Actual peak coordinates remain in DAR TSVs.
widen_for_display <- function(gr, pad = 5000) {
  mid <- start(gr) + (width(gr) %/% 2L)
  GRanges(seqnames(gr),
          IRanges(start = pmax(1L, mid - pad), end = mid + pad),
          strand  = strand(gr),
          mcols(gr))
}
dar_up_disp <- widen_for_display(dar_up)
dar_dn_disp <- widen_for_display(dar_dn)

features_gr <- GRangesList(
  DAR_up   = dar_up_disp,
  DAR_down = dar_dn_disp
)

# ---------- 2) Groups & genes ----------
USE_GROUPS <- c(
  "C1_x_wt",
  "C3_x_wt",
  "C1_x_heterozygous",
  "C3_x_heterozygous",
  "C4_x_homozygous"
)
cat("useGroups:", paste(USE_GROUPS, collapse = ", "), "\n")

# NKX2-5 symbol in hg38/GENCODE is "NKX2-5"; ArchR getGenes honors gene symbols.
GENES <- c("TBX20", "MEIS1", "TGIF1", "GATA2", "ISL1", "NKX2-5")

# ---------- 3) plotBrowserTrack ----------
cat("\n[browser] running plotBrowserTrack for", length(GENES), "genes...\n")

# Window: ±500 kb. Rationale — DAR peak density analysis shows that for 4 of 6
# loci (TBX20, MEIS1, TGIF1, GATA2) the nearest HOM-vs-WT DAR peak is 70-490 kb
# from the TSS. 500 kb ensures DAR shading is visible at every gene panel.
p_list <- plotBrowserTrack(
  ArchRProj = projSub,
  groupBy   = "ClusterByGenotype",
  useGroups = USE_GROUPS,
  geneSymbol = GENES,
  features   = features_gr,
  upstream   = 500000,
  downstream = 500000,
  loops      = getPeak2GeneLinks(projSub)  # harmless if NULL; ArchR tolerates missing links
)

# ---------- 4) Write one PDF per gene ----------
sanitize <- function(g) gsub("[^A-Za-z0-9._-]", "_", g)

for (g in GENES) {
  p <- p_list[[g]]
  if (is.null(p)) {
    cat(sprintf("  %-8s  SKIP (no plot returned)\n", g))
    next
  }
  out_path <- file.path(OUTDIR, sprintf("%s.pdf", sanitize(g)))
  pdf(out_path, width = 7, height = 9)
  grid::grid.newpage()
  grid::grid.draw(p)
  dev.off()
  cat(sprintf("  %-8s  → %s\n", g, out_path))
}

# Also save a combined multi-page PDF for easy flip-through
combined <- file.path(OUTDIR, "all_6_loci.pdf")
pdf(combined, width = 7, height = 9)
for (g in GENES) {
  p <- p_list[[g]]
  if (is.null(p)) next
  grid::grid.newpage()
  grid::grid.draw(p)
}
dev.off()
cat(sprintf("combined → %s\n", combined))

cat("\n===== Phase 10 summary =====\n")
cat(sprintf("  sig DAR peaks: %d up / %d down\n", length(dar_up), length(dar_dn)))
cat("  groups:", length(USE_GROUPS), "\n")
cat("  output dir:", OUTDIR, "\n")
cat("DONE.\n")
