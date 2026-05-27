# Created: 2026-04-20 23:45
# TBX20 locus browser track + quick locus QC (knockout effect check)
#
# Goal: visualize chromatin accessibility at TBX20 gene body and flanking region (±500kb) across
#       WT_BG (C1/C2 wt), HET_FG (C1/C2 het), HOM_FG (C3 hom)
# Output: outputs/dar_0420/TBX20_browser.pdf + overlapping DAR summary

suppressPackageStartupMessages({
  library(ArchR)
  library(data.table)
  library(ggplot2)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

WORK_DIR <- Sys.getenv("ARCHR_WORK_DIR", "/path/to/archr_project")
setwd(WORK_DIR)
addArchRThreads(threads = 4)
addArchRGenome("hg38")

OUTDIR <- "outputs/dar_0420"
projCM <- loadArchRProject("ArchRSubset_CM", showLogo = FALSE)

WT_BG  <- c("C1_x_wt", "C2_x_wt")
HET_FG <- c("C1_x_heterozygous", "C2_x_heterozygous")
HOM_FG <- "C3_x_homozygous"

target_groups <- c(WT_BG, HET_FG, HOM_FG)

# ---- 1) Browser track at TBX20 ----
# ArchR infers range from geneSymbol via gene annotation
pdf(file.path(OUTDIR, "TBX20_browser.pdf"), width = 8, height = 10)
for (ext in c(100000, 500000)){
  cat(sprintf("[browser] TBX20 ± %d bp\n", ext))
  tryCatch({
    pl <- plotBrowserTrack(
      ArchRProj  = projCM,
      groupBy    = "ClusterByGenotype",
      useGroups  = target_groups,
      geneSymbol = "TBX20",
      upstream   = ext,
      downstream = ext,
      loops      = NULL
    )
    grid::grid.newpage()
    grid::grid.draw(pl$TBX20)
  }, error = function(e) cat("  failed:", conditionMessage(e), "\n"))
}
dev.off()

# ---- 2) Peak overlap summary at TBX20 locus ----
ps <- getPeakSet(projCM)
tbx20_window <- GRanges("chr7", IRanges(34800000, 35800000))   # ±500kb
peaks_in_window <- subsetByOverlaps(ps, tbx20_window)
cat(sprintf("\n[QC] Total peaks in chr7:34.8-35.8Mb: %d\n", length(peaks_in_window)))

# Check DAR overlap
dar_hom_up <- fread(file.path(OUTDIR, "DAR_HOM_vs_WT_up.tsv"))
dar_hom_dn <- fread(file.path(OUTDIR, "DAR_HOM_vs_WT_down.tsv"))

locus_dar_up <- dar_hom_up[seqnames == "chr7" & start > 34800000 & end < 35800000]
locus_dar_dn <- dar_hom_dn[seqnames == "chr7" & start > 34800000 & end < 35800000]
cat(sprintf("HOM DAR up in TBX20 ±500kb: %d\n", nrow(locus_dar_up)))
cat(sprintf("HOM DAR down in TBX20 ±500kb: %d\n", nrow(locus_dar_dn)))
cat("\nDAR up at TBX20 locus:\n"); print(locus_dar_up)
cat("\nDAR down at TBX20 locus:\n"); print(locus_dar_dn)

# ---- 3) BigWig export (TBX20 region only, for IGV) ----
tryCatch({
  getGroupBW(
    ArchRProj  = projCM,
    groupBy    = "ClusterByGenotype",
    normMethod = "ReadsInTSS",
    tileSize   = 100,
    maxCells   = 2000,
    ceiling    = 4,
    threads    = 4
  )
  cat("\n[bw] group bigwig files exported to GroupBigWigs/\n")
}, error = function(e) cat("[bw] failed:", conditionMessage(e), "\n"))

cat("\nDONE.\n")
