# Created: 2026-04-23 13:25
# Updated: 2026-04-23 17:00
# ==============================================================================
# Phase 6 — DAR + motif enrichment on ArchRSubset_rnaCM.
# Mirrors 06_dar_peak_motif.R but writes to outputs/dar_rnaCM/.
#
# DAR axes:
#   Axis 1 (primary)  : HOM_FG vs WT_BG — amplified LVNC signature
#   Axis 2 (clinical) : HET_FG vs WT_BG — patient-state
#
# Cluster IDs below MUST be reviewed after 03_rnaCM_inject_leiden.R produces
# Leiden_vs_Genotype_rnaCM.tsv. Edit WT_BG / HET_FG / HOM_FG accordingly.
# ==============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(data.table)
  library(ggplot2)
  library(SummarizedExperiment)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

.libPaths(.libPaths()[!grepl("/wynton/home/.*/R/x86_64", .libPaths())])

WORK_DIR <- "/gladstone/theodoris/lab/bkim/multi_multi/archr_dar_rnaCM"
setwd(WORK_DIR)
addArchRThreads(threads = 8)
addArchRGenome("hg38")

OUTDIR <- "outputs/dar_rnaCM"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

MACS2 <- "/etc/code/MACS2-2.2.7.1/bin/macs2"
stopifnot(file.exists(MACS2))

projSub <- loadArchRProject("ArchRSubset_rnaCM", showLogo = FALSE)
cat("Loaded ArchRSubset_rnaCM:", length(projSub$cellNames), "cells\n")
cat("ClusterByGenotype levels:\n")
print(sort(table(projSub$ClusterByGenotype), decreasing = TRUE))

# ---------- 1) Pseudobulk coverage per ClusterByGenotype ----------
projSub <- addGroupCoverages(
  ArchRProj       = projSub,
  groupBy         = "ClusterByGenotype",
  useLabels       = TRUE,
  minCells        = 500,
  maxCells        = 1000,
  minReplicates   = 2,
  maxReplicates   = 2,
  sampleRatio     = 0.8,
  force           = TRUE
)
cat("\n[coverages] done\n")

# ---------- 2) Reproducible peak set via MACS2 ----------
projSub <- addReproduciblePeakSet(
  ArchRProj        = projSub,
  groupBy          = "ClusterByGenotype",
  pathToMacs2      = MACS2,
  reproducibility  = "2",
  peaksPerCell     = 500,
  maxPeaks         = 150000,
  excludeChr       = c("chrM", "chrY", "chrX"),
  shift            = -75,
  extsize          = 150,
  cutOff           = 0.1,
  extendSummits    = 250,
  promoterRegion   = c(2000, 100),
  force            = TRUE
)
projSub <- addPeakMatrix(projSub, force = TRUE)
saveArchRProject(projSub)

cat("\n[peaks] peak set size:\n")
print(getPeakSet(projSub))

# ---------- 3) DAR groups ----------
# Set per Leiden_vs_Genotype_rnaCM.tsv (2026-04-23 17:00):
#   C1 : Het 2360 / Hom   88 /  WT 1681  → WT/Het baseline (no Hom)
#   C2 : Het 1510 / Hom 1816 /  WT  762  → mixed, excluded
#   C3 : Het 1872 / Hom   97 /  WT 1745  → WT/Het baseline (no Hom)
#   C4 : Het  433 / Hom 2463 /  WT  262  → Hom-dominant (78%) → HOM signature
#   C5/C6 : mixed / tiny → excluded
WT_BG  <- c("C1_x_wt",            "C3_x_wt")              # 1,681 + 1,745 = 3,426
HET_FG <- c("C1_x_heterozygous",  "C3_x_heterozygous")    # 2,360 + 1,872 = 4,232
HOM_FG <- c("C4_x_homozygous")                            # 2,463

run_dar <- function(useGroups, bgdGroups, label){
  cat(sprintf("\n[DAR] %s : useGroups=%s | bgdGroups=%s\n",
              label, paste(useGroups, collapse=","), paste(bgdGroups, collapse=",")))
  mt <- getMarkerFeatures(
    ArchRProj   = projSub,
    useMatrix   = "PeakMatrix",
    groupBy     = "ClusterByGenotype",
    testMethod  = "wilcoxon",
    bias        = c("TSSEnrichment", "log10(nFrags)"),
    useGroups   = useGroups,
    bgdGroups   = bgdGroups
  )
  saveRDS(mt, file.path(OUTDIR, sprintf("markerTest_%s.rds", label)))
  up <- getMarkers(mt, cutOff = "FDR <= 0.1 & Log2FC >=  0.25", returnGR = TRUE)
  dn <- getMarkers(mt, cutOff = "FDR <= 0.1 & Log2FC <= -0.25", returnGR = TRUE)
  n_up <- if (length(up)) length(up[[1]]) else 0L
  n_dn <- if (length(dn)) length(dn[[1]]) else 0L
  cat(sprintf("  up: %d  down: %d\n", n_up, n_dn))

  if (n_up > 0) fwrite(as.data.frame(up[[1]]),
                       file.path(OUTDIR, sprintf("DAR_%s_up.tsv", label)), sep = "\t")
  if (n_dn > 0) fwrite(as.data.frame(dn[[1]]),
                       file.path(OUTDIR, sprintf("DAR_%s_down.tsv", label)), sep = "\t")
  list(mt = mt, up = up, dn = dn, n_up = n_up, n_dn = n_dn)
}

dar_hom <- run_dar(HOM_FG, WT_BG, "HOM_vs_WT")
dar_het <- run_dar(HET_FG, WT_BG, "HET_vs_WT")

# ---------- 4) Volcano plots ----------
pdf(file.path(OUTDIR, "DAR_volcano.pdf"), width = 7, height = 6)
tryCatch(
  print(plotMarkers(dar_hom$mt, name = HOM_FG[1], plotAs = "Volcano",
                    cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.25")),
  error = function(e) cat("volcano HOM failed:", conditionMessage(e), "\n")
)
for (g in HET_FG){
  tryCatch(
    print(plotMarkers(dar_het$mt, name = g, plotAs = "Volcano",
                      cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.25")),
    error = function(e) cat("volcano HET", g, "failed:", conditionMessage(e), "\n")
  )
}
dev.off()

# ---------- 5) Motif enrichment (cisbp) ----------
projSub <- addMotifAnnotations(projSub, motifSet = "cisbp", name = "Motif", force = TRUE)

save_motif <- function(se, path){
  if (is.null(se)) { cat("  (no motif SE)\n"); return(invisible()) }
  m <- assays(se)
  df <- data.frame(
    TF         = rownames(se),
    mlog10Padj = m[["mlog10Padj"]][, 1],
    mlog10p    = m[["mlog10p"]][, 1],
    Enrichment = m[["Enrichment"]][, 1]
  )
  df <- df[order(-df$mlog10Padj), ]
  fwrite(df, path, sep = "\t")
  cat("  top TFs:\n"); print(head(df, 10))
}

enrich <- function(mt, label, dir){
  cat(sprintf("\n[motif] %s %s\n", label, dir))
  cut_ <- if (dir == "up") "FDR <= 0.1 & Log2FC >=  0.25" else "FDR <= 0.1 & Log2FC <= -0.25"
  tryCatch({
    se <- peakAnnoEnrichment(mt, projSub, peakAnnotation = "Motif", cutOff = cut_)
    save_motif(se, file.path(OUTDIR, sprintf("motif_%s_%s.tsv", label, dir)))
  }, error = function(e) cat("  motif enrichment failed:", conditionMessage(e), "\n")
  )
}

enrich(dar_hom$mt, "HOM_vs_WT", "up")
enrich(dar_hom$mt, "HOM_vs_WT", "down")
enrich(dar_het$mt, "HET_vs_WT", "up")
enrich(dar_het$mt, "HET_vs_WT", "down")

saveArchRProject(projSub)

cat("\n===== Summary (rnaCM) =====\n")
cat(sprintf("HOM_vs_WT DAR: %d up / %d down\n", dar_hom$n_up, dar_hom$n_dn))
cat(sprintf("HET_vs_WT DAR: %d up / %d down\n", dar_het$n_up, dar_het$n_dn))
cat("Outputs:", OUTDIR, "\n")
cat("DONE.\n")
