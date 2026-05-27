# Created: 2026-04-20 22:45
# DAR: pseudobulk → MACS2 peak calling → PeakMatrix → DAR (Wilcoxon) → volcano + motif enrichment
#
# DAR contrasts:
#   HOM vs WT  (homozygous TBX20-KO signature)
#   HET vs WT  (heterozygous / patient-state)
#
# Prereq: 06_cm_inject_leiden.R complete (projCM$Clusters and ClusterByGenotype injected)

suppressPackageStartupMessages({
  library(ArchR)
  library(data.table)
  library(ggplot2)
  library(SummarizedExperiment)
  library(BSgenome.Hsapiens.UCSC.hg38)   # required for addKmerBiasToCoverage (eval(parse(text=genome)))
})

WORK_DIR <- Sys.getenv("ARCHR_WORK_DIR", "/path/to/archr_project")
setwd(WORK_DIR)
addArchRThreads(threads = 8)
addArchRGenome("hg38")

OUTDIR <- "outputs/dar_0420"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

MACS2 <- "/etc/code/MACS2-2.2.7.1/bin/macs2"
stopifnot(file.exists(MACS2))

projCM <- loadArchRProject("ArchRSubset_CM", showLogo = FALSE)
cat("Loaded projCM:", length(projCM$cellNames), "cells\n")
cat("ClusterByGenotype levels:\n")
print(sort(table(projCM$ClusterByGenotype), decreasing = TRUE))

# ---------- 1) Pseudobulk coverage per ClusterByGenotype ----------
projCM <- addGroupCoverages(
  ArchRProj       = projCM,
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
projCM <- addReproduciblePeakSet(
  ArchRProj        = projCM,
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
projCM <- addPeakMatrix(projCM, force = TRUE)
saveArchRProject(projCM)

cat("\n[peaks] peak set size:\n")
print(getPeakSet(projCM))

# ---------- 3) DAR groups ----------
WT_BG  <- c("C1_x_wt", "C2_x_wt")
HET_FG <- c("C1_x_heterozygous", "C2_x_heterozygous")
HOM_FG <- "C3_x_homozygous"

run_dar <- function(useGroups, bgdGroups, label){
  cat(sprintf("\n[DAR] %s : useGroups=%s | bgdGroups=%s\n",
              label, paste(useGroups, collapse=","), paste(bgdGroups, collapse=",")))
  mt <- getMarkerFeatures(
    ArchRProj   = projCM,
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

dar_hom <- run_dar(HOM_FG, WT_BG, "HOM_vs_WT")   # primary
dar_het <- run_dar(HET_FG, WT_BG, "HET_vs_WT")   # clinical

# ---------- 4) Volcano plots ----------
pdf(file.path(OUTDIR, "DAR_volcano.pdf"), width = 7, height = 6)
tryCatch(
  print(plotMarkers(dar_hom$mt, name = HOM_FG, plotAs = "Volcano",
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
projCM <- addMotifAnnotations(projCM, motifSet = "cisbp", name = "Motif", force = TRUE)

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
  cat("  top TFs:\n")
  print(head(df, 10))
}

enrich <- function(mt, label, dir){
  cat(sprintf("\n[motif] %s %s\n", label, dir))
  cut_ <- if (dir == "up") "FDR <= 0.1 & Log2FC >=  0.25" else "FDR <= 0.1 & Log2FC <= -0.25"
  tryCatch({
    se <- peakAnnoEnrichment(mt, projCM, peakAnnotation = "Motif", cutOff = cut_)
    save_motif(se, file.path(OUTDIR, sprintf("motif_%s_%s.tsv", label, dir)))
  }, error = function(e) cat("  motif enrichment failed:", conditionMessage(e), "\n"))
}

enrich(dar_hom$mt, "HOM_vs_WT", "up")
enrich(dar_hom$mt, "HOM_vs_WT", "down")
enrich(dar_het$mt, "HET_vs_WT", "up")
enrich(dar_het$mt, "HET_vs_WT", "down")

saveArchRProject(projCM)

cat("\n===== Summary =====\n")
cat(sprintf("HOM_vs_WT DAR: %d up / %d down\n", dar_hom$n_up, dar_hom$n_dn))
cat(sprintf("HET_vs_WT DAR: %d up / %d down\n", dar_het$n_up, dar_het$n_dn))
cat("Outputs:", OUTDIR, "\n")
cat("DONE.\n")
