# Created: 2026-04-23 13:25
# Updated: 2026-04-23 13:25
# ==============================================================================
# Phase 2 of DAR_analysis_plan_rnaCM_0423 — subset existing ArchRProject by
# RNA-defined CM barcodes (scanpy leiden_res0.1 cluster 0, 17,345 cells).
#
# Prereq: scp outputs/dar_rnaCM_inputs_0423/rnaCM_barcodes.tsv to
#   /gladstone/theodoris/lab/bkim/multi_multi/archr_dar_rnaCM/inputs/
# Input:  /gladstone/theodoris/lab/bkim/multi_multi/archr_dar/ArchRProject/
# Output: /gladstone/theodoris/lab/bkim/multi_multi/archr_dar_rnaCM/ArchRSubset_rnaCM/
#         inputs_logs/subset_summary.txt
# Next:   02_rnaCM_lsi_harmony.R
# ==============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(data.table)
})

.libPaths(.libPaths()[!grepl("/wynton/home/.*/R/x86_64", .libPaths())])

SRC_WORK  <- "/gladstone/theodoris/lab/bkim/multi_multi/archr_dar"
WORK_DIR  <- "/gladstone/theodoris/lab/bkim/multi_multi/archr_dar_rnaCM"
INPUT_DIR <- file.path(WORK_DIR, "inputs")
LOG_DIR   <- file.path(WORK_DIR, "inputs_logs")
dir.create(WORK_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR,   showWarnings = FALSE, recursive = TRUE)
setwd(WORK_DIR)

addArchRThreads(threads = 8)
addArchRGenome("hg38")

proj <- loadArchRProject(path = file.path(SRC_WORK, "ArchRProject"), showLogo = FALSE)
cat("Source ArchRProject cells:", length(proj$cellNames), "\n")

bc_path <- file.path(INPUT_DIR, "rnaCM_barcodes.tsv")
stopifnot(file.exists(bc_path))
bc <- fread(bc_path, header = TRUE)$cellName
cat("RNA CM barcodes from AnnData:", length(bc), "\n")

# 포맷 sanity check: G{1,2,3}#<barcode>
stopifnot(all(grepl("^G[123]#", bc)))

keep_mask <- proj$cellNames %in% bc
keep_cells <- proj$cellNames[keep_mask]
n_keep <- length(keep_cells)
n_drop <- length(bc) - n_keep
cat(sprintf("rnaCM cells in ArchR: %d / %d  (QC-dropped: %d, %.2f%%)\n",
            n_keep, length(bc), n_drop, 100 * n_drop / length(bc)))

projSub <- subsetArchRProject(
  ArchRProj       = proj,
  cells           = keep_cells,
  outputDirectory = file.path(WORK_DIR, "ArchRSubset_rnaCM"),
  dropCells       = TRUE,
  force           = TRUE
)
saveArchRProject(projSub)

# Summary log
sink(file.path(LOG_DIR, "subset_summary.txt"))
cat("== rnaCM subset summary ==\n")
cat("source project cells : ", length(proj$cellNames), "\n")
cat("rnaCM barcodes (RNA) : ", length(bc), "\n")
cat("kept in ArchR subset : ", n_keep, "\n")
cat("QC dropped           : ", n_drop, sprintf(" (%.2f%%)\n", 100 * n_drop / length(bc)))
cat("\n== Genotype breakdown (rnaCM subset) ==\n")
print(table(projSub$Genotype))
cat("\n== Sample breakdown ==\n")
print(table(projSub$Sample))
sink()

cat("DONE. subsetted project at:", file.path(WORK_DIR, "ArchRSubset_rnaCM"), "\n")
cat("Next: 02_rnaCM_lsi_harmony.R\n")
