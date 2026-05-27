# Created: 2026-04-21 00:30
# Extract HET vs WT (and HOM vs WT) stats at 2 TBX20 downstream peaks
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(S4Vectors)
})

OUT <- "/Users/bkim/vscode/V2_multiome_2026-04-14/outputs/dar_0420"
het <- readRDS(file.path(OUT, "markerTest_HET_vs_WT.rds"))
hom <- readRDS(file.path(OUT, "markerTest_HOM_vs_WT.rds"))

peaks <- data.frame(
  seqnames = c("chr7", "chr7"),
  start    = c(35418901, 35888942),
  end      = c(35419401, 35889442)
)

pick <- function(se, chr, s, e){
  rd <- rowData(se)
  idx <- which(as.character(rd$seqnames) == chr & rd$start == s & rd$end == e)
  if (length(idx) == 0) return(NULL)
  data.frame(
    Log2FC   = assay(se, "Log2FC")[idx, 1],
    FDR      = assay(se, "FDR")[idx, 1],
    Pval     = assay(se, "Pval")[idx, 1],
    MeanDiff = assay(se, "MeanDiff")[idx, 1]
  )
}

cat("=== Peak 1: chr7:35418901-35419401 (TBX20 +125kb) ===\n")
cat("HOM vs WT:\n"); print(pick(hom, "chr7", 35418901, 35419401))
cat("HET vs WT:\n"); print(pick(het, "chr7", 35418901, 35419401))

cat("\n=== Peak 2: chr7:35888942-35889442 (TBX20 +595kb) ===\n")
cat("HOM vs WT:\n"); print(pick(hom, "chr7", 35888942, 35889442))
cat("HET vs WT:\n"); print(pick(het, "chr7", 35888942, 35889442))
