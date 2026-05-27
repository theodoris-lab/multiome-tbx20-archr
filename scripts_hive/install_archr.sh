#!/bin/bash
# Created: 2026-04-17 14:31
#
# One-shot ArchR install on Wynton dev node (or compute via qrsh).
# Run inside tmux so SSH drops don't abort the 30-60min compile:
#
#   ssh dev1
#   tmux new -s archr_install
#   bash /gladstone/theodoris/lab/bkim/multi_multi/archr_dar/scripts/install_archr.sh
#   # Ctrl+b, d  → detach
#   # Later: tmux attach -t archr_install   OR   tmux kill-session -t archr_install
#
set -euo pipefail

module load CBI
module load r/4.4.3

# User-writable R library (system lib is read-only on Wynton).
# R 4.4.x uses x86_64-pc-linux-gnu-library/4.4 by convention.
USER_LIB="${HOME}/R/x86_64-pc-linux-gnu-library/4.4"
mkdir -p "$USER_LIB"
export R_LIBS_USER="$USER_LIB"

echo "===> R version + libPaths"
Rscript -e 'cat(R.version.string, "\n"); cat(.libPaths(), sep="\n")'

echo "===> Installing ArchR and its Bioc dependencies (this takes a while)"
Rscript -e '
  repos <- "https://cloud.r-project.org"
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos=repos)
  if (!requireNamespace("remotes",    quietly=TRUE)) install.packages("remotes",     repos=repos)

  BiocManager::install(
    c("chromVAR","motifmatchr","presto",
      "BSgenome.Hsapiens.UCSC.hg38","JASPAR2020","TFBSTools"),
    ask = FALSE, update = FALSE)

  # ArchR is github-only (not on CRAN/Bioc). Use Bioc + CRAN repos for deps.
  remotes::install_github("GreenleafLab/ArchR",
                          ref = "master",
                          repos = BiocManager::repositories(),
                          upgrade = "never")

  suppressMessages(library(ArchR))
  cat("\n>>> ArchR core:", as.character(packageVersion("ArchR")), "\n")

  ArchR::installExtraPackages()
  cat(">>> ArchR::installExtraPackages() OK\n")
'

echo "===> Final verification"
Rscript -e '
  suppressMessages({
    library(ArchR); library(chromVAR); library(motifmatchr)
    library(BSgenome.Hsapiens.UCSC.hg38)
  })
  cat("ArchR:", as.character(packageVersion("ArchR")), "\n")
  cat("R   :", R.version.string, "\n")
'
echo "===> DONE."
