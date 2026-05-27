#!/bin/bash
# Created: 2026-04-17 14:29
#
# Wynton uses SGE (Sun Grid Engine), not Slurm.
#   Submit:   qsub scripts/run_dar.sh
#   Status:   qstat -u bkim
#   Delete:   qdel <job_id>
#
# Memory (mem_free) is PER SLOT (per-core) in SGE, so total = mem_free × smp.
# 8G × 8 slots = 64G total.
#
#$ -N archr_dar
#$ -cwd
#$ -pe smp 8
#$ -l mem_free=8G
#$ -l h_rt=08:00:00
#$ -o outputs/logs/
#$ -e outputs/logs/
#$ -j n

set -euo pipefail

WORK_DIR="/gladstone/theodoris/lab/bkim/multi_multi/archr_dar"
cd "${WORK_DIR}"
mkdir -p outputs/logs

module load CBI
module load r/4.4.3

# User-writable R library (ArchR and Bioc deps installed here).
export R_LIBS_USER="${HOME}/R/x86_64-pc-linux-gnu-library/4.4"

# macs2 installed via `python3 -m pip install --user` → ~/.local/bin
export PATH="${HOME}/.local/bin:${PATH}"

Rscript scripts/02_archr_dar_analysis.R
