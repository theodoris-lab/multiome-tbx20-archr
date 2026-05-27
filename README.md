<!-- Created: 2026-05-26 17:20 -->
<!-- Updated: 2026-05-26 18:10 -->

# multiome-tbx20-archr

ATAC-seq analysis of TBX20 mutation iPSC-CMs — differentially accessible regions, motif enrichment, and chromVAR, using [ArchR](https://www.archrproject.com/).

## Background

TBX20 mutation is associated with left ventricular non-compaction (LVNC) cardiomyopathy. This analysis uses single-nucleus RNA + ATAC (10x Genomics Multiome) from iPSC-derived cardiomyocytes (iPSC-CMs) carrying a TBX20 mutation.

ArchR is used to identify differentially accessible regions (DARs), TF motif enrichment, and chromVAR deviation scores across genotypes (WT, Het, Hom).

Two parallel ATAC pipelines are included:
- **All-CM pipeline** (01–08): all cardiomyocyte clusters, ATAC peak-gated
- **rnaCM pipeline** (09–15): RNA-expression-gated CM subset with stricter identity criteria

## Data

- 10x Genomics Multiome kit (nuclear RNA + ATAC; nuclear only)
- Differentiation days D0 – D30
- Input: Seurat `.rds` object

## Repository Structure

```
.
├── scripts_hive/      # run on HPC (SGE)
└── scripts_local/     # run locally
```

### `scripts_hive/` — HPC scripts (SGE / apptainer)

**Environment**

| Script | What it does |
|---|---|
| `install_archr.sh` | Build ArchR apptainer container and install R packages |
| `run_dar_qsub.sh` | SGE job submission template (memory, slots, module loading) |

**All-CM pipeline**

| Script | What it does |
|---|---|
| `01_extract_barcodes.py` | Extract CM barcodes from Seurat UMAP clusters |
| `02_archr_dar_analysis.R` | ArchR project creation, peak calling, DAR (Wilcoxon) |
| `03_cm_identify.R` | Identify CM clusters within ArchR |
| `04_cm_subset_harmony.R` | Subset CMs, Harmony batch correction |
| `05_leiden_cm.py` | Leiden clustering (Python, leidenalg) |
| `06_cm_inject_leiden.R` | Inject Leiden cluster labels into ArchR object |
| `07_dar_peak_motif.R` | DAR per cluster + motif enrichment (chromVAR) |
| `08_tbx20_browser.R` | TBX20 locus ArchR browser track |

**rnaCM pipeline** (RNA-expression-gated CM subset)

| Script | What it does |
|---|---|
| `09_rnacm_subset.R` | Subset rnaCM cells from ArchR project |
| `10_rnacm_lsi_harmony.R` | LSI dimensionality reduction + Harmony |
| `11_rnacm_inject_leiden.R` | Inject Leiden labels into rnaCM ArchR object |
| `12_rnacm_dar_peak_motif.R` | DAR + motif enrichment for rnaCM |
| `13_rnacm_browser_tracks.R` | Browser tracks for rnaCM |
| `14_rnacm_chromvar_projection.R` | chromVAR deviations, WT → Hom projection axis |
| `15_rnacm_meis_validation.R` | MEIS binding-site Fisher test in Hom-DOWN DARs |

### `scripts_local/` — local scripts

| Script | What it does |
|---|---|
| `01_leiden_allcm.py` | Leiden clustering for all-CM (input for step 06) |
| `02_extract_rnacm_barcodes.py` | Extract rnaCM barcodes from RNA-gated cells |
| `03_leiden_rnacm.py` | Leiden clustering for rnaCM subset |
| `04_tbx20_locus_scatter.py` | TBX20 locus ATAC signal scatter |
| `04b_tbx20_het_extract.R` | Extract TBX20 Het allele read counts |
| `05_motif_rank_plot.py` | Motif enrichment rank plot |
| `06_rna_atac_concordance.py` | RNA–ATAC concordance at DAR loci |
| `07_chromvar_projection.py` | chromVAR WT→Hom axis projection plot |
| `08_meis_validation.py` | MEIS binding validation visualization |
| `09_dar_allcm_vs_rnacm.py` | Compare DAR results between all-CM and rnaCM |

## Dependencies

**R** (scripts_hive)
- `ArchR` (≥ 1.0.3), `Seurat`, `Signac`, `harmony`, `chromVAR`, `motifmatchr`, `BiocParallel`
- Run inside apptainer container — see `install_archr.sh`

**Python** (scripts_local + Leiden steps)
- `scanpy`, `leidenalg`, `igraph`, `anndata`, `matplotlib`, `pandas`, `numpy`

Tested with ArchR 1.0.3, R 4.4, Python 3.9.

## Notes

- This is analysis code for a specific dataset, not a general-purpose tool or package.
- Data paths are set at the top of each script — change them before running.
- Scripts are numbered by order but are not a fully automated pipeline. Check each step before proceeding.
- All-CM and rnaCM pipelines are independent from step 09 onward and can run in parallel after step 04.
- `run_dar_qsub.sh` shows the SGE job configuration used (8 cores, 8G per slot = 64G total, 8h wall time).
- `scripts_hive/` scripts are for Wynton HPC (SGE). Adapt as needed for other environments.

## Reference

ArchR: [Granja et al. 2021 *Nature Genetics*](https://doi.org/10.1038/s41588-021-00790-6)
