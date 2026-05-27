# Created: 2026-04-17 14:09
# ==============================================================================
# ArchR DAR pipeline — TBX20 dosage (WT / Het / Hom)
# ------------------------------------------------------------------------------
# Run on HPC cluster (SGE). Workflow adapted from Kathiriya 2026 ArchR pipeline.
#
# Inputs (placed alongside this script or in WORK_DIR):
#   valid_barcodes_G1.txt, _G2.txt, _G3.txt   (ATAC barcodes, doublet-free)
#   cell_metadata_for_archr.csv                (archr_cell_id,arc_library,
#                                               multiseq_group,multiseq_sample,...)
#   Fragment files: FRAG_DIR/G{1,2,3}/atac_fragments.tsv.gz
#
# Two-stage run (same script, edit STAGE + cluster vars between runs):
#   STAGE = "A"  → all-cell QC/LSI/Harmony/UMAP + CM marker feature plots
#                  ⇒ inspect plots, pick CM_CLUSTERS
#   STAGE = "B"  → CM subset re-cluster + genotype×cluster table
#                  ⇒ inspect, pick HET_FG/HOM_FG/WT_BG
#   STAGE = "C"  → pseudobulk + peaks + DAR + motif + browser tracks
# ==============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(parallel)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(ggplot2)
})

# ---- USER-EDITABLE CONFIG ----------------------------------------------------
STAGE       <- "A"   # "A" | "B" | "C"
WORK_DIR    <- "/gladstone/theodoris/lab/bkim/multi_multi/archr_dar"
FRAG_DIR    <- Sys.getenv("CELLRANGER_DIR", "/path/to/cellranger_output")
N_THREADS   <- 8

# CM cluster IDs (fill after STAGE A inspection of marker feature plots)
CM_CLUSTERS <- NULL    # e.g. c("C3","C4","C7")

# DAR group IDs (fill after STAGE B inspection of ClusterByGenotype cross-tab)
HET_FG <- NULL         # e.g. "C4_x_Het"
HOM_FG <- NULL         # e.g. "C3_x_Hom"
WT_BG  <- NULL         # e.g. c("C1_x_WT","C2_x_WT")
# -----------------------------------------------------------------------------

stopifnot(STAGE %in% c("A","B","C"))
setwd(WORK_DIR)
dir.create("outputs/plots",  recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/dar_0420",    recursive = TRUE, showWarnings = FALSE)

addArchRThreads(threads = N_THREADS)
addArchRGenome("hg38")

# ------------------------------------------------------------------------------
# STAGE A — all-cell Arrow → project → LSI → Harmony → clusters → CM markers
# ------------------------------------------------------------------------------
if (STAGE == "A") {

  sampleNames <- c("G1","G2","G3")
  inputFiles  <- setNames(
    file.path(FRAG_DIR, sampleNames, "atac_fragments.tsv.gz"),
    sampleNames)
  stopifnot(all(file.exists(inputFiles)))

  validBCs <- setNames(
    lapply(sampleNames, function(s) readLines(paste0("valid_barcodes_", s, ".txt"))),
    sampleNames)
  message("Valid barcodes per sample: ",
          paste(sampleNames, sapply(validBCs, length), sep="=", collapse=", "))

  # Arrow files
  if (!all(file.exists(paste0(sampleNames, ".arrow")))) {
    ArrowFiles <- createArrowFiles(
      inputFiles      = inputFiles,
      sampleNames     = sampleNames,
      validBarcodes   = validBCs,
      minTSS          = 4,
      minFrags        = 1000,
      addTileMat      = TRUE,
      addGeneScoreMat = TRUE,
      force           = FALSE)
  } else {
    ArrowFiles <- paste0(sampleNames, ".arrow")
    message("Reusing existing Arrow files.")
  }

  # Project
  proj <- ArchRProject(
    ArrowFiles        = ArrowFiles,
    outputDirectory   = "ArchRProject_all",
    copyArrows        = TRUE)

  # Metadata mapping — archr_cell_id matches proj$cellNames exactly ("G1#<bc>")
  meta <- read.csv("cell_metadata_for_archr.csv", stringsAsFactors = FALSE)
  idx  <- match(proj$cellNames, meta$archr_cell_id)
  stopifnot(all(!is.na(idx)))   # every ArchR cell must be in meta

  # Canonical Genotype labels: wt / Het / Hom (match reference "_x_" format)
  geno_map <- c(wt = "WT", heterozygous = "Het", homozygous = "Hom")
  proj$Genotype       <- factor(
    unname(geno_map[meta$multiseq_group[idx]]),
    levels = c("WT","Het","Hom"))
  proj$Library        <- meta$arc_library[idx]        # G1/G2/G3 (batch)
  proj$MultiseqSample <- meta$multiseq_sample[idx]    # e.g. htD6
  stopifnot(!any(is.na(proj$Genotype)))

  write.csv(
    as.data.frame.matrix(table(proj$Genotype, proj$Library)),
    "outputs/tables/all_genotype_by_library.csv")

  # Iterative LSI + Harmony (batch = Library) + clusters + UMAP
  proj <- addIterativeLSI(
    ArchRProj     = proj,
    useMatrix     = "TileMatrix",
    name          = "IterativeLSI",
    iterations    = 8,
    varFeatures   = 25000,
    dimsToUse     = 1:30,
    seed          = 1,
    clusterParams = list(resolution = c(0.1,0.2,0.4,0.8,1.0),
                         sampleCells = 10000, n.start = 10),
    force         = TRUE)

  proj <- addHarmony(proj, reducedDims = "IterativeLSI",
                     name = "Harmony", groupBy = "Library", force = TRUE)

  proj <- addClusters(proj, reducedDims = "Harmony", method = "Seurat",
                      name = "Clusters", resolution = 0.4, seed = 1, force = TRUE)

  proj <- addUMAP(proj, reducedDims = "Harmony", name = "UMAP",
                  nNeighbors = 20, minDist = 0.2, metric = "cosine",
                  seed = 1, force = TRUE)

  # Gene-score imputation for marker feature plots
  proj <- addImputeWeights(proj)

  # UMAPs: clusters / library / genotype
  geno_pal <- c(WT = "#0072B2", Het = "#E69F00", Hom = "#D55E00")  # project rule
  pC <- plotEmbedding(proj, colorBy = "cellColData", name = "Clusters",
                      embedding = "UMAP")
  pL <- plotEmbedding(proj, colorBy = "cellColData", name = "Library",
                      embedding = "UMAP")
  pG <- plotEmbedding(proj, colorBy = "cellColData", name = "Genotype",
                      embedding = "UMAP", pal = geno_pal)
  plotPDF(pC, pL, pG,
          name   = "StageA_UMAP_Clusters_Library_Genotype.pdf",
          ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

  # CM marker gene score feature plots
  cmMarkers <- c("TNNT2","MYH6","MYH7","TTN","NKX2-5","ACTN2","MYL7")
  pM <- plotEmbedding(proj, colorBy = "GeneScoreMatrix", name = cmMarkers,
                      embedding = "UMAP", rastr = FALSE, plotAs = "points",
                      quantCut = c(0.01, 0.95),
                      imputeWeights = getImputeWeights(proj))
  plotPDF(plotList = pM,
          name   = "StageA_GeneScore_CMmarkers.pdf",
          ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

  # Cluster × Genotype cross-tabs
  tab_cg <- as.data.frame.matrix(table(proj$Clusters, proj$Genotype))
  tab_cl <- as.data.frame.matrix(table(proj$Clusters, proj$Library))
  write.csv(tab_cg, "outputs/tables/StageA_cluster_x_genotype.csv")
  write.csv(tab_cl, "outputs/tables/StageA_cluster_x_library.csv")

  saveArchRProject(proj, outputDirectory = "ArchRProject_all", load = FALSE)

  message("STAGE A complete. Inspect PDFs under ArchRProject_all/Plots/",
          " and tables under outputs/tables/.",
          " Set CM_CLUSTERS and rerun with STAGE=\"B\".")
}

# ------------------------------------------------------------------------------
# STAGE B — CM subset → re-cluster → ClusterByGenotype cross-tab
# ------------------------------------------------------------------------------
if (STAGE == "B") {

  if (is.null(CM_CLUSTERS)) {
    stop("STAGE B requires CM_CLUSTERS. Edit the top of this script.")
  }

  proj <- loadArchRProject("ArchRProject_all", force = FALSE)
  cellsCM <- proj$cellNames[proj$Clusters %in% CM_CLUSTERS]
  message("CM cells selected: ", length(cellsCM),
          " (clusters: ", paste(CM_CLUSTERS, collapse=","), ")")

  projCM <- subsetArchRProject(
    ArchRProj       = proj,
    cells           = cellsCM,
    outputDirectory = "ArchRSubset_CM",
    dropCells       = TRUE,
    force           = TRUE)

  projCM <- addIterativeLSI(projCM, useMatrix = "TileMatrix",
                            name = "IterativeLSI",
                            iterations = 8, varFeatures = 25000,
                            dimsToUse = 1:30, seed = 1,
                            clusterParams = list(resolution = c(0.1,0.2,0.4,0.8,1.0),
                                                 sampleCells = 10000, n.start = 10),
                            force = TRUE)
  projCM <- addHarmony(projCM, reducedDims = "IterativeLSI",
                       name = "Harmony", groupBy = "Library", force = TRUE)
  projCM <- addClusters(projCM, reducedDims = "Harmony", method = "Seurat",
                        name = "Clusters", resolution = 0.4, seed = 1, force = TRUE)
  projCM <- addUMAP(projCM, reducedDims = "Harmony", name = "UMAP",
                    nNeighbors = 20, minDist = 0.2, metric = "cosine",
                    seed = 1, force = TRUE)

  # ClusterByGenotype: "C#_x_WT/Het/Hom" — matches reference pseudobulk labels
  projCM$ClusterByGenotype <- paste0(projCM$Clusters, "_x_",
                                     as.character(projCM$Genotype))

  projCM <- addImputeWeights(projCM)

  geno_pal <- c(WT = "#0072B2", Het = "#E69F00", Hom = "#D55E00")
  pC <- plotEmbedding(projCM, colorBy = "cellColData", name = "Clusters",
                      embedding = "UMAP")
  pL <- plotEmbedding(projCM, colorBy = "cellColData", name = "Library",
                      embedding = "UMAP")
  pG <- plotEmbedding(projCM, colorBy = "cellColData", name = "Genotype",
                      embedding = "UMAP", pal = geno_pal)
  plotPDF(pC, pL, pG,
          name   = "StageB_UMAP_CMsubset.pdf",
          ArchRProj = projCM, addDOC = FALSE, width = 5, height = 5)

  # cross-tabs
  tab_cg <- as.data.frame.matrix(table(projCM$Clusters, projCM$Genotype))
  tab_cbg_size <- as.data.frame(table(projCM$ClusterByGenotype))
  colnames(tab_cbg_size) <- c("ClusterByGenotype","nCells")
  tab_cl <- as.data.frame.matrix(table(projCM$ClusterByGenotype, projCM$Library))

  write.csv(tab_cg,       "outputs/tables/StageB_cluster_x_genotype.csv")
  write.csv(tab_cbg_size, "outputs/tables/StageB_ClusterByGenotype_sizes.csv",
            row.names = FALSE)
  write.csv(tab_cl,       "outputs/tables/StageB_ClusterByGenotype_by_library.csv")

  saveArchRProject(projCM, outputDirectory = "ArchRSubset_CM", load = FALSE)

  message("STAGE B complete. Pick HET_FG/HOM_FG/WT_BG from tables and rerun STAGE=\"C\".")
}

# ------------------------------------------------------------------------------
# STAGE C — pseudobulk, peaks, DAR (Het/Hom vs WT), motif, browser tracks
# ------------------------------------------------------------------------------
if (STAGE == "C") {

  if (is.null(HET_FG) || is.null(HOM_FG) || is.null(WT_BG)) {
    stop("STAGE C requires HET_FG, HOM_FG, WT_BG. Edit the top of this script.")
  }

  projCM <- loadArchRProject("ArchRSubset_CM", force = FALSE)

  # MACS2 path (module-provided on Wynton)
  pathToMacs2 <- unname(Sys.which("macs2"))
  if (!nzchar(pathToMacs2)) stop("macs2 not on PATH — load module before Rscript.")
  message("MACS2: ", pathToMacs2)

  # Pseudobulks — sample-aware via Library (G1/G2/G3 = 3 biological reps)
  projCM <- addGroupCoverages(
    ArchRProj     = projCM,
    groupBy       = "ClusterByGenotype",
    useLabels     = TRUE,
    minCells      = 500,
    maxCells      = 1000,
    maxFragments  = 25e6,
    minReplicates = 2,
    maxReplicates = 2,
    sampleRatio   = 0.8,
    force         = TRUE)

  # Reproducible peaks (2/2)
  projCM <- addReproduciblePeakSet(
    ArchRProj       = projCM,
    groupBy         = "ClusterByGenotype",
    reproducibility = "2",
    peaksPerCell    = 500,
    maxPeaks        = 150000,
    minCells        = 20,
    excludeChr      = c("chrM","chrY","chrX"),
    pathToMacs2     = pathToMacs2,
    shift           = -75,
    extsize         = 150,
    cutOff          = 0.1,
    extendSummits   = 250,
    promoterRegion  = c(2000, 100),
    plot            = TRUE,
    force           = TRUE)

  projCM <- addPeakMatrix(projCM)

  # -- DAR: Het vs WT ----------------------------------------------------------
  markerTest_Het <- getMarkerFeatures(
    ArchRProj   = projCM,
    useMatrix   = "PeakMatrix",
    groupBy     = "ClusterByGenotype",
    testMethod  = "wilcoxon",
    bias        = c("TSSEnrichment","log10(nFrags)"),
    useGroups   = HET_FG,
    bgdGroups   = WT_BG)

  up_Het   <- getMarkers(markerTest_Het, cutOff = "FDR <= 0.1 & Log2FC >= 0.25",  returnGR = TRUE)
  down_Het <- getMarkers(markerTest_Het, cutOff = "FDR <= 0.1 & Log2FC <= -0.25", returnGR = TRUE)
  write.csv(as.data.frame(up_Het[[HET_FG]]),
            sprintf("outputs/dar_0420/DAR_%s_vs_%s_up.csv",
                    HET_FG, paste(WT_BG, collapse="+")),
            row.names = FALSE)
  write.csv(as.data.frame(down_Het[[HET_FG]]),
            sprintf("outputs/dar_0420/DAR_%s_vs_%s_down.csv",
                    HET_FG, paste(WT_BG, collapse="+")),
            row.names = FALSE)

  vHet <- plotMarkers(seMarker = markerTest_Het, name = HET_FG,
                      cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.25", plotAs = "Volcano")
  plotPDF(vHet, name = sprintf("StageC_Volcano_%s_vs_WT.pdf", HET_FG),
          ArchRProj = projCM, addDOC = FALSE, width = 5, height = 5)

  # -- DAR: Hom vs WT ----------------------------------------------------------
  markerTest_Hom <- getMarkerFeatures(
    ArchRProj   = projCM,
    useMatrix   = "PeakMatrix",
    groupBy     = "ClusterByGenotype",
    testMethod  = "wilcoxon",
    bias        = c("TSSEnrichment","log10(nFrags)"),
    useGroups   = HOM_FG,
    bgdGroups   = WT_BG)

  up_Hom   <- getMarkers(markerTest_Hom, cutOff = "FDR <= 0.1 & Log2FC >= 0.25",  returnGR = TRUE)
  down_Hom <- getMarkers(markerTest_Hom, cutOff = "FDR <= 0.1 & Log2FC <= -0.25", returnGR = TRUE)
  write.csv(as.data.frame(up_Hom[[HOM_FG]]),
            sprintf("outputs/dar_0420/DAR_%s_vs_%s_up.csv",
                    HOM_FG, paste(WT_BG, collapse="+")),
            row.names = FALSE)
  write.csv(as.data.frame(down_Hom[[HOM_FG]]),
            sprintf("outputs/dar_0420/DAR_%s_vs_%s_down.csv",
                    HOM_FG, paste(WT_BG, collapse="+")),
            row.names = FALSE)

  vHom <- plotMarkers(seMarker = markerTest_Hom, name = HOM_FG,
                      cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.25", plotAs = "Volcano")
  plotPDF(vHom, name = sprintf("StageC_Volcano_%s_vs_WT.pdf", HOM_FG),
          ArchRProj = projCM, addDOC = FALSE, width = 5, height = 5)

  # -- Overlap Het ∩ Hom -------------------------------------------------------
  upHet_gr  <- up_Het[[HET_FG]];   downHet_gr <- down_Het[[HET_FG]]
  upHom_gr  <- up_Hom[[HOM_FG]];   downHom_gr <- down_Hom[[HOM_FG]]
  ov_up   <- findOverlaps(upHet_gr,   upHom_gr,   ignore.strand = FALSE)
  ov_down <- findOverlaps(downHet_gr, downHom_gr, ignore.strand = FALSE)
  write.csv(as.data.frame(ov_up),
            "outputs/dar_0420/Overlap_HetUp_HomUp.csv",     row.names = FALSE)
  write.csv(as.data.frame(ov_down),
            "outputs/dar_0420/Overlap_HetDown_HomDown.csv", row.names = FALSE)

  # -- Motif enrichment (cisbp) ------------------------------------------------
  projCM <- addMotifAnnotations(projCM, motifSet = "cisbp", name = "Motif", force = TRUE)

  motif_export <- function(seMarker, tag) {
    mUp <- peakAnnoEnrichment(seMarker, projCM, peakAnnotation = "Motif",
                              cutOff = "FDR <= 0.05 & Log2FC >= 0.25")
    mDo <- peakAnnoEnrichment(seMarker, projCM, peakAnnotation = "Motif",
                              cutOff = "FDR <= 0.05 & Log2FC <= -0.25")
    dfU <- data.frame(TF = rownames(mUp), mlog10Padj = assay(mUp)[,1])
    dfD <- data.frame(TF = rownames(mDo), mlog10Padj = assay(mDo)[,1])
    dfU <- dfU[order(-dfU$mlog10Padj), ];  dfU$rank <- seq_len(nrow(dfU))
    dfD <- dfD[order(-dfD$mlog10Padj), ];  dfD$rank <- seq_len(nrow(dfD))
    write.csv(dfU, sprintf("outputs/dar_0420/Motifs_%s_up.csv",   tag), row.names = FALSE)
    write.csv(dfD, sprintf("outputs/dar_0420/Motifs_%s_down.csv", tag), row.names = FALSE)
    list(up = dfU, down = dfD)
  }
  motif_Het <- motif_export(markerTest_Het, HET_FG)
  motif_Hom <- motif_export(markerTest_Hom, HOM_FG)

  # -- Browser tracks (TBX5, NPPA, NR2F1, GJA5, MYL7) --------------------------
  track_genes <- c("TBX5","NPPA","NR2F1","GJA5","MYL7")
  useGroups_track <- c(WT_BG, HET_FG, HOM_FG)
  pTrack <- plotBrowserTrack(
    ArchRProj   = projCM,
    groupBy     = "ClusterByGenotype",
    geneSymbol  = track_genes,
    useGroups   = useGroups_track)
  plotPDF(plotList = pTrack,
          name   = "StageC_BrowserTracks.pdf",
          ArchRProj = projCM, addDOC = FALSE, width = 6, height = 6)

  # Group-level BigWig for external visualization
  getGroupBW(projCM, groupBy = "ClusterByGenotype",
             normMethod = "ReadsInTSS", verbose = TRUE)

  saveArchRProject(projCM, outputDirectory = "ArchRSubset_CM", load = FALSE)

  message("STAGE C complete. DAR / motif / tracks written under outputs/.")
}
