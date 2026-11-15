#!/usr/bin/env Rscript
# GORG-Dark: assembly-derived 150bp fragments / 30bp overlap for 31-mers.
# Prepare UMAP inputs: CountSketch, Kraken taxon counts, informative 31-mer counts.
# No debarcoding, alignment, host filtering, read QC, or assembly workflow.
# GORG_ROOT/GORG_RUNNER/GORG_THREADS/GORG_MEM/GORG_KRAKEN_DB override defaults.
scriptFile <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
if(is.null(scriptFile)) {
  arg <- grep("^--file=", commandArgs(), value=TRUE)
  if(!length(arg)) stop("Use Rscript prep_gorg.R or source('prep_gorg.R')")
  scriptFile <- sub("^--file=", "", arg[1])
}
scriptDir <- dirname(normalizePath(scriptFile, mustWork=TRUE))
source(file.path(scriptDir, "prep_dna_common.R"))
config <- InitializeDNAWorkflow("GORG", "/husky/henriksson/atrandi/v9_wgs_gorg")
bascetRoot <- config$root

################################################################################
############################ Prepare sorted TIRP ###############################
################################################################################
inputName <- "filtered_nohost" # Existing requested basename; no host filtering.
if(!length(detectShardsForFile(bascetRoot, inputName))) {
  status <- system2("env", c(
    shQuote(paste0("GORG_DIR=", bascetRoot)),
    "bash", shQuote(file.path(scriptDir, "run_gorg_sag2tirp.sh"))
  ))
  if(status!=0) stop("GORG SAG-to-TIRP conversion failed")
}

################################################################################
############################ DNA feature workflows #############################
################################################################################
RunDNAFeatures(config, inputName)
message("GORG count/sketch inputs ready for downstream UMAP analysis: ", bascetRoot)
