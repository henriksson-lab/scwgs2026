#!/usr/bin/env Rscript
# Ling/Aarestrup: import barcode-A-split triplets, shard, prepare UMAP inputs.
# Methods: CountSketch, Kraken taxon counts, and informative 31-mer counts.
# No reference alignment, host filtering, read QC, or assembly workflow.
# LING_ROOT/LING_RUNNER/LING_THREADS/LING_MEM/LING_KRAKEN_DB override defaults.
# Cell QC follows prep_mock.R (95th read-count percentile); LING_KEEP_ALL=1
# retains every imported barcode instead.
scriptFile <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
if(is.null(scriptFile)) {
  arg <- grep("^--file=", commandArgs(), value=TRUE)
  if(!length(arg)) stop("Use Rscript prep_ling.R or source('prep_ling.R')")
  scriptFile <- sub("^--file=", "", arg[1])
}
scriptDir <- dirname(normalizePath(scriptFile, mustWork=TRUE))
source(file.path(scriptDir, "prep_dna_common.R"))
config <- InitializeDNAWorkflow("LING", "/husky/henriksson/atrandi/v9_ling")
bascetRoot <- config$root

################################################################################
############################ Import and shard ##################################
################################################################################
# The deposit already separates biological reads and the barcode FASTQ.
# The importer combines measured data by matched read IDs, without inventing
# spacers/qualities or pretending the files are raw Atrandi four-round reads.
if(!file.exists(file.path(bascetRoot, "debarcoded.meta"))) {
  status <- system2("bash", c(
    shQuote(file.path(scriptDir, "run_ling_import.sh")),
    "--root", shQuote(bascetRoot)
  ))
  if(status!=0) stop("Ling triplet import failed; ensure source downloads are complete")
}
minQuantile <- as.numeric(config$env("MIN_QUANTILE", "0.95"))
stopifnot(!is.na(minQuantile), minQuantile>=0, minQuantile<=1)
debstat <- PrepareSharding(
  bascetRoot, inputName="debarcoded", minQuantile=minQuantile,
  bascetInstance=config$instance
)
# Import always preserves every observation; this optional override also
# retains all barcodes in the downstream DNA analyses.
if(config$env("KEEP_ALL", "0")=="1") {
  debstat$list_picked_cells <- lapply(debstat$list_group_bc_count, function(x) x$cellid)
  debstat$all_knee$picked <- TRUE
}
message("Retaining ", sum(lengths(debstat$list_picked_cells)), " barcode combinations")
saveRDS(debstat, file.path(bascetRoot, "ling_sharding.rds"))
DebarcodedKneePlot(debstat, filename=file.path(bascetRoot, "kneeplot.pdf"))
WaitForJob(BascetShardify(
  debstat, numOutputShards=1, outputName="filtered",
  runner=config$runner, bascetInstance=config$instance
))

################################################################################
############################ DNA feature workflows #############################
################################################################################
RunDNAFeatures(config, "filtered")
message("Ling count/sketch inputs ready for downstream UMAP analysis: ", bascetRoot)
