#!/usr/bin/env Rscript
# Run after ./run_yanai_download.sh and the reference scripts complete.
# Environment overrides: YANAI_ROOT, YANAI_ZORN_DIR, YANAI_BASCET_DIR,
# YANAI_BASCET_TARGET, YANAI_RUNNER (local/slurm), YANAI_THREADS, YANAI_MEM.
# YANAI_ARRAY_WORKERS=1; threads/memory apply separately to each local task.
# All valid barcode combinations are retained; align RNA and count reads per gene.
# Gene-count matrices are inputs for downstream UMAP analysis; no QC/assembly stages.
# YANAI_GFF overrides the annotation next to YANAI_REFERENCE (all.gff by default).

scriptFile <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
if(is.null(scriptFile)) {
  arg <- grep("^--file=", commandArgs(), value=TRUE)
  if(!length(arg)) stop("Use Rscript prep_yanai.R or source('prep_yanai.R')")
  scriptFile <- sub("^--file=", "", arg[1])
}
scriptDir <- dirname(normalizePath(scriptFile, mustWork=TRUE))
source(file.path(scriptDir, "yanai_petriseq", "count_genes.R"))

bascetRoot <- normalizePath(
  Sys.getenv("YANAI_ROOT", "/husky/henriksson/atrandi/v9_rnaseq_yanai"),
  mustWork=TRUE
)
zornDir <- path.expand(Sys.getenv("YANAI_ZORN_DIR", "~/github/zorn"))

required <- c("stringr", "fs", "data.table", "readr", "ggplot2", "naturalsort", "processx")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) {
  stop("Missing R packages needed by Zorn: ", paste(missing, collapse=", "))
}

# Use the current checkout, like prep_mock.R, so the new PETRI-seq chemistry
# is available without requiring a new installed Zorn release.
for(file in c(
  "job_general.R", "job_local.R", "job_slurm.R", "jobscript.R",
  "bascet_file.R", "zorn.R", "shell.R", "zorn_aggr.R",
  "refgenome.R", "container.R", "debarcoding.R"
)) {
  source(file.path(zornDir, "R", file))
}
if(!"petriseq" %in% eval(formals(BascetDebarcode)$chemistry)) {
  stop("This Zorn checkout does not yet support chemistry='petriseq'")
}

setwd(bascetRoot)
manifest <- read.delim("fastq_manifest.tsv", stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(nrow(manifest)>0, !anyDuplicated(manifest$library_prefix))
rawmeta_dir <- readLines("rawdata.txt", warn=FALSE)
rawmeta_dir <- trimws(rawmeta_dir[nzchar(trimws(rawmeta_dir))])
if(length(rawmeta_dir)!=1L) stop("rawdata.txt must contain one BGZF FASTQ directory")

rawmeta <- DetectRawFileMeta(rawmeta_dir, verbose=TRUE)
if(nrow(rawmeta)==0) stop("No prepared FASTQs found; run run_yanai_download.sh first")
# DetectRawFileMeta deliberately omits the prefix when there is only one
# library. Keep it explicit here so subset and full runs use identical IDs.
rawmeta$prefix <- sub("_S[0-9].*$", "", rawmeta$r1)
row <- match(rawmeta$prefix, manifest$library_prefix)
if(anyNA(row)) stop("FASTQ library names do not match the GEO manifest")
if(!setequal(rawmeta$prefix, manifest$library_prefix)) {
  stop("Not all manifest libraries have prepared FASTQs; wait for the download/preparation to complete")
}
stopifnot(!anyDuplicated(rawmeta$prefix))
rawmeta$run_accession <- manifest$run_accession[row]
rawmeta$geo_accession <- manifest$geo_accession[row]
rawmeta$library_title <- manifest$library_title[row]
rawmeta$species <- manifest$scientific_name[row]
write.csv(rawmeta, "yanai_libraries.csv", row.names=FALSE)

is_beagle <- Sys.info()[["nodename"]]=="beagle"
runner_mode <- Sys.getenv("YANAI_RUNNER", if(is_beagle) "local" else "slurm")
numThreads <- as.integer(Sys.getenv("YANAI_THREADS", "10"))
stopifnot(!is.na(numThreads), numThreads>=6L)
default_bascet <- if(runner_mode=="local") "~/github/bascet" else "/home/m/mahogny/github/bascet"
bascetInstance.default <- getBascetDevDir(
  path.expand(Sys.getenv("YANAI_BASCET_DIR", default_bascet)),
  targetType=Sys.getenv("YANAI_BASCET_TARGET", "release")
)
if(runner_mode=="local") {
  bascetRunner.default <- LocalRunner(
    direct=TRUE, showScript=FALSE, ncpu=as.character(numThreads),
    mem=Sys.getenv("YANAI_MEM", "50g")
  )
  arrayWorkers <- as.integer(Sys.getenv("YANAI_ARRAY_WORKERS","1"))
  if(is.na(arrayWorkers)||arrayWorkers<1L) stop("YANAI_ARRAY_WORKERS must be positive")
  if(arrayWorkers>1L) {
    source(file.path(scriptDir,"workflow_local_parallel.R"))
    bascetRunner.default <- DatasetArrayRunner(bascetRunner.default,arrayWorkers)
  }
} else if(runner_mode=="slurm") {
  bascetRunner.default <- SlurmRunner(
    account=Sys.getenv("YANAI_SLURM_ACCOUNT", "hpc2n2026-186"),
    ncpu=as.character(numThreads), mem=Sys.getenv("YANAI_MEM", "40g"),
    logTime=TRUE, verbose=TRUE, time="0-14:00:00"
  )
} else {
  stop("YANAI_RUNNER must be local or slurm")
}

debarcode_job <- BascetDebarcode(
  bascetRoot, rawmeta,
  chemistry="petriseq",
  outputName="debarcoded",
  numThreads=numThreads,
  runner=bascetRunner.default,
  bascetInstance=bascetInstance.default
)
WaitForJob(debarcode_job)

# Build the normal Zorn library grouping/statistics, but retain ALL valid
# barcode combinations. minQuantile=0 alone would still drop minimum-count
# cells because PrepareSharding uses a strict greater-than comparison.
debstat <- PrepareSharding(
  bascetRoot, inputName="debarcoded", minQuantile=0,
  bascetInstance=bascetInstance.default
)
debstat$list_picked_cells <- lapply(debstat$list_group_bc_count, function(x) x$cellid)
debstat$all_knee$picked <- TRUE
message("Retaining all ", sum(lengths(debstat$list_picked_cells)),
        " valid barcode combinations across ", debstat$numgroup, " libraries")
saveRDS(debstat, "yanai_sharding.rds")
DebarcodedKneePlot(debstat, filename="kneeplot.pdf")

# One shard per independent library; Zorn merges only files with the same
# prefix, preserving cell identity. No WGS host filter or k-mer analysis.
shard_job <- BascetShardify(
  debstat, numOutputShards=1, outputName="filtered",
  runner=bascetRunner.default,
  bascetInstance=bascetInstance.default
)
WaitForJob(shard_job)
write.csv(data.frame(
  group=seq_along(debstat$group_prefix),
  library_prefix=debstat$group_prefix,
  shard=paste0("filtered.", seq_along(debstat$group_prefix), ".tirp.gz")
), "yanai_shard_libraries.csv", row.names=FALSE)

################################################################################
############################ BWA-MEM2 alignment ################################
################################################################################
reference <- normalizePath(Sys.getenv(
  "YANAI_REFERENCE", "/husky/henriksson/atrandi/bwamem_ref/ling/all.fa"
), mustWork=TRUE)
minimumScore <- as.integer(Sys.getenv("YANAI_BWAMEM2_MIN_SCORE", "20"))
stopifnot(!is.na(minimumScore), minimumScore>=0L)
alignment_job <- BascetAlignToReference(
  bascetRoot, useReference=reference, aligner="BWAMEM2",
  inputName="filtered", outputName="aligned",
  bwamem2MinScore=minimumScore,
  numThreads=numThreads,
  runner=bascetRunner.default,
  bascetInstance=bascetInstance.default
)
WaitForJob(alignment_job)

################################################################################
############################ Reads per cell and gene ###########################
################################################################################
annotation <- normalizePath(Sys.getenv(
  "YANAI_GFF", sub("\\.(fa|fasta|fna)$", ".gff", reference)
), mustWork=TRUE)
RunYanaiGeneCounts(
  bascetRoot, debstat, annotation,
  runner=bascetRunner.default, bascetInstance=bascetInstance.default
)
message("Yanai PETRI-seq preparation, alignment and gene read counting complete: ", bascetRoot)
