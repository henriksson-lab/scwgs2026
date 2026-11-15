# Shared count/sketch preparation for prep_gorg.R / prep_ling.R.
# Outputs feed downstream UMAP analysis; no read QC or assembly tools are run.
# Source these workflows with Rscript, or source() them in an R session.

InitializeDNAWorkflow <- function(dataset, defaultRoot) {
  env <- function(suffix, default="") Sys.getenv(paste0(dataset, "_", suffix), default)
  root <- normalizePath(env("ROOT", defaultRoot), mustWork=TRUE)
  required <- c("stringr", "fs", "data.table", "readr", "ggplot2", "naturalsort", "processx")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
  if(length(missing)) stop("Missing Zorn dependencies: ", paste(missing, collapse=", "))
  zornDir <- path.expand(env("ZORN_DIR", "~/github/zorn"))
  for(file in c(
    "job_general.R", "job_local.R", "job_slurm.R", "jobscript.R",
    "bascet_file.R", "zorn.R", "shell.R", "zorn_aggr.R",
    "container.R", "debarcoding.R", "count_kmer.R", "countsketch.R",
    "kraken.R"
  )) source(file.path(zornDir, "R", file), local=.GlobalEnv)

  mode <- env("RUNNER", if(Sys.info()[["nodename"]]=="beagle") "local" else "slurm")
  threads <- as.integer(env("THREADS", "10"))
  stopifnot(!is.na(threads), threads>=6L)
  defaultDev <- if(mode=="local") "~/github/bascet" else "/home/m/mahogny/github/bascet"
  setwd(root)
  instance <- getBascetDevDir(
    path.expand(env("BASCET_DIR", defaultDev)), targetType=env("BASCET_TARGET", "release")
  )
  if(mode=="local") {
    runner <- LocalRunner(direct=TRUE, showScript=FALSE,
                          ncpu=as.character(threads), mem=env("MEM", "50g"))
  } else if(mode=="slurm") {
    runner <- SlurmRunner(
      account=env("SLURM_ACCOUNT", "hpc2n2026-186"),
      ncpu=as.character(threads), mem=env("MEM", "40g"),
      logTime=TRUE, verbose=TRUE, time="0-14:00:00"
    )
  } else stop("Runner must be local or slurm")
  list(root=root, runner=runner, instance=instance,
       threads=threads, dataset=dataset, env=env)
}

RunDNAFeatures <- function(config, inputName) {
  root <- config$root
  runner <- config$runner
  instance <- config$instance
  input <- detectShardsForFile(root, inputName)
  if(!length(input)) stop("No ", inputName, " TIRP shards found")

  ##############################################################################
  ############################ Count sketch ####################################
  ##############################################################################
  WaitForJob(BascetRunCountsketch(
    root, inputName=inputName, kmerSize=31,
    numThreads=config$threads, runner=runner, bascetInstance=instance
  ))

  ##############################################################################
  ############################ KRAKEN ##########################################
  ##############################################################################
  if(config$env("RUN_KRAKEN", "1")=="1") {
    candidates <- c("/husky/henriksson/atrandi/kraken_ref/standard-8",
                    "/data/henlab/kraken/standard-8",
                    "/home/m/mahogny/mystore/atrandi/kraken_ref/standard-8")
    defaultDB <- candidates[dir.exists(candidates)][1]
    krakenDB <- config$env("KRAKEN_DB", if(is.na(defaultDB)) "" else defaultDB)
    if(!dir.exists(krakenDB)) stop("Set ", config$dataset, "_KRAKEN_DB to a Kraken database")
    WaitForJob(BascetRunKraken(
      root, inputName=inputName, useKrakenDB=krakenDB,
      numThreads=config$threads, runner=runner, bascetInstance=instance
    ))
  }

  ##############################################################################
  ############################ Informative 31-mers #############################
  ##############################################################################
  WaitForJob(BascetComputeMinhash(
    root, inputName=inputName, kmerSize=31,
    numThreads=config$threads, runner=runner, bascetInstance=instance
  ))
  WaitForJob(BascetMakeMinhashHistogram(root, runner=runner, bascetInstance=instance))
  kmerFile <- file.path(root, "use_kmers.txt")
  if(!file.exists(kmerFile)) {
    histogram <- BascetReadMinhashHistogram(root)
    if(!nrow(histogram)) stop("No 31-mers found in the DNA input")
    set.seed(666)
    histogram$tie_breaker <- sample.int(nrow(histogram))
    histogram <- histogram[order(histogram$cnt, histogram$tie_breaker, decreasing=TRUE),]
    writeLines(as.character(head(histogram$kmer, 100000)), kmerFile)
  }
  WaitForJob(BascetQueryFq(
    root, inputName=inputName, useKMERs=readLines(kmerFile),
    numThreads=config$threads, runner=runner, bascetInstance=instance
  ))
}
