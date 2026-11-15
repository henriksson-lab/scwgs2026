#!/usr/bin/env Rscript
# Verify one COMPLETE Yanai library using a locally built Bascet checkout.
# Run: Rscript /husky/henriksson/atrandi/revision_scripts/prep_yanai_dev.R
# Overrides: YANAI_LIBRARY=D1_lib1, YANAI_ROOT, YANAI_OUTPUT_DIR,
# YANAI_BASCET_DIR=~/github/bascet, YANAI_BASCET_TARGET=release|debug,
# YANAI_ZORN_DIR=~/github/zorn, YANAI_THREADS=10, YANAI_MEM=16g.
# YANAI_BUILD=0 skips the local cargo build. No remote version is fetched.
# Outputs must be new, so a changed binary cannot silently reuse an old test.
# This script debarcodes and verifies only; prep_yanai.R provides alignment.

scriptFile <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
if(is.null(scriptFile)) {
  arg <- grep("^--file=", commandArgs(), value=TRUE)
  if(!length(arg)) stop("Run with Rscript or source()")
  scriptFile <- sub("^--file=", "", arg[1])
}
scriptDir <- dirname(normalizePath(scriptFile, mustWork=TRUE))
required <- c("stringr", "fs", "data.table", "readr", "processx")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop("Missing R dependencies: ", paste(missing, collapse=", "))

root <- normalizePath(Sys.getenv("YANAI_ROOT", "/husky/henriksson/atrandi/v9_rnaseq_yanai"), mustWork=TRUE)
libraryName <- Sys.getenv("YANAI_LIBRARY", "D1_lib1")
manifest <- read.delim(file.path(root, "fastq_manifest.tsv"), stringsAsFactors=FALSE)
selected <- manifest[manifest$library_prefix==libraryName,,drop=FALSE]
if(nrow(selected)!=1L) stop("YANAI_LIBRARY must identify exactly one manifest library")
rawDir <- readLines(file.path(root, "rawdata.txt"), warn=FALSE)
rawDir <- trimws(rawDir[nzchar(trimws(rawDir))])
if(length(rawDir)!=1L) stop("rawdata.txt must contain one prepared BGZF directory")
rawDir <- normalizePath(rawDir, mustWork=TRUE)
r1 <- paste0(libraryName, "_S1_L001_R1_001.fastq.gz")
r2 <- paste0(libraryName, "_S1_L001_R2_001.fastq.gz")
stopifnot(all(file.exists(file.path(rawDir,c(r1,r2)))))
rawmeta <- data.frame(prefix=libraryName, r1=r1, r2=r2, dir=rawDir)

output <- path.expand(Sys.getenv("YANAI_OUTPUT_DIR", file.path(root,"dev_debarcode",libraryName)))
if(dir.exists(output) && length(list.files(output, all.files=TRUE, no..=TRUE))) {
  stop("Output directory is not empty; choose a fresh YANAI_OUTPUT_DIR: ", output)
}
dir.create(output, recursive=TRUE, showWarnings=FALSE)
output <- normalizePath(output, mustWork=TRUE)
zornDir <- path.expand(Sys.getenv("YANAI_ZORN_DIR", "~/github/zorn"))
for(file in c("job_general.R", "job_local.R", "job_slurm.R", "jobscript.R",
              "bascet_file.R", "zorn.R", "shell.R", "container.R", "debarcoding.R")) {
  source(file.path(zornDir,"R",file))
}
if(!"petriseq" %in% eval(formals(BascetDebarcode)$chemistry)) {
  stop("The local Zorn checkout needs PETRI-seq support")
}

devDir <- normalizePath(path.expand(Sys.getenv("YANAI_BASCET_DIR", "~/github/bascet")), mustWork=TRUE)
target <- Sys.getenv("YANAI_BASCET_TARGET", "release")
if(!target %in% c("release","debug")) stop("YANAI_BASCET_TARGET must be release or debug")
if(Sys.getenv("YANAI_BUILD", "1")=="1") {
  previous <- getwd()
  status <- tryCatch({
    setwd(devDir) # Load this checkout's .cargo/config.toml.
    args <- c("build", "--locked", if(target=="release") "--release", "--bin", "bascet")
    system2("cargo", args)
  }, finally=setwd(previous))
  if(status!=0) stop("Local Bascet build failed")
}
setwd(output)
# This is Zorn's DEV MODE: directly execute target/{release,debug}/bascet.
# No download, container, GitHub release, git pull, or git push is involved.
instance <- getBascetDevDir(devDir, targetType=target, containerMem="0B")
threads <- as.integer(Sys.getenv("YANAI_THREADS", "10"))
stopifnot(!is.na(threads), threads>=6L)
runner <- LocalRunner(direct=TRUE, showScript=TRUE, ncpu=as.character(threads),
                      mem=Sys.getenv("YANAI_MEM", "16g"))
write.csv(selected, "selected_library.csv", row.names=FALSE)
writeLines(c(paste("binary",instance@bin), paste("binary_md5",tools::md5sum(instance@bin)),
             paste("library",libraryName), paste("run",selected$run_accession),
             paste("R1",file.path(rawDir,r1)), paste("R2",file.path(rawDir,r2))),
           "local_bascet.txt")
message("Debarcoding complete library ",libraryName," (",selected$run_accession,") using ",instance@bin)
job <- BascetDebarcode(output, rawmeta, chemistry="petriseq", outputName="debarcoded",
                      maxShardSize="200g", numThreads=threads, compressionLevel=1L,
                      runner=runner, bascetInstance=instance)
WaitForJob(job)
tirp <- file.path(output,"debarcoded.1.tirp.gz")
stopifnot(file.exists(tirp), file.exists(paste0(tirp,".hist")))

# Independent streaming oracle checks EVERY input pair and EVERY output record.
# The multiset checksum is insensitive to sorting but preserves multiplicity.
status <- system2("bash", shQuote(c(
  file.path(scriptDir,"run_yanai_verify_debarcoded.sh"),
  "--r1",file.path(rawDir,r1), "--r2",file.path(rawDir,r2),
  "--library",libraryName, "--tirp",tirp, "--hist",paste0(tirp,".hist"),
  "--barcodes",file.path(devDir,"crates/bascet-cli/src/barcode/petriseq_barcodes.tsv"),
  "--report",file.path(output,"debarcoding_validation.json")
)))
if(status!=0) stop("Independent full-library debarcoding validation failed")
message("Full-library local Bascet validation passed: ",output)
