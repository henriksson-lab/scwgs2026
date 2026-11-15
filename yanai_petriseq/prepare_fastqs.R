#!/usr/bin/env Rscript
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
here <- dirname(normalizePath(script,mustWork=TRUE))
source(file.path(dirname(here),"dataset_download_common.R"))
options <- ParseDownloadOptions(list(root="/husky/henriksson/atrandi/v9_rnaseq_yanai",
  manifest=file.path(here,"fastq_manifest.tsv"),`raw-dir`="",`output-dir`="",
  download=FALSE,workers="4",threads="2",run=character()),
  flags="download",repeated="run")
threads <- as.integer(options$threads)
stopifnot(length(threads)==1L,!is.na(threads),threads>0L)
dir.create(options$root,recursive=TRUE,showWarnings=FALSE)
root <- normalizePath(options$root,mustWork=TRUE)
rawDir <- if(nzchar(options[["raw-dir"]])) options[["raw-dir"]] else file.path(root,"raw_gzip")
outputDir <- if(nzchar(options[["output-dir"]])) options[["output-dir"]] else file.path(root,"raw_bgzf")
for(path in c(rawDir,outputDir)) dir.create(path,recursive=TRUE,showWarnings=FALSE)
rawDir <- normalizePath(rawDir,mustWork=TRUE)
outputDir <- normalizePath(outputDir,mustWork=TRUE)
manifest <- read.delim(options$manifest,colClasses="character")
rows <- if(length(options$run)) manifest[manifest$run_accession %in% options$run,,drop=FALSE] else manifest
stopifnot(nrow(rows)>0L,!anyDuplicated(rows$library_prefix),
          all(options$run %in% manifest$run_accession))
CopyDownloadMetadata(options$manifest,file.path(root,"fastq_manifest.tsv"))
rawPointer <- file.path(root,"rawdata.txt")
if(file.exists(rawPointer)) stopifnot(identical(readLines(rawPointer),outputDir)) else
  writeLines(outputDir,rawPointer)
writeLines(as.character(Sys.getpid()),file.path(root,"download.pid"))
tasks <- expand.grid(row=seq_len(nrow(rows)),mate=1:2)
RunDownloadTasks(seq_len(nrow(tasks)),options$workers,function(task) {
  row <- rows[tasks$row[task],]
  mate <- tasks$mate[task]
  run <- row$run_accession;library <- row$library_prefix
  stopifnot(grepl("^SRR[0-9]+$",run),grepl("^D[0-9]+_lib[0-9]+$",library))
  urls <- strsplit(row$fastq_ftp,";",fixed=TRUE)[[1]]
  md5 <- strsplit(row$fastq_md5,";",fixed=TRUE)[[1]]
  bytes <- strsplit(row$fastq_bytes,";",fixed=TRUE)[[1]]
  stopifnot(length(urls)==2L,length(md5)==2L,length(bytes)==2L,
            endsWith(urls[mate],paste0("/",run,"_",mate,".fastq.gz")))
  input <- file.path(rawDir,paste0(run,"_",mate,".fastq.gz"))
  output <- file.path(outputDir,paste0(library,"_S1_L001_R",mate,"_001.fastq.gz"))
  receipt <- paste0(output,".json")
  identity <- list(run=run,library=library,mate=mate,source_md5=md5[mate])
  FetchDownload(input,paste0("https://",urls[mate]),bytes[mate],md5[mate],options$download)
  if(file.exists(output) && file.exists(receipt)) {
    old <- jsonlite::read_json(receipt)
    stopifnot(all(vapply(names(identity),function(k) {
      identical(as.character(old[[k]]),as.character(identity[[k]]))
    },logical(1))),file.info(output)$size==old$output_bytes)
    DownloadCommand("bgzip",c("-t",output))
    message("Verified existing ",basename(output))
    return(invisible())
  }
  if(file.exists(output)||file.exists(receipt)) stop("Incomplete existing output: ",output)
  temp <- tempfile(paste0(".",basename(output),"."),tmpdir=outputDir,fileext=".part")
  on.exit(unlink(temp))
  # pipefail includes gzip's CRC result; bgzip preserves the exact inflated bytes.
  command <- paste("gzip -dc --",shQuote(input),"| bgzip -@",threads,"-c >",shQuote(temp))
  DownloadCommand("bash",c("-o","pipefail","-c",command))
  DownloadCommand("bgzip",c("-t",temp))
  if(!file.link(temp,output)) stop("Cannot publish ",output)
  identity$output_bytes <- file.info(output)$size
  identity$completed_at <- as.numeric(Sys.time())
  PublishDownloadJSON(identity,receipt)
  message("Prepared ",basename(output))
})
message("All ",nrow(rows)," selected libraries ready; rawdata.txt -> ",outputDir)
