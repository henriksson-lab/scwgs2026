#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
root <- normalizePath(args[1],mustWork=TRUE)
binary <- normalizePath(args[2],mustWork=TRUE)
sha <- function(path) digest::digest(file=path,algo="sha256",serialize=FALSE)
suffixes <- c(".0123",".amb",".ann",".bwt.2bit.64",".pac")
artifacts <- setNames(lapply(suffixes,function(suffix) {
  path <- file.path(root,paste0("all.fa",suffix))
  stopifnot(file.exists(path),file.info(path)$size>0)
  list(bytes=file.info(path)$size,sha256=sha(path))
}),suffixes)
report <- list(reference=file.path(root,"all.fa"),reference_sha256=sha(file.path(root,"all.fa")),
  bascet_binary=binary,bascet_version=system2(binary,"--version",stdout=TRUE),
  command=c("exttool","bwa-mem2","index",file.path(root,"all.fa")),artifacts=artifacts)
path <- file.path(root,"index.json")
if(file.exists(path)) {
  old <- jsonlite::read_json(path)
  stopifnot(old$reference_sha256==report$reference_sha256)
  for(suffix in suffixes) stopifnot(old$artifacts[[suffix]]$bytes==artifacts[[suffix]]$bytes,
                                   old$artifacts[[suffix]]$sha256==artifacts[[suffix]]$sha256)
  message("Existing reference and all five index artifact checksums verified")
} else jsonlite::write_json(report,path,auto_unbox=TRUE,pretty=TRUE)
