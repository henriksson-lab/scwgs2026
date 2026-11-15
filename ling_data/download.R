#!/usr/bin/env Rscript
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
here <- dirname(normalizePath(script,mustWork=TRUE))
source(file.path(dirname(here),"dataset_download_common.R"))
options <- ParseDownloadOptions(list(root="/husky/henriksson/atrandi/v9_ling",
  manifest=file.path(here,"manifest.tsv"),workers="6",library=character()),repeated="library")
dir.create(options$root,recursive=TRUE,showWarnings=FALSE)
root <- normalizePath(options$root,mustWork=TRUE)
manifest <- read.delim(options$manifest,colClasses="character")
rows <- if(length(options$library)) manifest[manifest$library %in% options$library,,drop=FALSE] else manifest
stopifnot(nrow(rows)>0L,all(options$library %in% manifest$library))
CopyDownloadMetadata(options$manifest,file.path(root,"manifest.tsv"))
# Preserve the earlier human-readable snapshot when repository docs evolve.
# Source manifests and deposited metadata remain strict checksum comparisons.
if(!file.exists(file.path(root,"README.md"))) {
  CopyDownloadMetadata(file.path(here,"README.md"),file.path(root,"README.md"))
}
for(path in list.files(here,pattern="^zenodo_.*[.]json$",full.names=TRUE)) {
  CopyDownloadMetadata(path,file.path(root,basename(path)))
}
message("Preparing ",nrow(rows)," source files (",sum(as.numeric(rows$bytes))/1e9," GB)")
RunDownloadTasks(seq_len(nrow(rows)),options$workers,function(i) {
  row <- rows[i,]
  folder <- file.path(root,"source",row$library)
  dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  path <- file.path(folder,row$filename)
  FetchDownload(path,row$url,row$bytes,row$md5)
  receipt <- paste0(path,".verified.json")
  if(file.exists(receipt)) {
    old <- jsonlite::read_json(receipt)
    stopifnot(old$md5==row$md5,old$bytes==as.numeric(row$bytes))
  } else PublishDownloadJSON(list(zenodo_record=row$record,source_url=row$url,
    bytes=as.numeric(row$bytes),md5=row$md5,verified_at=as.numeric(Sys.time()),
    library=row$library,barcode_a=row$barcode_a,role=row$role),receipt)
  message("Verified ",row$library," / ",row$filename)
})
message("All ",nrow(rows)," selected source FASTQs verified in ",file.path(root,"source"))
