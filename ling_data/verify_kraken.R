#!/usr/bin/env Rscript
# Verify sparse structure and per-barcode conservation of selected read pairs.
VerifyLingKrakenInput <- function(root,sharding=NULL) {
  if(is.null(sharding)) sharding <- readRDS(file.path(root,"ling_sharding.rds"))
  stopifnot(length(sharding$group_prefix)==4L)
  report <- lapply(seq_along(sharding$group_prefix),function(i) {
    path <- file.path(root,paste0("kraken_counts.",i,".h5"))
    ids <- as.character(rhdf5::h5read(path,"/obs/_index"))
    shape <- rhdf5::h5readAttributes(path,"/X")$shape
    ptr <- as.integer(rhdf5::h5read(path,"/X/indptr"))
    values <- as.numeric(rhdf5::h5read(path,"/X/data"))
    indices <- as.integer(rhdf5::h5read(path,"/X/indices"))
    features <- as.character(rhdf5::h5read(path,"/var/_index"))
    unmapped <- as.numeric(rhdf5::h5read(path,"/obs/_unmapped"))
    stopifnot(setequal(ids,sharding$list_picked_cells[[i]]),!anyDuplicated(ids),
              shape[1]==length(ids),shape[2]==length(features),!anyDuplicated(features),
              length(ptr)==length(ids)+1L,ptr[1]==0,tail(ptr,1)==length(values),
              length(indices)==length(values),all(diff(ptr)>=0),all(values>=0),
              all(indices>=0 & indices<shape[2]),length(unmapped)==length(ids))
    cumulative <- c(0,cumsum(values))
    counts <- cumulative[ptr[-1]+1L]-cumulative[ptr[-length(ptr)]+1L]
    bc <- sharding$list_group_bc_count[[i]]
    expected <- bc$count[match(ids,bc$cellid)]
    stopifnot(!anyNA(expected),all(counts+unmapped==expected))
    data.frame(library=sharding$group_prefix[i],barcodes=length(ids),features=shape[2],
               input_read_pairs=sum(expected),classified_pairs=sum(counts),
               unclassified_pairs=sum(unmapped),
               exact_barcode_coverage=TRUE,exact_read_conservation=TRUE)
  })
  report <- do.call(rbind,report)
  write.table(report,file.path(root,"kraken_input_validation.tsv"),
              sep="\t",quote=FALSE,row.names=FALSE)
  report
}
if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  root <- if(length(args)) args[1] else
    Sys.getenv("LING_ROOT","/husky/henriksson/atrandi/v9_ling")
  print(VerifyLingKrakenInput(root),row.names=FALSE)
}
