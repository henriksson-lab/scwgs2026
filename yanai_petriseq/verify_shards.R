#!/usr/bin/env Rscript
# Check all library labels, sampling only the first record of each TIRP shard.
# Full records are read by alignment; this is not a whole-file sequence audit.
VerifyYanaiShards <- function(root) {
  mapping <- read.csv(file.path(root,"yanai_shard_libraries.csv"))
  manifest <- read.delim(file.path(root,"fastq_manifest.tsv"))
  stopifnot(nrow(mapping)==88L,!anyDuplicated(mapping$library_prefix),
            setequal(mapping$library_prefix,manifest$library_prefix))
  report <- lapply(seq_len(nrow(mapping)),function(i) {
    con <- gzfile(file.path(root,mapping$shard[i]),"rt")
    line <- tryCatch(readLines(con,n=1L,warn=FALSE),finally=close(con))
    stopifnot(length(line)==1L)
    fields <- strsplit(line,"\t",fixed=TRUE)[[1]]
    stopifnot(length(fields)==8L,
              startsWith(fields[1],paste0(mapping$library_prefix[i],"_")),
              fields[4]=="",fields[6]=="",nchar(fields[5])==26L,
              nchar(fields[7])==26L,nchar(fields[8])==7L)
    data.frame(library=mapping$library_prefix[i],shard=mapping$shard[i],
               first_cell=fields[1],library_identity_valid=TRUE,
               first_record_layout_valid=TRUE)
  })
  report <- do.call(rbind,report)
  write.table(report,file.path(root,"yanai_shard_validation.tsv"),
              sep="\t",quote=FALSE,row.names=FALSE)
  report
}
if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  root <- if(length(args)) args[1] else
    Sys.getenv("YANAI_ROOT","/husky/henriksson/atrandi/v9_rnaseq_yanai")
  invisible(VerifyYanaiShards(root))
  message("All 88 library labels and sampled PETRI-seq record layouts verified")
}
