#!/usr/bin/env Rscript
# Compare the complete merged matrix with every original library feature/count.
VerifyLingMergedKraken <- function(root) {
  zorn <- path.expand(Sys.getenv("LING_ZORN_DIR","~/github/zorn"))
  sourceFile <- file.path(zorn,"R","countmatrix.R")
  if(!exists("ReadBascetCountMatrix_one",mode="function")) source(sourceFile)
  inputs <- lapply(1:4,function(i) {
    ReadBascetCountMatrix_one(file.path(root,paste0("kraken_counts.",i,".h5")))
  })
  elapsed <- system.time(merged <- MergeBascetCountMatrix(inputs))
  stopifnot(!anyDuplicated(rownames(merged@X)),all(merged@X@x>=0),
            nrow(merged@X)==sum(vapply(inputs,function(x)nrow(x@X),integer(1))))
  for(input in inputs) {
    rows <- rownames(input@X);cols <- colnames(input@X)
    stopifnot(Matrix::nnzero(merged@X[rows,cols,drop=FALSE]-input@X)==0L)
    extra <- setdiff(colnames(merged@X),cols)
    stopifnot(Matrix::nnzero(merged@X[rows,extra,drop=FALSE])==0L)
    index <- match(rows,rownames(merged@X))
    stopifnot(identical(as.character(merged@obs$`_index`[index]),as.character(input@obs$`_index`)),
              identical(as.numeric(merged@obs$`_unmapped`[index]),as.numeric(input@obs$`_unmapped`)))
  }
  report <- list(status="verified",barcodes=nrow(merged@X),features=ncol(merged@X),
    classified_pairs=sum(merged@X),exact_per_feature_counts=TRUE,
    zero_extra_feature_counts=TRUE,exact_metadata=TRUE,
    merge_elapsed=unname(elapsed["elapsed"]),source_md5=unname(tools::md5sum(sourceFile)))
  jsonlite::write_json(report,file.path(root,"kraken_merged_validation.json"),auto_unbox=TRUE,pretty=TRUE)
  report
}
if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  root <- if(length(args)) args[1] else
    Sys.getenv("LING_ROOT","/husky/henriksson/atrandi/v9_ling")
  print(VerifyLingMergedKraken(root))
}
