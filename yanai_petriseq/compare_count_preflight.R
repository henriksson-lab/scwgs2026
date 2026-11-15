#!/usr/bin/env Rscript
# Compare the production matrix with the complete-library independent preflight.
CompareYanaiCountPreflight <- function(root,library="D3_lib1") {
  stopifnot(requireNamespace("Matrix",quietly=TRUE))
  mapping <- read.csv(file.path(root,"yanai_shard_libraries.csv"))
  group <- mapping$group[mapping$library_prefix==library]
  stopifnot(length(group)==1L)
  preflight <- file.path(root,paste0("validation_counts_",library))
  oracle <- jsonlite::read_json(file.path(preflight,"bam_gene_counts.validation.json"))
  original <- file.path(preflight,"genecount.1.h5")
  production <- file.path(root,paste0("genecount.",group,".h5"))
  sha <- function(path) digest::digest(file=path,algo="sha256",serialize=FALSE)
  stopifnot(oracle$status=="verified",oracle$exact_all_cell_gene_counts,
            oracle$counts_sha256==sha(original))
  readMatrix <- function(path) {
    cells <- as.character(rhdf5::h5read(path,"/obs/_index"))
    genes <- as.character(rhdf5::h5read(path,"/var/_index"))
    stopifnot(!anyDuplicated(cells),!anyDuplicated(genes))
    mat <- methods::new("dgCMatrix",Dim=as.integer(c(length(genes),length(cells))),
      Dimnames=list(genes,cells),p=as.integer(rhdf5::h5read(path,"/X/indptr")),
      i=as.integer(rhdf5::h5read(path,"/X/indices")),
      x=as.numeric(rhdf5::h5read(path,"/X/data")))
    methods::validObject(mat)
    mat
  }
  expected <- readMatrix(original)
  actual <- readMatrix(production)
  stopifnot(setequal(rownames(expected),rownames(actual)),
            setequal(colnames(expected),colnames(actual)))
  actual <- actual[rownames(expected),colnames(expected),drop=FALSE]
  stopifnot(Matrix::nnzero(actual-expected)==0L)
  report <- list(status="verified",scope="every production cell/gene entry versus independently verified preflight",
    library=library,group=group,cells=ncol(actual),genes=nrow(actual),
    gene_read_counts=sum(actual),exact_all_cell_gene_counts=TRUE,
    production_counts_sha256=sha(production),preflight_counts_sha256=oracle$counts_sha256,
    checked_at=format(Sys.time(),tz="UTC",usetz=TRUE))
  jsonlite::write_json(report,file.path(root,paste0("gene_preflight_comparison_",library,".json")),
                       auto_unbox=TRUE,pretty=TRUE)
  message("Production counts match every independently verified preflight entry: ",library)
  invisible(report)
}
if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  root <- if(length(args)) args[1] else
    Sys.getenv("YANAI_ROOT","/husky/henriksson/atrandi/v9_rnaseq_yanai")
  CompareYanaiCountPreflight(root)
}
