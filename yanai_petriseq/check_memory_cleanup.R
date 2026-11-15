#!/usr/bin/env Rscript
# Re-run the complete validated pilot with earlier release of temporary objects.
# Preserve all genes in scale.data and compare every matrix and both embeddings.
# The active production script and its outputs are never modified here.
args <- grep("^--file=",commandArgs(),value=TRUE)
scriptDir <- dirname(dirname(normalizePath(sub("^--file=","",args[1]))))
root <- Sys.getenv("YANAI_ROOT","/husky/henriksson/atrandi/v9_rnaseq_yanai")
pilot <- file.path(root,"validation_full_D1_lib1")
trial <- file.path(root,"validation_memory_cleanup_D1_lib1")
if(dir.exists(trial)) stop("Trial directory already exists; preserve or move it before rerunning: ",trial)
dir.create(trial,recursive=TRUE)
for(name in c("genecount.1.h5","fastq_manifest.tsv","yanai_gene_count_libraries.csv")) {
  stopifnot(file.symlink(normalizePath(file.path(pilot,name),mustWork=TRUE),file.path(trial,name)))
}
lines <- readLines(file.path(scriptDir,"post_yanai.R"))
replaceLine <- function(old,new) {
  starts <- which(lines==new[1])
  if(any(vapply(starts,function(i) identical(lines[i+seq_along(new)-1L],new),logical(1))))
    return(invisible(NULL))
  hit <- which(lines==old)
  stopifnot(length(hit)==1L)
  lines <<- append(lines[-hit],new,after=hit-1L)
}
replaceLine("rm(counts);gc()",c(
  "# The full-cell selection has been written; release its original cell names.",
  "rm(counts,retained,libraryIds,index);gc()"))
replaceLine("adata <- ScaleData(adata,features=rownames(adata))",c(
  "adata <- ScaleData(adata,features=rownames(adata))",
  "gc() # Release temporary scaling blocks before PCA allocates dense matrices."))
replaceLine("adata <- RunPCA(adata,features=VariableFeatures(adata),npcs=npcs,seed.use=seed)",c(
  "adata <- RunPCA(adata,features=VariableFeatures(adata),npcs=npcs,seed.use=seed)",
  "gc() # Release temporary PCA matrices before nearest-neighbor construction."))
candidate <- file.path(trial,"post_yanai_memory_cleanup.R")
writeLines(lines,candidate)
Sys.setenv(YANAI_ROOT=trial,YANAI_PLOT_DIR=file.path(trial,"plots"))
status <- system2(file.path(R.home("bin"),"Rscript"),shQuote(candidate))
if(status!=0L) stop("Pilot candidate failed with status ",status)
suppressPackageStartupMessages(library(Seurat))
expected <- readRDS(file.path(pilot,"cache_adata_gene.RDS"))
actual <- readRDS(file.path(trial,"cache_adata_gene.RDS"))
stopifnot(identical(rownames(expected),rownames(actual)),
          identical(colnames(expected),colnames(actual)),
          identical(expected@meta.data,actual@meta.data),
          identical(VariableFeatures(expected),VariableFeatures(actual)))
for(layer in c("counts","data","scale.data")) {
  stopifnot(identical(GetAssayData(expected,layer=layer),GetAssayData(actual,layer=layer)))
}
differences <- list()
for(reduction in c("pca","gene_umap")) {
  baseline <- Embeddings(expected,reduction)
  result <- Embeddings(actual,reduction)
  stopifnot(identical(dimnames(baseline),dimnames(result)),
            isTRUE(all.equal(unname(baseline),unname(result),tolerance=1e-10)))
  differences[[reduction]] <- max(abs(baseline-result))
}
stopifnot(identical(readLines(file.path(pilot,"umap","gene_cell_selection.csv")),
                    readLines(file.path(trial,"umap","gene_cell_selection.csv"))))
report <- list(status="verified",scope="complete validated D1_lib1 pilot",
  cells=ncol(actual),genes=nrow(actual),scaled_genes=nrow(GetAssayData(actual,layer="scale.data")),
  exact_all_counts=TRUE,exact_all_normalized_expression=TRUE,
  exact_all_scaled_expression=TRUE,exact_metadata=TRUE,exact_full_cell_selection=TRUE,
  maximum_embedding_difference=differences,production_workflow_modified=FALSE,
  source_sha256=digest::digest(file=file.path(scriptDir,"post_yanai.R"),algo="sha256"),
  candidate_sha256=digest::digest(file=candidate,algo="sha256"),
  checked_at=format(Sys.time(),tz="UTC",usetz=TRUE))
jsonlite::write_json(report,file.path(trial,"memory_cleanup.validation.json"),
                     auto_unbox=TRUE,pretty=TRUE)
print(report)
