#!/usr/bin/env Rscript
# Read-only diagnostic on the COMPLETE validated D1_lib1 pilot.
# Does not change the production run or any saved Seurat object.
suppressPackageStartupMessages(library(Seurat))
if(requireNamespace("future",quietly=TRUE)) future::plan("sequential")
root <- Sys.getenv("YANAI_ROOT","/husky/henriksson/atrandi/v9_rnaseq_yanai")
pilot <- file.path(root,"validation_full_D1_lib1")
object <- readRDS(file.path(pilot,"cache_adata_gene.RDS"))
features <- VariableFeatures(object)
expected <- GetAssayData(object,layer="scale.data")[features,,drop=FALSE]
expectedPca <- Embeddings(object,"pca")
candidate <- ScaleData(object,features=features,verbose=FALSE)
actual <- GetAssayData(candidate,layer="scale.data")[features,,drop=FALSE]
stopifnot(identical(dimnames(expected),dimnames(actual)),
          identical(as.vector(expected),as.vector(actual)),
          nrow(candidate)==nrow(object),ncol(candidate)==ncol(object))
candidate <- RunPCA(candidate,features=features,npcs=ncol(expectedPca),
                    seed.use=42L,verbose=FALSE)
actualPca <- Embeddings(candidate,"pca")
stopifnot(identical(dimnames(expectedPca),dimnames(actualPca)),
  isTRUE(all.equal(unname(expectedPca),unname(actualPca),
                   tolerance=1e-10,check.attributes=FALSE)))
report <- list(status="verified",scope="complete validated D1_lib1 pilot",
  cells=ncol(object),all_genes=nrow(object),pca_genes=length(features),
  exact_scaled_pca_inputs=TRUE,pca_embeddings_equal=TRUE,
  maximum_pca_difference=max(abs(expectedPca-actualPca)),
  production_workflow_modified=FALSE,
  checked_at=format(Sys.time(),tz="UTC",usetz=TRUE))
jsonlite::write_json(report,file.path(pilot,"scaling_equivalence.validation.json"),
                     auto_unbox=TRUE,pretty=TRUE)
print(report)
