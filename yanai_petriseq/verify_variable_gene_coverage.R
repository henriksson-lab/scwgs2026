# Full-cache verification of the optional native-count warning investigation.
VerifyYanaiVariableGeneCoverage <- function(root,object,scriptDir) {
  out <- file.path(root,"umap")
  checked <- list()
  # Check the warning investigation against the actual full saved object.
  featureFile <- file.path(out,"gene_variable_features.reconstructed.txt")
  profileCurrent <- FALSE
  if(file.exists(featureFile)) {
    profileFile <- file.path(out,"gene_variable_gene_coverage.json")
    profile <- jsonlite::read_json(profileFile,simplifyVector=TRUE)
    source(file.path(scriptDir,"yanai_petriseq","gene_profile_fingerprint.R"))
    profileCurrent <- identical(profile$input_fingerprint,YanaiGeneProfileFingerprint(root))
    if(!profileCurrent) {
      message("Optional Yanai feature-coverage diagnostic is stale for these inputs; regenerate it with profile_variable_gene_coverage.R")
      checked$gene_feature_coverage <- list(status="stale_optional_diagnostic")
    }
  }
  if(profileCurrent) {
    features <- readLines(featureFile)
    stopifnot(identical(unname(features),unname(VariableFeatures(object))))
    coverage <- data.table::fread(file.path(out,"gene_variable_gene_coverage.tsv"))
    stopifnot(nrow(coverage)==ncol(object),!anyDuplicated(coverage$cell_id),
              setequal(coverage$cell_id,colnames(object)))
    coverage <- coverage[match(colnames(object),coverage$cell_id)]
    values <- GetAssayData(object,layer="counts")[features,,drop=FALSE]
    stopifnot(all(Matrix::colSums(values>0)==coverage$variable_genes),
              all(Matrix::colSums(values)==coverage$variable_gene_reads))
    # Inspect stored dimensions only; LayerData may subset/copy the dense layer
    # even with fast=TRUE. The all-gene matrix is about 81 GB for this dataset.
    assay <- object[["RNA"]]
    scaled <- if("layers" %in% methods::slotNames(assay))
      methods::slot(assay,"layers")[["scale.data"]] else methods::slot(assay,"scale.data")
    stopifnot(nrow(object)==profile$expressed_genes,nrow(scaled)==nrow(object))
    rm(scaled,assay)
    profile$status <- "verified_against_actual_full_cache"
    profile$actual_umap_features_checked <- TRUE
    profile$actual_umap_variable_counts_checked <- TRUE
    jsonlite::write_json(profile,profileFile,auto_unbox=TRUE,pretty=TRUE)
    checked$gene_feature_coverage <- list(features=length(features),cells=nrow(coverage),
      exact_actual_feature_list=TRUE,exact_actual_variable_gene_counts=TRUE,
      feature_genomes=profile$feature_genomes)
    annotationFile <- file.path(out,"gene_variable_features.annotations.json")
    annotationTable <- data.table::fread(file.path(out,"gene_variable_features.annotated.tsv"))
    stopifnot(identical(annotationTable$feature,features))
    annotation <- jsonlite::read_json(annotationFile,simplifyVector=TRUE)
    annotation$status <- "verified_against_actual_full_cache"
    annotation$actual_umap_features_checked <- TRUE
    jsonlite::write_json(annotation,annotationFile,auto_unbox=TRUE,pretty=TRUE)
    checked$gene_feature_coverage$gene_biotypes <- annotation$genome_biotype_counts
    log <- readLines(file.path(root,"all_workflow.log"),warn=FALSE)
    failure <- grep("smooth knn distance failures",log,value=TRUE)
    recall <- grep("Annoy recall =",log,value=TRUE)
    fallback <- grep("connected components, falling back to",log,value=TRUE)
    diagnostics <- list(cells=ncol(object),expressed_genes=nrow(object),
      input_fingerprint=profile$input_fingerprint,
      all_expressed_genes_scaled=TRUE,variable_genes=length(features),
      reported_smooth_knn_distance_failures=if(length(failure))
        as.integer(sub(".* ([0-9]+) smooth knn distance failures.*","\\1",tail(failure,1))) else NULL,
      reported_annoy_recall_percent=if(length(recall))
        as.numeric(sub(".*Annoy recall = ([0-9.]+)%.*","\\1",tail(recall,1))) else NULL,
      reported_initialization_fallback=if(length(fallback)) tail(fallback,1) else NULL,
      feature_genomes=profile$feature_genomes,
      note="Global variable-gene selection represents genomes unevenly; warnings are recorded without changing cell selection or the normal workflow.")
    jsonlite::write_json(diagnostics,file.path(out,"gene_umap_diagnostics.json"),
                         auto_unbox=TRUE,pretty=TRUE)
    checked$gene_umap_diagnostics <- diagnostics
    rm(annotationTable,annotation,log,diagnostics);gc()
    rm(coverage,values,profile);gc()
  }
  checked
}
