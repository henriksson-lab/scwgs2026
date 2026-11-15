#!/usr/bin/env Rscript
# Reconstruct only the normal variable-gene selection to diagnose UMAP warnings.
# All production inputs/outputs remain untouched. The actual cache must later
# confirm this feature list before using it to explain the production embedding.
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(data.table)})
if(requireNamespace("future",quietly=TRUE)) future::plan("sequential")
root <- Sys.getenv("YANAI_ROOT","/husky/henriksson/atrandi/v9_rnaseq_yanai")
scriptArg <- grep("^--file=",commandArgs(),value=TRUE)
source(file.path(dirname(normalizePath(sub("^--file=","",scriptArg[1]))),"gene_profile_fingerprint.R"))
out <- file.path(root,"umap")
set.seed(42L)
selection <- fread(file.path(out,"gene_cell_selection.csv"),
  select=c("cell_id","gene_read_counts","genes_detected","retained"))
selection <- selection[retained==TRUE]
stopifnot(nrow(selection)==864570L,!anyDuplicated(selection$cell_id))
gc()
mapping <- read.csv(file.path(root,"yanai_gene_count_libraries.csv"))
manifest <- read.delim(file.path(root,"fastq_manifest.tsv"))
stopifnot(nrow(mapping)==88L,nrow(manifest)==88L)
parts <- vector("list",nrow(mapping))
genes <- NULL
for(i in seq_len(nrow(mapping))) {
  path <- file.path(root,mapping$counts[i])
  thisGenes <- as.character(rhdf5::h5read(path,"/var/_index"))
  cells <- as.character(rhdf5::h5read(path,"/obs/_index"))
  if(is.null(genes)) genes <- thisGenes
  stopifnot(identical(genes,thisGenes))
  counts <- methods::new("dgCMatrix",Dim=as.integer(c(length(genes),length(cells))),
    Dimnames=list(genes,cells),p=as.integer(rhdf5::h5read(path,"/X/indptr")),
    i=as.integer(rhdf5::h5read(path,"/X/indices")),
    x=as.numeric(rhdf5::h5read(path,"/X/data")))
  methods::validObject(counts)
  parts[[i]] <- counts[,cells %in% selection$cell_id,drop=FALSE]
  message("Selected counts loaded: ",mapping$library_prefix[i]," (",ncol(parts[[i]])," cells)")
}
rm(counts,cells,thisGenes);gc()
# Pairwise concatenation avoids repeatedly copying the entire growing matrix.
while(length(parts)>1L) {
  parts <- lapply(split(seq_along(parts),ceiling(seq_along(parts)/2)),function(index) {
    if(length(index)==1L) parts[[index]] else Matrix::cbind2(parts[[index[1]]],parts[[index[2]]])
  })
}
counts <- parts[[1]]
rm(parts);gc()
stopifnot(setequal(colnames(counts),selection$cell_id))
# Zorn's MergeBascetCountMatrix sorts raw feature IDs before Seurat replaces
# underscores. Match that ordering as well as the exact cell ordering: ties in
# VST/loess and variable-feature ranks otherwise need not resolve identically.
counts <- counts[as.character(sort(rownames(counts))),
                 match(selection$cell_id,colnames(counts)),drop=FALSE]
stopifnot(identical(colnames(counts),selection$cell_id),
          all(Matrix::colSums(counts)==selection$gene_read_counts),
          all(Matrix::colSums(counts>0)==selection$genes_detected))
object <- CreateSeuratObject(counts=counts,project="Yanai",min.cells=0,min.features=0)
rm(counts);gc()
object <- object[Matrix::rowSums(GetAssayData(object,layer="counts"))>0,]
object <- NormalizeData(object)
object <- FindVariableFeatures(object,selection.method="vst",nfeatures=min(2000L,nrow(object)))
features <- VariableFeatures(object)
write.csv(object[["RNA"]][[]],file.path(out,"gene_feature_metadata.reconstructed.csv"))
counts <- GetAssayData(object,layer="counts")[features,,drop=FALSE]
coverage <- data.table(cell_id=colnames(counts),variable_genes=as.numeric(Matrix::colSums(counts>0)),
                       variable_gene_reads=as.numeric(Matrix::colSums(counts)))
coverage[,library:=sub("^(D[0-9]+_lib[0-9]+)_.*$","\\1",cell_id)]
coverage[,species:=manifest$scientific_name[match(library,manifest$library_prefix)]]
stopifnot(!anyNA(coverage$species),uniqueN(coverage$library)==88L)
summarize <- function(x) list(cells=nrow(x),cells_with_zero_variable_genes=sum(x$variable_genes==0),
  cells_with_at_most_one_variable_gene=sum(x$variable_genes<=1),
  detected_variable_gene_quantiles=as.list(quantile(x$variable_genes,c(0,.01,.1,.25,.5,.75,.9,1))))
featureGenome <- ifelse(grepl("^b[0-9]+$",features),"Escherichia coli",
  ifelse(grepl("^SAUSA300-RS",features),"Staphylococcus aureus",
    ifelse(grepl("^BSU-",features),"Bacillus subtilis","unclassified_gene_id")))
fwrite(coverage,file.path(out,"gene_variable_gene_coverage.tsv"),sep="\t")
fwrite(data.table(feature=features,reference_species=featureGenome),
       file.path(out,"gene_variable_features.tsv"),sep="\t")
writeLines(features,file.path(out,"gene_variable_features.reconstructed.txt"))
report <- list(status="profiled_pending_actual_cache_comparison",cells=ncol(object),
  expressed_genes=nrow(object),variable_features=length(features),
  exact_selected_input_count_totals=TRUE,exact_selected_input_detected_genes=TRUE,
  all=summarize(coverage),species=lapply(split(coverage,coverage$species),summarize),
  feature_genomes=as.list(table(featureGenome)),actual_umap_features_checked=FALSE,
  selection_modified=FALSE,Seurat_version=as.character(packageVersion("Seurat")),
  gene_order_policy="Zorn sort(raw feature IDs) before Seurat name normalization",
  gene_sort_locale=Sys.getlocale("LC_COLLATE"),
  input_fingerprint=YanaiGeneProfileFingerprint(root),
  checked_at=format(Sys.time(),tz="UTC",usetz=TRUE))
jsonlite::write_json(report,file.path(out,"gene_variable_gene_coverage.json"),
                     auto_unbox=TRUE,pretty=TRUE)
print(report)
