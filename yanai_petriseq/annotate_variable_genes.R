#!/usr/bin/env Rscript
# Annotate the reconstructed feature inventory; final audit checks the actual cache.
root <- Sys.getenv("YANAI_ROOT","/husky/henriksson/atrandi/v9_rnaseq_yanai")
suppressPackageStartupMessages(library(data.table))
features <- fread(file.path(root,"umap","gene_variable_features.tsv"))
mapping <- read.csv(file.path(root,"yanai_gene_count_libraries.csv"))
paths <- unique(mapping$annotation)
stopifnot(length(paths)==1L)
gff <- read.delim(paths,header=FALSE,comment.char="#",quote="",colClasses="character")
gff <- gff[gff$V3=="gene",]
attributes <- strsplit(gff$V9,";",fixed=TRUE)
attribute <- function(key) vapply(attributes,function(fields) {
  prefix <- paste0(key,"=")
  value <- fields[startsWith(fields,prefix)]
  if(length(value)==1L) utils::URLdecode(substring(value,nchar(prefix)+1L)) else NA_character_
},character(1))
locus <- attribute("locus_tag")
namesInCache <- gsub("_","-",locus,fixed=TRUE)
stopifnot(!anyNA(locus),!anyDuplicated(namesInCache))
index <- match(features$feature,namesInCache)
stopifnot(!anyNA(index))
features[,locus_tag:=locus[index]]
features[,gene_name:=attribute("Name")[index]]
features[,gene_biotype:=attribute("gene_biotype")[index]]
features[is.na(gene_biotype),gene_biotype:="unspecified"]
fwrite(features,file.path(root,"umap","gene_variable_features.annotated.tsv"),sep="\t")
summary <- features[,.(features=.N),by=.(reference_species,gene_biotype)]
jsonlite::write_json(list(status="annotated_pending_actual_cache_comparison",
  feature_count=nrow(features),reference_gff_sha256=digest::digest(file=paths,algo="sha256",serialize=FALSE),
  genome_biotype_counts=summary,actual_umap_features_checked=FALSE),
  file.path(root,"umap","gene_variable_features.annotations.json"),auto_unbox=TRUE,pretty=TRUE)
print(summary)
