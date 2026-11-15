#!/usr/bin/env Rscript
# Full-dataset completion audit. Run after run_{gorg,ling,yanai}_all.sh.
# Optional positional arguments: gorg ling yanai (default: all three).
suppressPackageStartupMessages(library(Seurat))
scriptArg <- grep("^--file=",commandArgs(),value=TRUE)
scriptDir <- dirname(normalizePath(sub("^--file=","",scriptArg[1]),mustWork=TRUE))
requested <- commandArgs(trailingOnly=TRUE)
if(!length(requested)) requested <- c("gorg","ling","yanai")
stopifnot(all(requested %in% c("gorg","ling","yanai")))
defaults <- c(gorg="/husky/henriksson/atrandi/v9_wgs_gorg",
              ling="/husky/henriksson/atrandi/v9_ling",
              yanai="/husky/henriksson/atrandi/v9_rnaseq_yanai")
requireFiles <- function(paths) {
  missing <- paths[!file.exists(paths) | is.na(file.info(paths)$size) | file.info(paths)$size==0]
  if(length(missing)) stop("Missing/empty artifacts: ",paste(missing,collapse=", "))
}
audit <- function(dataset) {
  root <- Sys.getenv(paste0(toupper(dataset),"_ROOT"),defaults[[dataset]])
  out <- file.path(root,"umap")
  requireFiles(file.path(out,c("run_summary.json","sessionInfo.txt")))
  summary <- jsonlite::read_json(file.path(out,"run_summary.json"))
  stopifnot(identical(summary$status,"complete"))
  if(dataset=="gorg") stopifnot(summary$kraken$min_counts==0)
  plotDir <- Sys.getenv(paste0(toupper(dataset),"_PLOT_DIR"),file.path(root,"plots"))
  methods <- if(dataset=="yanai") "gene" else c("kraken","infokmer","cs")
  if(dataset!="yanai") {
    shards <- if(dataset=="ling") 1:4 else 1L
    requireFiles(c(file.path(root,"countsketch_mat.feather"),
                   file.path(root,paste0("kraken_counts.",shards,".h5")),
                   file.path(root,paste0("kmer_counts.",shards,".h5"))))
  }
  reductions <- c(gene="gene_umap",kraken="kraken_umap",
                  infokmer="infokmers_umap",cs="kmersketch_umap")
  checked <- list()
  if(dataset!="yanai") {
    kmers <- readLines(file.path(root,"use_kmers.txt"))
    stopifnot(length(kmers)==100000L,!anyDuplicated(kmers),all(grepl("^[ACGT]{31}$",kmers)))
    checked$informative_kmer_pool <- list(features=length(kmers),k=31L)
    checked$informative_inputs <- lapply(shards,function(shard) {
      file <- file.path(root,paste0("kmer_counts.",shard,".h5"))
      cells <- as.character(rhdf5::h5read(file,"/obs/_index"))
      features <- as.character(rhdf5::h5read(file,"/var/_index"))
      originalCells <- as.character(rhdf5::h5read(
        file.path(root,paste0("kraken_counts.",shard,".h5")),"/obs/_index"))
      shape <- as.integer(rhdf5::h5readAttributes(file,"/X")$shape)
      stopifnot(!anyDuplicated(cells),setequal(cells,originalCells),
                !anyDuplicated(features),setequal(features,kmers),
                identical(shape,c(length(cells),100000L)))
      list(shard=shard,cells=length(cells),features=length(features),
           exact_barcode_coverage=TRUE,exact_feature_pool=TRUE)
    })
  }
  for(method in methods) {
    cache <- file.path(root,paste0("cache_adata_",method,".RDS"))
    coordinates <- file.path(out,paste0(method,"_umap.csv"))
    requireFiles(c(cache,coordinates))
    object <- readRDS(cache)
    xy <- Embeddings(object,reductions[[method]])
    table <- read.csv(coordinates,check.names=FALSE)
    stopifnot(nrow(xy)>=4L,ncol(xy)==2L,all(is.finite(xy)),
              !anyDuplicated(rownames(xy)),identical(rownames(xy),colnames(object)),
              identical(table$cell_id,rownames(xy)),
              isTRUE(all.equal(unname(as.matrix(table[,c("UMAP1","UMAP2")])),
                               unname(xy),tolerance=1e-10,check.attributes=FALSE)))
    selectionMethod <- if(method=="cs") "countsketch" else method
    selection <- read.csv(file.path(out,paste0(selectionMethod,"_cell_selection.csv")))
    stopifnot(!anyDuplicated(selection$cell_id),
              setequal(selection$cell_id[selection$retained],rownames(xy)))
    groups <- if(dataset=="yanai") c("species","dataset","library","depth") else
      c("library","sample","sequencing","phylum","genus","species","depth")
    requireFiles(as.vector(outer(file.path(plotDir,paste0("umap_",method,"_",groups)),
                                 c("pdf","png"),paste,sep=".")))
    recordedCells <- if(dataset=="yanai") summary$cells else
      summary[[selectionMethod]]$cells
    stopifnot(recordedCells==nrow(xy))
    if(dataset=="gorg" && method %in% c("kraken","cs")) {
      stopifnot(nrow(selection)==9698L)
      if(method=="kraken" && summary$kraken$min_counts==0) stopifnot(nrow(xy)==9698L)
    }
    checked[[method]] <- list(cells=nrow(xy),input_cells=nrow(selection),
                              libraries=sort(unique(object$library)))
    if(dataset=="yanai") {
      source(file.path(scriptDir,"yanai_petriseq","verify_cache_native_counts.R"))
      checked$gene_cached_counts <- VerifyYanaiCacheNativeCounts(root,object)
      source(file.path(scriptDir,"yanai_petriseq","verify_variable_gene_coverage.R"))
      checked <- c(checked,VerifyYanaiVariableGeneCoverage(root,object,scriptDir))
    }
    rm(object,xy,table,selection);gc()
  }
  if(dataset=="ling") {
    receipts <- file.path(root,paste0("debarcoded.",1:32,".tirp.gz.json"))
    requireFiles(receipts)
    for(receipt in receipts) {
      record <- jsonlite::read_json(receipt)
      tirp <- sub("[.]json$","",receipt)
      histogram <- paste0(tirp,".hist")
      requireFiles(c(tirp,histogram))
      stopifnot(file.info(tirp)$size==record$tirp_bytes,
                file.info(histogram)$size==record$hist_bytes,
                record$records>0,record$cells>0)
    }
    meta <- read.csv(file.path(root,"debarcoded.meta"))
    stopifnot(nrow(meta)==32L,length(unique(meta$prefix))==4L)
    sharding <- readRDS(file.path(root,"ling_sharding.rds"))
    expectedCells <- unlist(sharding$list_picked_cells,use.names=FALSE)
    sketchSelection <- read.csv(file.path(out,"countsketch_cell_selection.csv"))
    stopifnot(!anyDuplicated(expectedCells),setequal(expectedCells,sketchSelection$cell_id))
    source(file.path(scriptDir,"ling_data","verify_kraken.R"))
    checked$kraken_input <- VerifyLingKrakenInput(root,sharding)
    source(file.path(scriptDir,"ling_data","verify_merge.R"))
    checked$merged_kraken <- VerifyLingMergedKraken(root)
    requireFiles(file.path(root,paste0("filtered.",1:4,".tirp.gz")))
    for(method in methods) stopifnot(length(checked[[method]]$libraries)==4L)
  } else if(dataset=="yanai") {
    source(file.path(scriptDir,"yanai_petriseq","verify_shards.R"))
    checked$shard_headers <- VerifyYanaiShards(root)
    manifest <- read.delim(file.path(root,"fastq_manifest.tsv"))
    mapping <- read.csv(file.path(root,"yanai_gene_count_libraries.csv"))
    stopifnot(nrow(manifest)==88L,nrow(mapping)==88L,
              !anyDuplicated(mapping$library_prefix),
              setequal(manifest$library_prefix,mapping$library_prefix),
              summary$libraries_processed==88L)
    requireFiles(c(file.path(root,mapping$counts),file.path(root,mapping$alignment),
                   file.path(root,paste0(mapping$alignment,".bai"))))
    selection <- read.csv(file.path(out,"gene_cell_selection.csv"))
    stopifnot(nrow(selection)==sum(mapping$retained_cells),
              setequal(unique(selection$library),manifest$library_prefix))
    for(i in seq_len(nrow(mapping))) {
      stopifnot(sum(selection$library==mapping$library_prefix[i])==mapping$retained_cells[i])
    }
    source(file.path(scriptDir,"yanai_petriseq","verify_gene_counts.R"))
    checked$gene_inputs <- VerifyYanaiGeneCounts(root,selection)
    # A library can have no cells above the declared UMAP-only depth threshold.
    # Report that explicitly; every valid barcode must still be counted above.
    checked$processed_libraries <- mapping$library_prefix
    checked$libraries_below_umap_threshold <- setdiff(mapping$library_prefix,checked$gene$libraries)
  }
  report <- list(status="verified",dataset=dataset,root=root,
                 checked_at=format(Sys.time(),tz="UTC",usetz=TRUE),methods=checked)
  jsonlite::write_json(report,file.path(out,"validation.json"),auto_unbox=TRUE,pretty=TRUE)
  message(dataset,": full-dataset UMAP artifacts verified")
}
failed <- character()
for(dataset in requested) tryCatch(audit(dataset),error=function(e) {
  message(dataset,": INCOMPLETE: ",conditionMessage(e))
  failed <<- c(failed,dataset)
})
if(length(failed)) quit(status=1L)
