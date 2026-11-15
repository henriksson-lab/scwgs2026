#!/usr/bin/env Rscript
# Yanai gene-read-count UMAP, using the normal RNA branch of post_unified.R.
# Normalization -> variable genes -> scaling -> PCA -> UMAP. No k-mer analysis.
# Overrides: YANAI_ROOT, YANAI_ZORN_DIR, YANAI_PLOT_DIR,
# YANAI_MIN_GENE_COUNTS=1000, YANAI_UMAP_SEED=42.
# YANAI_PLOT_ONLY=1 redraws saved coordinates without recomputing the analysis.
# YANAI_PLOT_GROUPS=species,dataset,library,depth selects views to redraw.
scriptFile <- tryCatch(sys.frame(1)$ofile,error=function(e) NULL)
if(is.null(scriptFile)) {
  arg <- grep("^--file=",commandArgs(),value=TRUE)
  if(!length(arg)) stop("Run with Rscript or source()")
  scriptFile <- sub("^--file=","",arg[1])
}
scriptDir <- dirname(normalizePath(scriptFile,mustWork=TRUE))
root <- normalizePath(Sys.getenv("YANAI_ROOT","/husky/henriksson/atrandi/v9_rnaseq_yanai"),mustWork=TRUE)
zornDir <- normalizePath(path.expand(Sys.getenv("YANAI_ZORN_DIR","~/github/zorn")),mustWork=TRUE)
plotDir <- path.expand(Sys.getenv("YANAI_PLOT_DIR",file.path(root,"plots")))
outDir <- file.path(root,"umap")
dir.create(plotDir,recursive=TRUE,showWarnings=FALSE)
dir.create(outDir,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(Seurat);library(ggplot2);library(Matrix)})
for(file in c("zorn.R","bascet_file.R","countmatrix.R")) source(file.path(zornDir,"R",file))
if(requireNamespace("future",quietly=TRUE)) future::plan("sequential")
seed <- as.integer(Sys.getenv("YANAI_UMAP_SEED","42"))
minimum <- as.numeric(Sys.getenv("YANAI_MIN_GENE_COUNTS","1000"))
stopifnot(!is.na(seed),is.finite(minimum),minimum>=0)
set.seed(seed)
StyleYanaiUMAP <- function(p) {
  p + labs(title=NULL,x="RNAseq counts 1",y="RNAseq counts 2") +
    theme(plot.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())
}
SaveYanaiPlots <- function(adata,groups=c("species","dataset","library","depth")) {
  stopifnot(length(groups)>0L,all(groups %in% c("species","dataset","library","depth")))
  for(group in intersect(groups,c("species","dataset","library"))) {
    p <- DimPlot(adata,reduction="gene_umap",group.by=group,pt.size=1,raster=TRUE,
                 raster.dpi=c(800,800),shuffle=TRUE,seed=seed,combine=FALSE)[[1]]
    p <- StyleYanaiUMAP(p)
    if(group=="species") p <- LabelClusters(p,id="species",repel=TRUE,size=4)
    if(group %in% c("species","library")) p <- p+NoLegend()
    for(ext in c("pdf","png")) ggsave(file.path(plotDir,paste0("umap_gene_",group,".",ext)),plot=p,width=9,height=6)
  }
  if(!"depth" %in% groups) return(invisible(adata))
  p <- FeaturePlot(adata,reduction="gene_umap",features="log10_gene_reads",pt.size=1,
                   raster=TRUE,raster.dpi=c(800,800),combine=FALSE)[[1]]
  p <- StyleYanaiUMAP(p)
  for(ext in c("pdf","png")) ggsave(file.path(plotDir,paste0("umap_gene_depth.",ext)),plot=p,width=7,height=6)
}
if(identical(Sys.getenv("YANAI_PLOT_ONLY","0"),"1")) {
  saved <- read.csv(file.path(outDir,"gene_umap.csv"),row.names=1,check.names=FALSE)
  # A small plotting object preserves the original embedding and metadata.
  counts <- Matrix::sparseMatrix(i=integer(),j=integer(),x=numeric(),
    dims=c(1L,nrow(saved)),dimnames=list("plot-placeholder",rownames(saved)))
  adata <- CreateSeuratObject(counts=counts,meta.data=saved)
  coords <- as.matrix(saved[,c("UMAP1","UMAP2")])
  colnames(coords) <- c("geneUMAP_1","geneUMAP_2")
  adata[["gene_umap"]] <- CreateDimReducObject(embeddings=coords,key="geneUMAP_",assay="RNA")
  groups <- strsplit(Sys.getenv("YANAI_PLOT_GROUPS","species,dataset,library,depth"),",",fixed=TRUE)[[1]]
  SaveYanaiPlots(adata,groups=groups)
  message("Yanai: redrew saved UMAPs in ",plotDir)
} else {
manifest <- read.delim(file.path(root,"fastq_manifest.tsv"),stringsAsFactors=FALSE)
mapping <- read.csv(file.path(root,"yanai_gene_count_libraries.csv"),stringsAsFactors=FALSE)
stopifnot(nrow(mapping)==nrow(manifest),setequal(mapping$library_prefix,manifest$library_prefix))
message("Reading gene read counts for all ",nrow(mapping)," independent Yanai libraries")
mat <- ReadBascetCountMatrix(root,"genecount")
adata <- CreateSeuratObject(counts=mat,project="Yanai",min.cells=0,min.features=0)
rm(mat);gc()
libraryIds <- sub("^(D[0-9]+_lib[0-9]+)_.*$","\\1",colnames(adata))
index <- match(libraryIds,manifest$library_prefix)
if(anyNA(index)) stop("Gene-count cell IDs do not map to independent GEO libraries")
adata$library <- libraryIds
adata$dataset <- sub("_lib.*$","",libraryIds)
adata$species <- manifest$scientific_name[index]
adata$geo_sample <- manifest$geo_accession[index]
retained <- adata$nCount_RNA>minimum
# Preserve every valid combination in the original HDF5 matrices. UMAP-only
# depth selection follows the >1,000 count threshold used in post_unified.R.
write.csv(data.frame(cell_id=colnames(adata),library=adata$library,
                     gene_read_counts=adata$nCount_RNA,genes_detected=adata$nFeature_RNA,
                     retained=retained),file.path(outDir,"gene_cell_selection.csv"),row.names=FALSE)
before <- table(adata$library)
adata <- adata[,retained]
if(ncol(adata)<4L) stop("Too few cells above the gene-count threshold")
counts <- GetAssayData(adata,layer="counts")
adata <- adata[Matrix::rowSums(counts)>0,]
# The full-cell selection has been written; release its original cell names.
rm(counts,retained,libraryIds,index);gc()
adata$log10_gene_reads <- log10(1+adata$nCount_RNA)
adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata,selection.method="vst",nfeatures=min(2000L,nrow(adata)))
adata <- ScaleData(adata,features=rownames(adata))
gc() # Release temporary scaling blocks before PCA allocates dense matrices.
npcs <- min(50L,ncol(adata)-1L,length(VariableFeatures(adata))-1L)
if(npcs<2L) stop("Too few variable genes for PCA")
adata <- RunPCA(adata,features=VariableFeatures(adata),npcs=npcs,seed.use=seed)
gc() # Release temporary PCA matrices before nearest-neighbor construction.
ndim <- min(30L,npcs)
adata <- RunUMAP(adata,dims=seq_len(ndim),reduction="pca",
                 reduction.name="gene_umap",seed.use=seed,n.neighbors=min(30L,ncol(adata)-1L))
embedding <- Embeddings(adata,"gene_umap")
stopifnot(nrow(embedding)==ncol(adata),ncol(embedding)==2L,all(is.finite(embedding)))
saveRDS(adata,file.path(root,"cache_adata_gene.RDS"))
write.csv(data.frame(cell_id=rownames(embedding),UMAP1=embedding[,1],UMAP2=embedding[,2],
                     adata@meta.data[rownames(embedding),,drop=FALSE],check.names=FALSE),
            file.path(outDir,"gene_umap.csv"),row.names=FALSE)
after <- table(factor(adata$library,levels=names(before)))
write.csv(data.frame(library=names(before),valid_barcodes=as.integer(before),
                     umap_cells=as.integer(after)),file.path(outDir,"gene_library_cell_counts.csv"),row.names=FALSE)
SaveYanaiPlots(adata)
jsonlite::write_json(list(status="complete",method="gene_read_counts",libraries_processed=nrow(mapping),
                         libraries_in_umap=length(unique(adata$library)),cells=ncol(adata),genes=nrow(adata),
                         minimum_gene_counts=minimum,seed=seed,pca_dimensions=ndim,
                         workflow="NormalizeData/FindVariableFeatures(vst,2000)/ScaleData/RunPCA/RunUMAP",
                         count_unit="primary mapped read alignments, no UMI deduplication"),
                     file.path(outDir,"run_summary.json"),auto_unbox=TRUE,pretty=TRUE)
capture.output(sessionInfo(),file=file.path(outDir,"sessionInfo.txt"))
message("Yanai gene-count UMAP complete: ",plotDir)

}
