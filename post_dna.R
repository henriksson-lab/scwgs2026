#!/usr/bin/env Rscript
# DNA UMAPs following post_unified.R, using the local Zorn checkout.
# Source via post_ling.R or post_gorg.R; DNA_DATASET chooses LING or GORG.
# Dataset-prefixed overrides: ROOT, ZORN_DIR, PLOT_DIR, MIN_KRAKEN,
# MIN_SKETCH_KMERS=40000, UMAP_SEED=42, UMAP_THREADS=20.
# UMAP_METHODS=kraken,infokmer,countsketch selects methods to compute now.
# PLOT_ONLY=1 redraws saved reductions without recomputing or rewriting caches.
# GENUS_LABEL_MIN_CELLS=50 selects genus labels for readable cluster annotation.
# Omit it to resume all methods, reusing caches with matching inputs/settings.
# Count UMAPs: RunTFIDF -> FindTopFeatures(q0) -> RunSVD -> RunUMAP(LSI 1:30).
# Sketch UMAP: BascetLoadCountSketchMatrix -> CountSketchUMAP (cosine).

scriptFile <- tryCatch(sys.frame(1)$ofile,error=function(e) NULL)
if(is.null(scriptFile)) {
  arg <- grep("^--file=",commandArgs(),value=TRUE)
  if(!length(arg)) stop("Run with Rscript or source()")
  scriptFile <- sub("^--file=","",arg[1])
}
scriptDir <- dirname(normalizePath(scriptFile,mustWork=TRUE))
dataset <- toupper(Sys.getenv("DNA_DATASET","LING"))
stopifnot(dataset %in% c("LING","GORG"))
env <- function(suffix,default="") Sys.getenv(paste0(dataset,"_",suffix),default)
label <- if(dataset=="LING") "Ling" else "GORG-Dark"
defaultRoot <- if(dataset=="LING") "/husky/henriksson/atrandi/v9_ling" else "/husky/henriksson/atrandi/v9_wgs_gorg"
root <- normalizePath(env("ROOT",defaultRoot),mustWork=TRUE)
zornDir <- normalizePath(path.expand(env("ZORN_DIR","~/github/zorn")),mustWork=TRUE)
plotDir <- path.expand(env("PLOT_DIR",file.path(root,"plots")))
outDir <- file.path(root,"umap")
dir.create(plotDir,recursive=TRUE,showWarnings=FALSE)
dir.create(outDir,recursive=TRUE,showWarnings=FALSE)

required <- c("Seurat","Signac","uwot","irlba","arrow","rhdf5","taxonomizr",
              "ggplot2","patchwork","data.table","naturalsort","progress","jsonlite","MatrixExtra")
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop("Missing R UMAP dependencies: ",paste(missing,collapse=", "))
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(Matrix)
  library(patchwork)
})
for(file in c("zorn.R","bascet_file.R","countmatrix.R","kraken.R","countsketch.R")) {
  source(file.path(zornDir,"R",file))
}
setwd(root)
seed <- as.integer(env("UMAP_SEED","42"))
threads <- as.integer(env("UMAP_THREADS",env("THREADS","20")))
stopifnot(!is.na(seed),!is.na(threads),threads>=1L)
if(requireNamespace("future",quietly=TRUE)) future::plan("sequential")
options(mc.cores=threads)
set.seed(seed)
# SAG identities are established by their assemblies, not called from raw-read
# depth. Retain every SAG with classified fragments; the raw-WGS threshold
# would discard most marine SAGs because reference coverage is sparse.
minKraken <- as.numeric(env("MIN_KRAKEN",if(dataset=="GORG") "0" else "1000"))
minSketch <- as.numeric(env("MIN_SKETCH_KMERS","40000"))
stopifnot(is.finite(minKraken),minKraken>=0,is.finite(minSketch),minSketch>=0)

AddLingMetadata <- function(adata) {
  if(dataset=="GORG") {
    adata$library <- "GORG-Dark"
    adata$sample <- "GORG-Dark"
    adata$sequencing <- "SAG assemblies"
    return(adata)
  }
  ids <- colnames(adata)
  library <- sub("^((feces|sewage)_(deep|shallow))_.*$","\\1",ids)
  if(!all(library %in% c("feces_deep","feces_shallow","sewage_deep","sewage_shallow"))) {
    stop("Unexpected Ling cell ID prefix")
  }
  adata$library <- library
  adata$sample <- sub("_.*$","",library)
  adata$sequencing <- sub("^.*_","",library)
  adata
}

SaveSelection <- function(method,allCells,keptCells,rule) {
  write.csv(data.frame(cell_id=allCells,retained=allCells %in% keptCells,rule=rule),
            file.path(outDir,paste0(method,"_cell_selection.csv")),row.names=FALSE)
}

StyleDNAUMAP <- function(p,method) {
  prefix <- c(kraken="KRAKEN",infokmer="IK",cs="CS")[[method]]
  p + labs(title=NULL,x=paste0(prefix,"1"),y=paste0(prefix,"2")) +
    theme(plot.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())
}

SavePlots <- function(adata,reduction,method,depthFeature,saveData=TRUE) {
  # DimPlot/FeaturePlot as in post_unified.R, with Ling library labels.
  pointSize <- if(ncol(adata)<1000L) 2 else 1
  for(group in c("library","sample","sequencing","phylum","genus","species")) {
    if(!group %in% colnames(adata@meta.data)) next
    p <- DimPlot(adata,reduction=reduction,group.by=group,pt.size=pointSize,
                 raster=TRUE,raster.dpi=c(800,800),shuffle=TRUE,seed=seed,
                 combine=FALSE)[[1]]
    if(group=="genus") {
      minLabels <- as.integer(env("GENUS_LABEL_MIN_CELLS","50"))
      stopifnot(!is.na(minLabels),minLabels>=0L)
      genusCounts <- table(adata$genus)
      genera <- names(genusCounts)[genusCounts>=minLabels &
                                     !names(genusCounts) %in% c("","NA","Unclassified","unclassified")]
      if(length(genera)) p <- LabelClusters(p,id="genus",clusters=genera,repel=TRUE,size=3)
    }
    p <- StyleDNAUMAP(p,method)
    if(group %in% c("phylum","genus","species")) p <- p + NoLegend()
    for(ext in c("pdf","png")) {
      ggsave(file.path(plotDir,paste0("umap_",method,"_",group,".",ext)),
             plot=p,width=9,height=6,limitsize=FALSE)
    }
  }
  p <- FeaturePlot(adata,reduction=reduction,features=depthFeature,pt.size=pointSize,
                   raster=TRUE,raster.dpi=c(800,800),combine=FALSE)[[1]]
  p <- StyleDNAUMAP(p,method)
  for(ext in c("pdf","png")) {
    ggsave(file.path(plotDir,paste0("umap_",method,"_depth.",ext)),plot=p,width=7,height=6)
  }
  if(!saveData) return(invisible(adata))
  embedding <- Embeddings(adata,reduction)
  stopifnot(ncol(embedding)==2L,all(is.finite(embedding)),!anyDuplicated(rownames(embedding)))
  write.csv(data.frame(cell_id=rownames(embedding),UMAP1=embedding[,1],UMAP2=embedding[,2],
                       adata@meta.data[rownames(embedding),,drop=FALSE],check.names=FALSE),
            file.path(outDir,paste0(method,"_umap.csv")),row.names=FALSE)
  saveRDS(adata,file.path(root,paste0("cache_adata_",method,".RDS")))
  invisible(adata)
}

RunCountUMAP <- function(adata,reduction,method) {
  # Keep the same operations/defaults as post_unified.R. Adapt dimension count
  # only for small validation datasets where 30/50 dimensions are impossible.
  counts <- GetAssayData(adata,layer="counts")
  keepFeatures <- Matrix::rowSums(counts)>0
  keepCells <- Matrix::colSums(counts)>0
  adata <- adata[keepFeatures,keepCells]
  if(ncol(adata)<4L || nrow(adata)<3L) stop("Too few informative cells/features for ",method)
  adata <- RunTFIDF(adata)
  adata <- FindTopFeatures(adata,min.cutoff="q0")
  nsvd <- min(50L,ncol(adata)-1L,length(VariableFeatures(adata))-1L)
  if(nsvd<2L) stop("Too few variable features for ",method)
  adata <- RunSVD(adata,n=nsvd,seed.use=seed)
  p <- DepthCor(adata)
  ggsave(file.path(plotDir,paste0(method,"_depthcor.pdf")),plot=p,width=4,height=3)
  ndim <- min(30L,ncol(Embeddings(adata,"lsi")))
  adata <- RunUMAP(adata,reduction="lsi",dims=seq_len(ndim),
                   reduction.name=reduction,seed.use=seed,
                   n.neighbors=min(30L,ncol(adata)-1L))
  adata
}

allMethods <- c("kraken","infokmer","countsketch")
methodCache <- c(kraken="kraken",infokmer="infokmer",countsketch="cs")
requested <- env("UMAP_METHODS","")
methods <- if(nzchar(requested)) strsplit(requested,",",fixed=TRUE)[[1]] else allMethods
stopifnot(length(methods)>0L,all(methods %in% allMethods))
if(identical(env("PLOT_ONLY","0"),"1")) {
  reductions <- c(kraken="kraken_umap",infokmer="infokmers_umap",countsketch="kmersketch_umap")
  depths <- c(kraken="log_cnt",infokmer="log10_nCount_infokmer",countsketch="log10_celldepth")
  for(method in methods) {
    cache <- file.path(root,paste0("cache_adata_",methodCache[[method]],".RDS"))
    if(!file.exists(cache)) stop("No saved UMAP object: ",cache)
    adata <- readRDS(cache)
    SavePlots(adata,reductions[[method]],methodCache[[method]],depths[[method]],saveData=FALSE)
    rm(adata);gc()
  }
  message(label,": redrew cached UMAPs in ",plotDir)
} else {
summaryFile <- file.path(outDir,"run_summary.json")
summary <- if(file.exists(summaryFile)) jsonlite::read_json(summaryFile) else list()
if(is.null(summary$cache_inputs)) summary$cache_inputs <- list()
shards <- if(dataset=="LING") 1:4 else 1L
sourceSignature <- as.character(tools::md5sum(c(file.path(scriptDir,"post_dna.R"),
  file.path(zornDir,"R",c("zorn.R","countmatrix.R","kraken.R","countsketch.R")))))
CacheFiles <- function(method) {
  files <- c(file.path(root,paste0("cache_adata_",methodCache[[method]],".RDS")),
             file.path(outDir,paste0(methodCache[[method]],"_umap.csv")),
             file.path(outDir,paste0(method,"_cell_selection.csv")))
  groups <- c("library","sample","sequencing","phylum","genus","species","depth")
  files <- c(files,as.vector(outer(file.path(plotDir,paste0("umap_",methodCache[[method]],"_",groups)),
                                  c("pdf","png"),paste,sep=".")))
  if(method=="kraken") files <- c(files,file.path(root,c("kraken_metadata.RDS","adata_kraken_metadata.RDS")))
  files
}
InputSignature <- function(method) {
  files <- c(file.path(root,paste0("kraken_counts.",shards,".h5")),zornTaxonomizrSqlFile())
  if(method=="infokmer") files <- c(files,file.path(root,paste0("kmer_counts.",shards,".h5")),
                                    file.path(root,"use_kmers.txt"))
  if(method=="countsketch") files <- c(files,file.path(root,"countsketch_mat.feather"))
  if(!all(file.exists(files))) return(NULL)
  info <- file.info(files)
  as.character(jsonlite::toJSON(list(paths=normalizePath(files),bytes=info$size,
    modified=as.numeric(info$mtime),seed=seed,min_kraken=minKraken,min_sketch=minSketch,
    source=sourceSignature,Seurat=as.character(packageVersion("Seurat")),
    Signac=as.character(packageVersion("Signac")),uwot=as.character(packageVersion("uwot"))),
    auto_unbox=TRUE,digits=15))
}
ShouldRun <- function(method) {
  if(!method %in% methods) return(FALSE)
  signature <- InputSignature(method)
  reusable <- !nzchar(requested) && !is.null(signature) && !is.null(summary[[method]]) &&
    identical(summary$cache_inputs[[method]],signature) && all(file.exists(CacheFiles(method)))
  if(reusable) message(label,": reusing saved ",method," UMAP (inputs/settings unchanged)")
  !reusable
}
SaveProgress <- function() {
  summary$seed <<- seed
  summary$dataset <<- dataset
  summary$plot_directory <<- plotDir
  summary$coordinate_directory <<- outDir
  complete <- all(vapply(allMethods,function(method) {
    !is.null(summary[[method]]) && all(file.exists(CacheFiles(method)))
  },logical(1)))
  summary$status <<- if(complete) "complete" else "partial"
  tmp <- tempfile(".run_summary.",tmpdir=outDir)
  jsonlite::write_json(summary,tmp,auto_unbox=TRUE,pretty=TRUE)
  stopifnot(file.rename(tmp,summaryFile))
}
if(ShouldRun("kraken")) {
message(label," Kraken UMAP: reading counts and consensus taxonomy")
mat <- ReadBascetCountMatrix(root,"kraken_counts",verbose=FALSE)
allKrakenCells <- rownames(mat)
kraken_taxid <- KrakenFindConsensusTaxonomy(mat)
saveRDS(kraken_taxid,file.path(root,"kraken_metadata.RDS"))
adata <- CreateSeuratObject(counts=mat,project=label,min.cells=0,min.features=0)
adata <- BascetAddMetaData(adata,kraken_taxid)
adata <- AddLingMetadata(adata)
saveRDS(adata@meta.data,file.path(root,"adata_kraken_metadata.RDS"))
adata <- adata[,adata$nCount_RNA>minKraken]
adata$log_cnt <- log10(1+adata$nCount_RNA)
adata <- RunCountUMAP(adata,"kraken_umap","kraken")
SaveSelection("kraken",allKrakenCells,colnames(adata),paste0("nCount_RNA > ",minKraken,"; nonempty features"))
SavePlots(adata,"kraken_umap","kraken","log_cnt")
summary$kraken <- list(cells=ncol(adata),features=nrow(adata),min_counts=minKraken,
                       workflow="RunTFIDF/FindTopFeatures(q0)/RunSVD/RunUMAP(lsi,1:30)")
summary$cache_inputs$kraken <- InputSignature("kraken")
SaveProgress()
rm(mat,adata);gc()
} else kraken_taxid <- readRDS(file.path(root,"kraken_metadata.RDS"))

if(ShouldRun("infokmer")) {
message(label," informative 31-mer UMAP")
cnt <- ReadBascetCountMatrix(root,"kmer_counts")
cnt <- cnt[order(rowSums(cnt),decreasing=TRUE),]
allInfoCells <- rownames(cnt)
adata <- CreateSeuratObject(counts=cnt,assay="infokmer",min.cells=0,min.features=0)
adata <- BascetAddMetaData(adata,kraken_taxid,subsetCommon=TRUE)
adata <- AddLingMetadata(adata)
adata$log10_nCount_infokmer <- log10(1+adata$nCount_infokmer)
adata <- RunCountUMAP(adata,"infokmers_umap","infokmer")
SaveSelection("infokmer",allInfoCells,colnames(adata),"Cells with Kraken metadata and nonzero informative 31-mer counts")
SavePlots(adata,"infokmers_umap","infokmer","log10_nCount_infokmer")
summary$infokmer <- list(cells=ncol(adata),features=nrow(adata),
                         workflow="RunTFIDF/FindTopFeatures(q0)/RunSVD/RunUMAP(lsi,1:30)")
summary$cache_inputs$infokmer <- InputSignature("infokmer")
SaveProgress()
rm(cnt,adata);gc()
}

if(ShouldRun("countsketch")) {
message(label," CountSketch UMAP")
adata <- BascetLoadCountSketchMatrix(root,inputName="countsketch_mat.feather")
allSketchCells <- colnames(adata)
adata$log10_celldepth <- log10(1+adata$celldepth)
adata <- adata[,adata$celldepth>minSketch]
adata <- BascetAddMetaData(adata,kraken_taxid,subsetCommon=TRUE)
adata <- AddLingMetadata(adata)
if(ncol(adata)<4L) stop("Too few cells above CountSketch depth threshold")
adata <- CountSketchUMAP(adata,n_neighbors=min(30L,ncol(adata)-1L),seed=seed)
SaveSelection("countsketch",allSketchCells,colnames(adata),paste0("celldepth > ",minSketch," and Kraken metadata"))
SavePlots(adata,"kmersketch_umap","cs","log10_celldepth")
summary$countsketch <- list(cells=ncol(adata),dimensions=ncol(Embeddings(adata,"kmersketch")),
                           min_kmers=minSketch,workflow="BascetLoadCountSketchMatrix/CountSketchUMAP(cosine)")
summary$cache_inputs$countsketch <- InputSignature("countsketch")
SaveProgress()
}

# Show available reductions during a partial run; write the standard overview
# once all three methods are complete. Existing reductions are retained.
available <- allMethods[vapply(allMethods,function(method) {
  !is.null(summary[[method]]) && all(file.exists(CacheFiles(method)))
},logical(1))]
objects <- lapply(available,function(method)readRDS(file.path(root,paste0("cache_adata_",methodCache[[method]],".RDS"))))
reductions <- c(kraken="kraken_umap",infokmer="infokmers_umap",countsketch="kmersketch_umap")
panels <- lapply(seq_along(objects),function(i) {
  p <- DimPlot(objects[[i]],reduction=reductions[[available[i]]],group.by="library",
               pt.size=if(ncol(objects[[i]])<1000L) 2 else 1,raster=TRUE,raster.dpi=c(800,800))
  StyleDNAUMAP(p,methodCache[[available[i]]])
})
p <- wrap_plots(panels,nrow=1,guides="collect")
overview <- if(length(available)==3L) "all_umaps_library" else "available_umaps_library"
for(ext in c("pdf","png")) ggsave(file.path(plotDir,paste0(overview,".",ext)),plot=p,width=6*length(available),height=5)
SaveProgress()
capture.output(sessionInfo(),file=file.path(outDir,"sessionInfo.txt"))
message(label," UMAP status: ",summary$status,"; saved plots: ",plotDir)
}
