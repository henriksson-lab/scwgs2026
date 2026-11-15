#!/usr/bin/env Rscript
# Integration fixture for the standard R/Zorn Ling UMAP workflow.
args <- grep('^--file=',commandArgs(),value=TRUE)
stopifnot(length(args)==1L)
scriptDir <- dirname(normalizePath(sub('^--file=','',args[1]),mustWork=TRUE))
suppressPackageStartupMessages({library(Matrix);library(rhdf5);library(arrow)})
root <- tempfile('ling-r-umap-validation-')
dir.create(root,showWarnings=FALSE)
set.seed(42)
libraries <- c('feces_deep','feces_shallow','sewage_deep','sewage_shallow')
cells <- unlist(lapply(libraries,function(lib)paste0(lib,'_ACTGACTG_',sprintf('%024d',1:10))))
taxids <- paste0('taxid_',c(562,1280,1423,1351,1301,287,573,28901,199310,1496,1678,851))
expectedKraken <- matrix(0,nrow=length(cells),ncol=length(taxids),dimnames=list(cells,taxids))
expectedKmer <- matrix(0,nrow=length(cells),ncol=60,dimnames=list(cells,paste0('kmer',1:60)))
WriteCounts <- function(path,x,ids,features) {
  if(file.exists(path)) unlink(path)
  stopifnot(nrow(x)==length(ids),ncol(x)==length(features))
  x <- as(as(Matrix(x,sparse=TRUE),'RsparseMatrix'),'dgRMatrix')
  h5createFile(path)
  for(group in c('obs','var','X')) h5createGroup(path,group)
  h5write(ids,path,'obs/_index');h5write(features,path,'var/_index')
  h5write(rep(5L,length(ids)),path,'obs/_unmapped')
  h5write(x@j,path,'X/indices');h5write(x@p,path,'X/indptr');h5write(x@x,path,'X/data')
  h5write(as.integer(dim(x)),path,'X/shape')
}
for(i in seq_along(libraries)) {
  ids <- cells[((i-1)*10+1):(i*10)]
  x <- matrix(rpois(10*length(taxids),500),nrow=10)
  x[,i] <- x[,i]+rpois(10,1000)
  # Real libraries have different feature sets; exercise the sparse union path.
  keep <- setdiff(seq_along(taxids),i)
  expectedKraken[ids,taxids[keep]] <- x[,keep]
  WriteCounts(file.path(root,paste0('kraken_counts.',i,'.h5')),x[,keep],ids,taxids[keep])
  x <- matrix(rpois(10*60,15),nrow=10)
  x[,((i-1)*15+1):(i*15)] <- x[,((i-1)*15+1):(i*15)]+10
  keep <- setdiff(1:60,i)
  expectedKmer[ids,paste0('kmer',keep)] <- x[,keep]
  WriteCounts(file.path(root,paste0('kmer_counts.',i,'.h5')),x[,keep],ids,paste0('kmer',keep))
}
writeLines(colnames(expectedKmer),file.path(root,'use_kmers.txt'))
cs <- data.frame(cell_id=cells,depth=rep(100000,40))
for(i in 0:63) cs[[paste0('cs_',i)]] <- rnorm(40)
arrow::write_feather(cs,file.path(root,'countsketch_mat.feather'))
Sys.setenv(LING_ROOT=root,LING_PLOT_DIR=file.path(root,'plots'),LING_UMAP_THREADS='2')
source(file.path(scriptDir,'..','post_ling.R'))
for(method in c('kraken','kmer')) {
  expected <- if(method=='kraken') expectedKraken else expectedKmer
  actual <- ReadBascetCountMatrix(root,paste0(method,'_counts'))@X
  stopifnot(setequal(rownames(actual),rownames(expected)),setequal(colnames(actual),colnames(expected)),
            all(as.matrix(actual[rownames(expected),colnames(expected)])==expected))
}
for(method in c('kraken','infokmer','cs')) {
  adata <- readRDS(file.path(root,paste0('cache_adata_',method,'.RDS')))
  stopifnot(ncol(adata)==40L,setequal(adata$library,libraries))
  coordinates <- read.csv(file.path(root,'umap',paste0(method,'_umap.csv')))
  stopifnot(nrow(coordinates)==40L,all(is.finite(as.matrix(coordinates[,c('UMAP1','UMAP2')]))))
}
stopifnot(file.exists(file.path(root,'plots','all_umaps_library.pdf')))
cat('All three R/Zorn Ling UMAPs passed on 40-cell, four-library fixture\n')
