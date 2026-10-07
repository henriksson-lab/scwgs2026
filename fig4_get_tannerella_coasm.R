library(Zorn)
library(Seurat)
library(Signac)
library(sqldf)
library(stringr)
library(ggplot2)
library(Matrix)
library(GenomicRanges)

library(SingleCellExperiment)
library(scDblFinder)

library(future)
plan("multicore", workers = 10)


bascet_runner <- LocalRunner(direct = TRUE)

dataset_name <- "v9_wgs_saliva2_mda"
bascetRoot <- file.path("/husky/henriksson/atrandi/",dataset_name)

plotDirAll <- file.path("~/github/scmetag-revision/plots/all")
if(!file.exists(plotDirAll)){
  dir.create(plotDirAll, recursive = TRUE)
  print("made plot dir all")
}


setwd("/home/mahogny/github/zorn")


################################################################################
############### Alignment-based analysis (tannerella only) #####################
################################################################################


##########
########## Read the aligned counts
##########
cnt <- ReadBascetCountMatrix(bascetRoot,"countchrom_tannerella", verbose=FALSE)

adata <- CreateSeuratObject(
  counts = cnt,
  project = "proj", min.cells = 0, min.features = 0, assay = "chrom_cnt"
) 


min_nCount_RNA <- 1
## Dimensional reduction using kraken
sum(adata$nCount_chrom_cnt>min_nCount_RNA)
keep_cells <- adata$nCount_chrom_cnt > min_nCount_RNA
#keep_cells <- adata$nCount_chrom_cnt > sort(adata$nCount_chrom_cnt, decreasing = TRUE)[10000]
sum(keep_cells)
adata <- adata[,keep_cells] ## Reduce to sensible number
adata


#Figure out which species has most reads in which cell
adata$log_cnt <- log10(1+adata$nCount_chrom_cnt)


##########
########## Generate UMAP based on alignment
##########

if(TRUE){
  #ATAC-seq style. think not the best way here
  adata <- RunTFIDF(adata)
  adata <- FindTopFeatures(adata, min.cutoff = 'q0')
  ndim <- min(50, nrow(adata)-1)
  adata <- RunSVD(adata, n = ndim)
  DepthCor(adata)
  #ggsave(file.path("~/github/scmetag-revision/plots/umap_mock_atac_depthcor.png"), width = 7, height = 5)
  adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 1:ndim, reduction.name = "chrom_umap")  ## depth seems to be less of a problem here
} else {
  #RNA-seq style
  adata <- NormalizeData(adata)
  adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = nrow(adata))
  adata <- ScaleData(adata, features = rownames(adata))
  adata <- RunPCA(adata, features = VariableFeatures(object = adata))
  adata <- RunUMAP(adata, dims = 1:min(50, nrow(adata)-2), reduction.name = "chrom_umap")  #  why -2 here?
}

saveRDS(adata, "/husky/henriksson/atrandi/v9_wgs_saliva2_mda/seurat_tannerella.RDS")

adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata, selection.method = "vst")
adata <- ScaleData(adata)
adata <- RunPCA(adata, features = VariableFeatures(object = adata))

adata <- FindNeighbors(adata, dims = 1:10)
adata <- FindClusters(adata, resolution = 0.5)
DimPlot(object = adata, label = TRUE)


##########
########## Plot tanerella count UMAPs
##########

FeaturePlot(adata, features = "log_cnt")

#kneeplot
plot(sort(adata@meta.data$log_cnt))




################################################################################
################## Get tanerella genomes for co-assembly #######################
################################################################################


beagle_filename_base <- "/husky/henriksson/atrandi/v9_wgs_saliva2_mda/aligned_justtannerella_tirp"

find_cell_in_shards <- function() {
  #### We should have this command for ZIP and TIRP. so add TIRP support in the same place?

  alldf <- list()  
  for(i in 1:20) { ### detect num shards TODO
    fname <- paste0(beagle_filename_base ,".",i,".tirp.gz") ################# tannerella reads
    cmd <- paste("tabix ",fname,"--list-chroms")
    print(cmd)
    df_one <- data.frame(
      shard = i,
      shardname = fname,
      cellid = system(cmd, intern = TRUE)
    )
    alldf[[i]] <- df_one
  }
  do.call(rbind, alldf)
}


extract_fq <- function(fname_tabix, cellid, fname_fq_r1, fname_fq_r2) {
  cell_fasta <- system(paste("tabix ",fname_tabix, cellid), intern = TRUE)
  cell_fasta <- stringr::str_split_fixed(cell_fasta,"\t",8)
  
  r1 <- cell_fasta[,4]
  q1 <- cell_fasta[,6]

  r2 <- cell_fasta[,5]
  q2 <- cell_fasta[,7]
  
  ######### R1
  tou <- data.frame(
    "@someread",
    r1,
    "+",
    q1
  )
#  print(dim(tou))
  
  writeLines(
    as.character(t(as.matrix(tou))),
    fname_fq_r1
  )  
  
  ######### R2
  tou <- data.frame(
    "@someread",
    r2,
    "+",
    q2
  )
  #  print(dim(tou))
  
  writeLines(
    as.character(t(as.matrix(tou))),
    fname_fq_r2
  )  
  
}



extract_all <- function(list_get_cells, outdir_fq) {
#  rawread_pos <- find_cell_in_shards()
#  list_get_cells <- colnames(adata)[adata$log_cnt>3]
  rawread_pos <- rawread_pos[rawread_pos$cellid %in% list_get_cells,]
  dir.create(outdir_fq)
  for(i in 1:nrow(rawread_pos)) {
    print(i)
    extract_fq(
      rawread_pos$shardname[i],
      rawread_pos$cellid[i],
      file.path(outdir_fq, paste0("cell__",rawread_pos$cellid[i],".R1.fq")),
      file.path(outdir_fq, paste0("cell__",rawread_pos$cellid[i],".R2.fq"))
    )
  }
}


############## Extract all cells with >10k tannerella reads for coassembly

rawread_pos <- find_cell_in_shards()
list_get_cells <- colnames(adata)[adata$log_cnt>4]
outdir_fq <- "/husky/henriksson/atrandi/v9_wgs_saliva2_mda/outset_log4"
extract_all(list_get_cells, outdir_fq)

