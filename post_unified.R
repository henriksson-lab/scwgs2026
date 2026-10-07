### preprocessing of data done in separate file, prep_*.R
if(FALSE) {
  BiocManager::install("scDblFinder")
}

library(Zorn)

list_color_mock <- c(khroma::color("muted")(n = 9),"#DDDDDD")

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
#plan("multiprocess", workers = 20)


bascet_runner <- LocalRunner(direct = TRUE)

#pick any of these
if(FALSE) {
  dataset_name <- "v9_wgs_mock1_mda"
  dataset_name <- "v9_wgs_mock2_mda"
  dataset_name <- "v9_wgs_mock3_mda"
  dataset_name <- "v9_wgs_saliva1_mda"
  
  dataset_name <- "v9_wgs_chaitra2"
  
  
  dataset_name <- "v9_wgs_mock4_mda"
  dataset_name <- "v9_wgs_mock5_pta"
  dataset_name <- "v9_simulated5"
  dataset_name <- "v9_wgs_saliva2_mda"
  dataset_name <- "zheng2022"
}
dataset_name <- "v9_wgs_mock5_pta"

bascetRoot <- file.path("/husky/henriksson/atrandi/",dataset_name)

plotDir <- file.path("~/github/scmetag-revision/plots/",dataset_name)
if(!file.exists(plotDir)){
  dir.create(plotDir, recursive = TRUE)
  print("made plot dir")
}

plotDirAll <- file.path("~/github/scmetag-revision/plots/all")
if(!file.exists(plotDirAll)){
  dir.create(plotDirAll, recursive = TRUE)
  print("made plot dir all")
}


setwd("/home/mahogny/github/zorn") #to put SQL in the right place

## Mock community: which chromosome is which strain
if(stringr::str_detect(dataset_name, "sim")) {
  #Cannot use ATCC public genomes for sim
  mapSeq2strain <- read.csv("~/github/scmetag-revision/map_seq2strain.sim.csv")
  mapSeq2strain$id <- stringr::str_replace_all(mapSeq2strain$id,"_","-")
  strain_genomesize <- sqldf::sqldf("select 0 as len, strain from mapSeq2strain group by strain") ### clean up names more?
  mapSeq2strain$strain <- stringr::str_split_i(mapSeq2strain$strain,stringr::fixed(" ("),1)
} else {
  #ATCC for real data
  mapSeq2strain <- read.csv("/husky/henriksson/atrandi/bwa_ref/atcc/map_seq2strain.csv")
  mapSeq2strain$id <- stringr::str_replace_all(mapSeq2strain$id,"_","-")
  strain_genomesize <- sqldf::sqldf("select 0 as len, strain from mapSeq2strain group by strain")
}






################################################################################
################## Alignment-based analysis (mock only) ########################
################################################################################

if(str_detect(bascetRoot,"mock")){
  
  ##########
  ########## Read the aligned counts
  ##########
  
  cnt <- ReadBascetCountMatrix(bascetRoot,"chrom_counts", verbose=FALSE)
  
  #Remove phiX (do we need to?)
  if(TRUE) {
    cnt <- cnt[,colnames(cnt)!="phix"] 
  }
  
  
  adata <- CreateSeuratObject(
    counts = cnt,
    project = "proj", min.cells = 0, min.features = 0, assay = "chrom_cnt"
  ) 
  
  
  #min_nCount_RNA <- 5000
  ## Dimensional reduction using kraken
  #sum(adata$nCount_chrom_cnt>min_nCount_RNA)
  #keep_cells <- adata$nCount_chrom_cnt > min_nCount_RNA
  #keep_cells <- adata$nCount_chrom_cnt > sort(adata$nCount_chrom_cnt, decreasing = TRUE)[10000]
  #sum(keep_cells)
  #adata <- adata[,keep_cells] ## Reduce to sensible number
  adata
  
  ########## Produce a count matrix on strain level
  #adata[["species_cnt"]] <- ChromToSpeciesCount(adata, mapSeq2strain)  #gives warning. coerce ourselves to dgCMatrix 
  mapSeq2strain_10strain <- mapSeq2strain[mapSeq2strain$strain!="phix" & mapSeq2strain$strain!="Xanthomonas campestris",]
  adata[["species_cnt"]] <- ChromToSpeciesCount(adata, mapSeq2strain_10strain)  #gives warning. coerce ourselves to dgCMatrix
  
  
  
  
  #Figure out which species has most reads in which cell
  cnt_spec <- GetAssayData(adata, assay = "species_cnt", layer = "counts")
  adata$species_aln <- rownames(cnt_spec)[apply(cnt_spec, 2, which.max)]
  adata$log_cnt <- log10(1+adata$nCount_chrom_cnt)
  adata$species_aln_short <- stringr::str_split_i(adata$species_aln,stringr::fixed(" ("),1)
  
  ##########
  ########## Generate UMAP
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
  
  
  ##########
  ########## Doublet detection using scDblFinder
  ##########
  
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts=GetAssayData(adata,layer = "counts"))) 
  
  sce <- scDblFinder::scDblFinder(sce)
  #14% doublets, in mock1
  adata$scDblFinder.score <- sce$scDblFinder.score  ### performs not so great
  adata$scDblFinder.class <- sce$scDblFinder.class
  
  if(FALSE){
    #does not appear to be reliable
    FeaturePlot(adata, "scDblFinder.score")
  }
  
  ### Cache results  
  saveRDS(adata, file.path(bascetRoot,"cache_adata_align.RDS"))
  adata <- readRDS(file.path(bascetRoot,"cache_adata_align.RDS"))
  
  ##########
  ########## Barnyard plot
  ##########
  
  #this assumes no phix or xantham left
  bp <- data.frame(
    maxc = MatrixGenerics::colMaxs(adata@assays$species_cnt$counts),
    totc = colSums(adata@assays$species_cnt$counts)
  )
  bp$restc <- bp$totc - bp$maxc
  mean(bp$maxc/bp$totc>0.8)
  adata$barnyard_doublet_score <- bp$maxc/bp$totc < 0.8
  adata$by_isdoublet <- bp$maxc/bp$totc < 0.95
  
  if(FALSE) {
    ggplot(bp, aes(maxc,restc)) + 
      geom_point() +
      theme_bw() +
      xlab("Dominant species count") +
      ylab("Other species count") #+ scale_x_log10() + scale_y_log10()
    ggsave(file.path(plotDir,"alignment_barnyard.pdf"), width = 5, height = 5)
  } else {
    #hist(log10(bp$maxc/bp$totc), breaks=100)
    
    ### This version highlights doublets better
    ggplot(bp, aes(maxc,restc)) + 
      geom_point() +
      theme_bw() +
      scale_x_log10()+
      scale_y_log10()+
      coord_cartesian(xlim=c(10,1e7),ylim=c(10,1e7))+
      geom_abline(slope=1, color="blue")+
      xlab("Dominant species count") +
      ylab("Other species count")
    #hist(bp$maxc/bp$totc)
    
  }
  
  
  if(FALSE) {
    #Use the ratio to remove doublets
    #adata$by_isdoublet <- bp$maxc/bp$totc < 0.95
    FeaturePlot(adata, "by_isdoublet")
    FeaturePlot(adata, "log_cnt")
  }
  
  
  ##########
  ########## Barnyard vs doublet prediction
  ##########
  
  FeaturePlot(adata, c("scDblFinder.score","barnyard_doublet_score"))
  
  ##########
  ########## Plot UMAPs
  ##########
  
  
  if(FALSE){
    
    #Deinococcus radiodurans  
    FeaturePlot(adata, features = mapSeq2strain[mapSeq2strain$strain=="Deinococcus radiodurans",]$id) 
    mapSeq2strain[mapSeq2strain$strain=="Deinococcus radiodurans",]
    ## possibly small plasmids?  ==> more simulated reads would help clustering. or higher connectivity
    
    #Bacillus
    FeaturePlot(adata, features = mapSeq2strain[mapSeq2strain$strain=="Bacillus pacificus",]$id)
    mean(adata@meta.data$species_aln_short=="Bacillus pacificus")
    
    FeaturePlot(adata[,!adata$by_isdoublet], features = mapSeq2strain[mapSeq2strain$strain=="Bacillus pacificus",]$id)
    
    
    ##### Bacillus contig sizes
    #cat /husky/henriksson/atrandi/bwa_ref/atcc/ATCC-AGP-24087822e79a496f/genome_ATCC_10987_687931d9b06b4cb4/assembly/Bacillus_pacificus_ATCC_10987.fasta.fai
    # ATCC_10987_contig_1	5235074	81	60	61
    # ATCC_10987_contig_2	147664	5322488	60	61
    # ATCC_10987_contig_3	33603	5472695	60	61
    # ATCC_10987_contig_4	19072	5506940	60	61
    # ATCC_10987_contig_5	3582	5526411	60	61
    # ATCC_10987_contig_6	2377	5530134	60	61
    # ATCC_10987_contig_7	1447	5532632	60	61
    
    FeaturePlot(adata, features = rownames(adata)[1:16]) #28 total
    FeaturePlot(adata, features = rownames(adata)[17:28]) #28 total
  }
  
  
  
  
  ### Plotting - max species per cell  --- barnyard filter!
  DimPlot(object = adata[,!adata$barnyard_doublet_score], label = TRUE, group.by = "species_aln_short") + 
    xlab("BWA1") + 
    ylab("BWA2") + ggtitle("") + NoLegend()
  #ggsave(file.path("~/github/scmetag-revision/plots/umap_mock_rna.png"), width = 7, height = 5)
  #ggsave(file.path("~/github/scmetag-revision/plots/umap_mock_atac.png"), width = 7, height = 5)
  ggsave(file.path(plotDir,"umap_alignment_species.pdf"), width = 7, height = 5)
  
  ### Plotting - max species per cell, doublets removed  [not as reliable]
  DimPlot(object = adata, label = TRUE, group.by = "species_aln_short") + 
    xlab("BWA1") + 
    ylab("BWA2") + 
    ggtitle("") 
  #ggsave(file.path(plotDir,"umap_alignment_nodoublet.pdf"), width = 7, height = 5)
  
  ### Plotting - max species per cell, doublets removed  [not as reliable]
  DimPlot(object = adata[,adata$scDblFinder.class=="singlet"], label = TRUE, group.by = "species_aln_short") + 
    xlab("BWA1") + 
    ylab("BWA2") + 
    ggtitle("") 
  ggsave(file.path(plotDir,"umap_alignment_nodoublet.pdf"), width = 7, height = 5)
  
  # Save distribution for simulation etc
  cnt <- adata@assays$species_cnt$counts
  df_align_bc_dist <- data.frame(
    species=adata$species_aln, 
    cnt=apply(cnt, 2, max)
  )
  saveRDS(df_align_bc_dist, file.path(bascetRoot,"align_metadata.RDS"))
  
  ##########
  ########## Kneeplot per species 
  ##########
  
  DefaultAssay(adata) <- "species_cnt"
  KneeplotPerSpecies(adata) #+ scale_y_log10(limits=c(1e3,1e7)) 
  ggsave(file.path(plotDir,"alignment_kneeplot.pdf"), width = 6, height = 4)
  
  ##########
  ########## Bias plot
  ##########
  
  df <- adata@meta.data
  df <- sqldf("select count(*) as cnt, species_aln_short from df group by species_aln_short")
  df$dataset <- dataset_name
  df
  saveRDS(df, file.path(bascetRoot,"count_per_species_aln.RDS"))
}







################################################################################
##################### Kraken first analysis ####################################
################################################################################

### Read matrix, rename 
mat <- ReadBascetCountMatrix(bascetRoot,"kraken_counts", verbose=FALSE)
dim(mat)
rownames(mat) <- stringr::str_remove(rownames(mat),stringr::fixed("BASCET_"))

adata <- CreateSeuratObject(mat, project = "proj", min.cells = 0, min.features = 0)

## Add KRAKEN consensus taxonomy to metadata
kraken_taxid <- KrakenFindConsensusTaxonomy(mat)
adata <- BascetAddMetaData(adata, kraken_taxid)

#### Save metadata for overlay in other sections (support merging different objects too)
saveRDS(kraken_taxid, file.path(bascetRoot,"kraken_metadata.RDS"))
saveRDS(adata@meta.data, file.path(bascetRoot,"adata_kraken_metadata.RDS"))

## Subset data
min_nCount_RNA <- 1000
sum(adata$nCount_RNA>min_nCount_RNA)
adata <- adata[,adata$nCount_RNA>min_nCount_RNA]
adata

#adata <- adata[,adata$species!="Homo sapiens"]  #clearly background!
adata$log_cnt <- log10(1+adata$nCount_RNA)
adata$perc_human <- adata@assays$RNA$counts["taxid-9607",]/colSums(adata@assays$RNA$counts) 

## Dimensional reduction using kraken
if(TRUE){
  #ATAC-seq style
  adata <- RunTFIDF(adata)
  adata <- FindTopFeatures(adata, min.cutoff = 'q0')
  adata <- RunSVD(adata)
  DepthCor(adata)
  adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 1:30, reduction.name = "kraken_umap")  ## depth seems to be less of a problem here
} else {
  #RNA-seq style
  adata <- NormalizeData(adata)
  adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
  adata <- ScaleData(adata, features = rownames(adata))
  adata <- RunPCA(adata, features = VariableFeatures(object = adata))
  adata <- RunUMAP(adata, dims = 1:20, reduction.name = "kraken_umap")
}

#saveRDS(adata, "/home/mahogny/github/rbiscvi/adata.RDS")



################################################################################
################## kraken - doublet analysis ###################################
################################################################################

sce <- SingleCellExperiment::SingleCellExperiment(list(counts=GetAssayData(adata,layer = "counts"))) 
sce <- scDblFinder::scDblFinder(sce)
adata$scDblFinder.score <- sce$scDblFinder.score
adata$scDblFinder.class <- sce$scDblFinder.class

################################################################################
################## kraken - UMAPs ##############################################
################################################################################

saveRDS(adata, file.path(bascetRoot,"cache_adata_kraken.RDS"))
adata <- readRDS(file.path(bascetRoot,"cache_adata_kraken.RDS"))
  
FeaturePlot(adata, features = c("scDblFinder.score")) + xlab("KRAKEN1") + ylab("KRAKEN2")
ggsave(file.path(plotDir, "umap_kraken_doublet_score.pdf"), width=7, height=7, limitsize=FALSE) 

DimPlot(adata, group.by = c("scDblFinder.class")) + xlab("KRAKEN1") + ylab("KRAKEN2")
ggsave(file.path(plotDir, "umap_kraken_doublet_class.pdf"), width=7, height=7, limitsize=FALSE) 

VlnPlot(adata, "log_cnt", group.by = "scDblFinder.class")
ggsave(file.path(plotDir, "kraken_doublet_violin.pdf"), width=3, height=7, limitsize=FALSE)


#### Save metadata for overlay in other sections (support merging different objects too)
saveRDS(kraken_taxid, file.path(bascetRoot,"kraken_metadata.RDS"))
saveRDS(adata@meta.data, file.path(bascetRoot,"adata_kraken_metadata.RDS"))  ### could subset all cells to these?


################################################################################
################## kraken - plots ##############################################
################################################################################

#####
##### UMAP plots
#####

## UMAP of phylum abundance according to KRAKEN
DimPlot(object = adata, label = TRUE, group.by = "phylum", reduction = "kraken_umap") + 
  NoLegend() + 
  xlab("KRAKEN1") + 
  ylab("KRAKEN2")
ggsave(file.path(plotDir, "umap_kraken_phylum.pdf"), width=7, height=7)

## UMAP of genus abundance according to KRAKEN
DimPlot(object = adata, label = TRUE, group.by = "genus", reduction = "kraken_umap") + 
  NoLegend() + 
  xlab("KRAKEN1") + 
  ylab("KRAKEN2")
ggsave(file.path(plotDir, "umap_kraken_genus.pdf"), width=7, height=7)

#For grants etc
DimPlot(object = adata, group.by = "genus", reduction = "kraken_umap") + 
  xlab("KRAKEN1") + 
  ylab("KRAKEN2")


## UMAP of species abundance according to KRAKEN
DimPlot(object = adata, label = TRUE, group.by = "species", reduction = "kraken_umap") + 
  NoLegend() + 
  xlab("KRAKEN1") + 
  ylab("KRAKEN2")
ggsave(file.path(plotDir, "umap_kraken_species.pdf"), width=7, height=7)


if(FALSE) {

  DimPlot(object = adata, label = TRUE, group.by = "species", reduction = "kraken_umap") + 
    NoLegend() + 
    xlab("KRAKEN1") + 
    ylab("KRAKEN2")
  
  adata$species[stringr::str_starts(adata$species,"Tan")]
  
}



#Compare with depth. "human" cells got few counts and end up in the middle
if(!stringr::str_detect(dataset_name,"mock")){
  FeaturePlot(adata, features = "taxid-339") ## "Xanthomonas campestris" 
  ggsave(file.path(plotDir, "umap_kraken_xanthomonas.pdf"), width=7, height=7, limitsize=FALSE) 
}

FeaturePlot(adata, features = c("log_cnt")) + xlab("KRAKEN1") + ylab("KRAKEN2")
ggsave(file.path(plotDir, "umap_kraken_logcnt.pdf"), width=7, height=7, limitsize=FALSE)

FeaturePlot(adata, features = c("perc_human")) + xlab("KRAKEN1") + ylab("KRAKEN2")        ## it really seems gone in saliva data!
ggsave(file.path(plotDir, "umap_kraken_human.pdf"), width=7, height=7, limitsize=FALSE)


#####
##### Other plots
#####

############# Barplot of species abundance according to KRAKEN
df <- adata@meta.data
df <- sqldf::sqldf("select species, count(*) as cnt from df group by species")
df <- df[order(df$cnt),]
df$species <- factor(df$species, levels = df$species)
ggplot(df[df$cnt>5,], aes(species, cnt)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  ylab("Cell count")
ggsave(file.path(plotDir, "kraken_hist_species.pdf"), width=7, height=7, limitsize=FALSE)

############# Barplot of phylum abundance according to KRAKEN
df <- adata@meta.data
df <- sqldf::sqldf("select phylum, count(*) as cnt from df group by phylum")
df <- df[order(df$cnt),]
df$phylum <- factor(df$phylum, levels = df$phylum)
ggplot(df[df$cnt>5,], aes(phylum, cnt)) + 
  geom_bar(stat="identity") +
  coord_flip() +
  theme_bw() +
  ylab("Cell count")
ggsave(file.path(plotDir, "kraken_hist_phylum.pdf"), width=7, height=7, limitsize=FALSE)

############# Barplot of genus abundance according to KRAKEN
df <- adata@meta.data
df <- sqldf::sqldf("select genus, count(*) as cnt from df group by genus")
df <- df[order(df$cnt),]
df$genus <- factor(df$genus, levels = df$genus)
ggplot(df[df$cnt>5,], aes(genus, cnt)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() +
  ylab("Cell count")
ggsave(file.path(plotDir, "kraken_hist_genus.pdf"), width=7, height=7, limitsize=FALSE) 

############### Kneeplots for kraken
showNumSpecies <- 10
if(stringr::str_detect(bascetRoot,"saliva")){
  showNumSpecies <- 30
}

KrakenKneePlot(adata, groupby = "phylum", showNumSpecies=showNumSpecies) 
ggsave(file.path(plotDir, "kraken_kneeplot_per_phylum.pdf"), width=7, height=7, limitsize=FALSE)

KrakenKneePlot(adata, groupby = "genus", showNumSpecies=showNumSpecies)
ggsave(file.path(plotDir, "kraken_kneeplot_per_genus.pdf"), width=7, height=7, limitsize=FALSE)

KrakenKneePlot(adata, groupby = "species", showNumSpecies=showNumSpecies)
ggsave(file.path(plotDir, "kraken_kneeplot_per_species.pdf"), width=7, height=7, limitsize=FALSE)





### Ground truth align data if mock
if(str_detect(bascetRoot,"mock")){
  
  adata_aligned <- readRDS(file.path(bascetRoot,"cache_adata_align.RDS"))
  #colnames(adata_aligned) <- stringr::str_remove(colnames(adata_aligned),stringr::fixed("BASCET_"))
  
  adata_with_aligned <- BascetAddMetaData(
    adata, 
    adata_aligned@meta.data[,"species_aln",drop=FALSE], 
    subsetCommon = TRUE) 
  
  DimPlot(object = adata_with_aligned, label = TRUE, group.by = "species_aln") + NoLegend() + xlab("KRAKEN1") + ylab("KRAKEN2")
  ggsave(file.path(plotDir,"umap_infokmer_vs_aln.pdf"), width = 7, height = 7)
  
}










################################################################################
################## informative KMER-analysis - distribution ####################
################################################################################

p_minh <- file.path(bascetRoot,"minhash_hist.csv")
p_use_kmers <- file.path(bascetRoot,"use_kmers.txt")
if(file.exists(p_minh)) {
  kmerHist <- BascetReadMinhashHistogram(bascetRoot)
  
  ### KMER count histogram
  kmerHist$rank <- 1:nrow(kmerHist)
  ggplot(kmerHist[sample(1:nrow(kmerHist),min(30000, nrow(kmerHist))),], aes(rank, cnt)) +   # not>2!!
    #ggplot(kmerHist[sample(which(kmerHist$cnt>10), 5000),], aes(rank, cnt)) +   # not>2!!
    geom_point() + 
    scale_x_log10() + 
    scale_y_log10() +
    theme_bw() +
    ylab("Count")+
    xlab("Rank")
  ggsave(file.path(plotDir,"info_kmerHist.pdf"), width = 5, height = 5)
  ggsave(file.path(plotDir,"info_kmerHist.png"), width = 5, height = 5)
}



################################################################################
################## informative KMER-analysis - umap etc ########################
################################################################################

if(file.exists(file.path(bascetRoot,paste0("kmer_counts.1.h5")))) {

  
  ### Read count matrix
  cnt <- ReadBascetCountMatrix(  
    bascetRoot,
    "kmer_counts"
  ) 
 
  cnt <- cnt[order(rowSums(cnt), decreasing = TRUE),]
  
  ## Load subset of real cells
  #keep_cells <- rowSums(cnt)>20000
  #sum(keep_cells)
  adata <- CreateSeuratObject(
    counts = cnt,
    assay = "infokmer"
  )
  
  
  ##### Add kraken metadata. keep cells based on kraken so we have annotation for all
  kraken_taxid <- readRDS(file.path(bascetRoot,"kraken_metadata.RDS"))
  adata <- BascetAddMetaData(adata, kraken_taxid, subsetCommon=TRUE)  
  
  ##### Perform ATAC-seq style analysis
  adata$log10_nCount_infokmer <- log10(adata$nCount_infokmer)
  adata <- RunTFIDF(adata)
  adata <- FindTopFeatures(adata, min.cutoff = 'q0')
  adata <- RunSVD(adata)
  
  DepthCor(adata)
  ggsave(file.path(plotDir,"infokmer_depthcor.pdf"), width = 4, height = 3)
  
  adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 1:30, reduction.name = "infokmers_umap")  ## dim 1 affected plenty if binarizing matrix
  
  saveRDS(adata, file.path(bascetRoot,"cache_adata_infokmer.RDS"))
  
  ##### 
  ##### Plots
  ##### 
  DimPlot(object = adata, label = TRUE, reduction = "infokmers_umap") + NoLegend() + xlab("IK1") + ylab("IK2")
  
  ### KRAKEN
  if(TRUE) {
    #Add kraken metadata
    #kraken_taxid <- readRDS(file.path(bascetRoot,"kraken_metadata.RDS"))
    #adata <- BascetAddMetaData(adata, kraken_taxid, subsetCommon = TRUE)

    DimPlot(object = adata, label = TRUE, group.by = "genus", reduction = "infokmers_umap") + NoLegend() + xlab("IK1") + ylab("IK2")
    ggsave(file.path(plotDir,"infokmer_umap_genus.pdf"), width = 5, height = 5)
    
    DimPlot(object = adata, label = TRUE, group.by = "phylum", reduction = "infokmers_umap") + NoLegend() + xlab("IK1") + ylab("IK2")
    ggsave(file.path(plotDir,"infokmer_umap_phylum.pdf"), width = 5, height = 5)

    DimPlot(
      object = adata[,adata$genus %in% names(which(table(adata$genus)>50))], 
      label = TRUE, group.by = "genus", reduction = "infokmers_umap")  + xlab("IK1") + ylab("IK2")
    #NoLegend()
    ggsave(file.path(plotDir,"infokmer_umap_genus_dominant.pdf"), width = 7, height = 7)

  }
  
  ### Ground truth align data if mock
  if(!str_detect(bascetRoot,"saliva")){
    adata_aligned <- readRDS(file.path(bascetRoot,"cache_adata_align.RDS"))
    colnames(adata_aligned) <- stringr::str_remove(colnames(adata_aligned),stringr::fixed("BASCET_"))
    
    adata_with_aligned <- BascetAddMetaData(
      adata, 
      adata_aligned@meta.data[,"species_aln",drop=FALSE], 
      subsetCommon = TRUE) 
    
    DimPlot(object = adata_with_aligned, label = TRUE, group.by = "species_aln") + NoLegend() + xlab("IK1") + ylab("IK2")
    ggsave(file.path(plotDir,"umap_infokmer_vs_aln.pdf"), width = 7, height = 7)
  }
  
}






################################################################################
################## Count sketch-based analysis #################################
################################################################################

#Load data as seurat object 
adata <- BascetLoadCountSketchMatrix(bascetRoot,inputName = "countsketch_mat.feather")
adata$log10_celldepth <- log10(adata$celldepth)

#Subset cells
#min_kmers <- 400000  #40k for saliva => 8k cells; option: pick top N cells for comparability
min_kmers <- 40000  #40k for saliva => 8k cells; option: pick top N cells for comparability
keep_cells <- adata$celldepth > min_kmers
sum(keep_cells)
adata <- adata[,keep_cells]

#if(FALSE) {
#  adata <- adata[,1:5000] ##################################### can we speed up if needed? #######################################
#}

#Add kraken metadata
kraken_taxid <- readRDS(file.path(bascetRoot,"kraken_metadata.RDS"))
adata <- BascetAddMetaData(adata, kraken_taxid, subsetCommon = TRUE) 

#Non-linear dimensional reduction
adata <- CountSketchUMAP(adata)

saveRDS(adata, file.path(bascetRoot,"cache_adata_cs.RDS"))

####################################################
##### Plotting UMAPS using kraken labels ###########
####################################################

DimPlot(object = adata, label = TRUE) + NoLegend() + xlab("CS1") + ylab("CS2")

if(TRUE) {
  DimPlot(object = adata, label = TRUE, group.by = "phylum") + NoLegend() + xlab("CS1") + ylab("CS2")
  ggsave(file.path(plotDir,"umap_countsketch_phylum.pdf"), width = 7, height = 7)
  
  DimPlot(object = adata[,adata$genus!="Streptococcus"], label = TRUE, group.by = "phylum") + NoLegend() + xlab("CS1") + ylab("CS2")
  ggsave(file.path(plotDir,"umap_countsketch_phylum_nostrepto.pdf"), width = 7, height = 7)
  
  DimPlot(object = adata, label = TRUE, group.by = "genus") + NoLegend() + xlab("CS1") + ylab("CS2")
  ggsave(file.path(plotDir,"umap_countsketch_genus.pdf"), width = 7, height = 7)
  
  DimPlot(object = adata, label = TRUE, group.by = "species") + NoLegend() + xlab("CS1") + ylab("CS2")
  ggsave(file.path(plotDir,"umap_countsketch_species.pdf"), width = 7, height = 7)
  
  FeaturePlot(adata, features = "log10_celldepth") + xlab("CS1") + ylab("CS2")
  ggsave(file.path(plotDir,"umap_countsketch_logkmer.pdf"), width = 7, height = 7)
  
}


##########################################################
##### Plotting UMAPS using ground truth labels ###########
##########################################################

if(str_detect(bascetRoot,"mock")){
  
  adata_aligned <- readRDS(file.path(bascetRoot,"cache_adata_align.RDS"))

  adata_with_aligned <- BascetAddMetaData(
    adata, 
    adata_aligned@meta.data[,"species_aln",drop=FALSE], 
    subsetCommon = TRUE) 
  
  DimPlot(object = adata_with_aligned, label = TRUE, group.by = "species_aln") + NoLegend() + xlab("CS1") + ylab("CS2")
  ggsave(file.path(plotDir,"umap_countsketch_vs_aln.pdf"), width = 7, height = 7)
  
}

