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

dataset_name <- "zheng2022"

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


#Not sure how to interpret
zheng_meta <- read.csv("/husky/henriksson/atrandi/zheng2022/tofetch.runinfo.csv")
nrow(zheng_meta) #20,915 cells


################################################################################
##################### Processing time ##########################################
################################################################################







################################################################################
##################### Kraken first analysis ####################################
################################################################################

### Read matrix, rename 
mat <- ReadBascetCountMatrix(bascetRoot,"kraken_counts", verbose=FALSE)
#dim(mat)
#rownames(mat) <- stringr::str_remove(rownames(mat),stringr::fixed("BASCET_"))

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
  #ATAC-seq style. think not the best way here
  adata <- RunTFIDF(adata)
  adata <- FindTopFeatures(adata, min.cutoff = 'q0')
  adata <- RunSVD(adata)
  DepthCor(adata)
  adata <- RunUMAP(object = adata, reduction = 'lsi', dims = 1:30, reduction.name = "kraken_umap")  ## depth seems to be less of a problem here
  
} else {
  
  adata <- NormalizeData(adata)
  adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
  adata <- ScaleData(adata, features = rownames(adata))
  adata <- RunPCA(adata, features = VariableFeatures(object = adata))
  adata <- RunUMAP(adata, dims = 1:20, reduction.name = "kraken_umap")
}



################################################################################
################## kraken - doublet analysis ###################################
################################################################################

## Using scDblFinder.
## need to use bioconductors classes
sce <- SingleCellExperiment::SingleCellExperiment(list(counts=adata@assays$RNA$counts)) #cols should be cell names, rows gene names
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
ggsave(file.path(plotDir, "umap_kraken_phylum.png"), width=7, height=7)

## UMAP of genus abundance according to KRAKEN
DimPlot(object = adata, label = TRUE, group.by = "genus", reduction = "kraken_umap") + 
  NoLegend() + 
  xlab("KRAKEN1") + 
  ylab("KRAKEN2")
ggsave(file.path(plotDir, "umap_kraken_genus.pdf"), width=14, height=14)
ggsave(file.path(plotDir, "umap_kraken_genus.png"), width=14, height=14)

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



#Compare with depth. "human" cells got few counts and end up in the middle
if(stringr::str_detect(dataset_name,"saliva")){
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


################################################################################
################## informative KMER-analysis - distribution ####################
################################################################################

p_minh <- file.path(bascetRoot,"minhash_hist.csv")
p_use_kmers <- file.path(bascetRoot,"use_kmers.txt")
if(file.exists(p_minh)) {
  kmerHist <- BascetReadMinhashHistogram(bascetRoot)
  if(!file.exists(p_use_kmers)) {
    print(dim(kmerHist))    
    kmerHist$rand_index <- sample(1:nrow(kmerHist)) #add a tie breaker
    kmerHist <- kmerHist[order(kmerHist$cnt, kmerHist$rand_index, decreasing=TRUE),]

    #Pick the most common KMERs    
    picked_kmers <- kmerHist$kmer[1:100000]

#     picked_kmers <- ChooseInformativeKMERs( ### 0.002 => 11212 this killed conda!
#       kmerHist,
# #      minFreq = 0.2 #simulated2; much higher cutoff needed for 34k 
# #      minFreq = 0.002 #novaseq1, really few!
# #      minFreq = 0.0015 #saliva, really few!
#       minFreq = 0.001 #novaseq2, really few!
# #      minFreq = 0.0005 #novaseq3, really few!
#       ) ################################################### option of max count?
    writeLines(picked_kmers, p_use_kmers)
    bascetRoot
  }
  
  
  
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


### Read count matrix
cnt <- ReadBascetCountMatrix(  
  bascetRoot,
  "kmer_counts"
) 
dim(cnt)

cnt <- cnt[order(rowSums(cnt), decreasing = TRUE),]

## Load subset of real cells
#keep_cells <- rowSums(cnt)>20000
#sum(keep_cells)
adata <- CreateSeuratObject(
  counts = cnt,#[1:10000,],  #cutoff on cell count 
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
#DimPlot(object = adata, label = TRUE,reduction = "infokmers_umap") + NoLegend() + xlab("IK1") + ylab("IK2")

saveRDS(adata, file.path(bascetRoot,"cache_adata_infokmer.RDS"))

##### 
##### Plots
##### 

DimPlot(object = adata, label = TRUE, group.by = "genus", reduction = "infokmers_umap") + NoLegend() + xlab("IK1") + ylab("IK2")
ggsave(file.path(plotDir,"infokmer_umap_genus.pdf"), width = 5, height = 5)

DimPlot(object = adata, label = TRUE, group.by = "phylum", reduction = "infokmers_umap") + NoLegend() + xlab("IK1") + ylab("IK2")
ggsave(file.path(plotDir,"infokmer_umap_phylum.pdf"), width = 5, height = 5)

#FeaturePlot(adata, "log10_nCount_infokmer", reduction = "infokmers_umap") + xlab("IK1") + ylab("IK2")
#ggsave(file.path(plotDir,"infokmer_umap_kmercount.pdf"), width = 5, height = 5)


#DimPlot(object = adata[,adata$genus=="Bacillus"], label = TRUE, group.by = "genus", reduction = "infokmers_umap") + NoLegend()
#DimPlot(object = adata, label = TRUE, reduction = "infokmers_umap") + NoLegend()
#table(adata$genus) #MOCK: cereibacter dominates plenty!!

DimPlot(
  object = adata[,adata$genus %in% names(which(table(adata$genus)>50))], 
  label = TRUE, group.by = "genus", reduction = "infokmers_umap")  + xlab("IK1") + ylab("IK2")
  #NoLegend()
ggsave(file.path(plotDir,"infokmer_umap_genus_dominant.pdf"), width = 7, height = 7)



################################################################################
################## Count sketch-based analysis #################################
################################################################################

#Load data as seurat object 
adata <- BascetLoadCountSketchMatrix(bascetRoot,inputName = "countsketch_mat.csv")
#adata <- BascetLoadCountSketchMatrix(bascetRoot,inputName = "countsketch.0.tsv") ############# update!
adata$log10_celldepth <- log10(adata$celldepth)

#Subset cells
#min_kmers <- 400000  #40k for saliva => 8k cells; option: pick top N cells for comparability
min_kmers <- 40000  #40k for saliva => 8k cells; option: pick top N cells for comparability
keep_cells <- adata$celldepth > min_kmers
sum(keep_cells)
adata <- adata[,keep_cells]
if(FALSE) {
  adata <- adata[,1:5000] ##################################### can we speed up if needed? #######################################
}

#Add kraken metadata
kraken_taxid <- readRDS(file.path(bascetRoot,"kraken_metadata.RDS"))
adata <- BascetAddMetaData(adata, kraken_taxid, subsetCommon = TRUE) 

#Non-linear dimensional reduction
adata <- CountSketchUMAP(adata, nn.method = "hnsw")

CountSketchUMAP <- function(adata, reduction="kmersketch", outReduction="kmersketch_umap", metric="cosine", n_neighbors = 30, doFast=FALSE, seed = 42) {
  dims <- 1:ncol(adata@reductions[[reduction]]@cell.embeddings)

  n_sgd_threads <- 0
  if(doFast) {
    print("Fast UMAP mode. Note that this means the UMAP will be different between runs (non-deterministic SGD updates)")
    n_sgd_threads <- parallel::detectCores()
  }
  
  emb <- adata@reductions[[reduction]]@cell.embeddings
  um  <- uwot::umap(
    emb, 
    metric = "cosine", 
    n_neighbors = n_neighbors,
    nn_method = "hnsw",
    n_threads = parallel::detectCores(),
    n_sgd_threads = n_sgd_threads, 
    seed = 42
  )
  colnames(um) <- c("CS_1", "CS_2")
  adata[[outReduction]] <- CreateDimReducObject(
    embeddings = um, key = "CS_", assay = DefaultAssay(adata)
  )  
  adata
}
adata <- CountSketchUMAP(adata, doFast=FALSE)



### CountSketchUMAP




saveRDS(adata, file.path(bascetRoot,"cache_adata_cs.RDS"))

####################################################
##### Plotting UMAPS using kraken labels ###########
####################################################

DimPlot(object = adata, label = TRUE) + NoLegend() + xlab("CS1") + ylab("CS2")

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


##########################################################
##### Plotting UMAPS using ground truth labels ###########
##########################################################


adata_aligned <- readRDS(file.path(bascetRoot,"cache_adata_align.RDS"))
colnames(adata_aligned) <- stringr::str_remove(colnames(adata_aligned),stringr::fixed("BASCET_"))

adata_with_aligned <- BascetAddMetaData(
  adata, 
  adata_aligned@meta.data[,"species_aln",drop=FALSE], 
  subsetCommon = TRUE) 

DimPlot(object = adata_with_aligned, label = TRUE, group.by = "species_aln") + NoLegend() + xlab("CS1") + ylab("CS2")
ggsave(file.path(plotDir,"umap_countsketch_vs_aln.pdf"), width = 7, height = 7)


