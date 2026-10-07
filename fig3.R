library(Zorn)


list_color_mock <- c(khroma::color("muted")(n = 9),"#DDDDDD")
plotDirAll <- file.path("~/github/scmetag-revision/plots/all")

################################################################################
########################## Fig 3 umap matrix ###################################
################################################################################

datasets_dir <- "/husky/henriksson/atrandi"

######### 
######### simulated mock data
######### 

bdata_dataset_name <- "v9_simulated5"
bdata1 <- readRDS(file.path(datasets_dir,bdata_dataset_name,"cache_adata_align.RDS"))
bdata2 <- readRDS(file.path(datasets_dir,bdata_dataset_name,"cache_adata_kraken.RDS"))
bdata3 <- readRDS(file.path(datasets_dir,bdata_dataset_name,"cache_adata_infokmer.RDS"))
bdata4 <- readRDS(file.path(datasets_dir,bdata_dataset_name,"cache_adata_cs.RDS"))

#take cell names from bdata1
bdata2$species_aln_short <- bdata1@meta.data[colnames(bdata2),]$species_aln_short
bdata3$species_aln_short <- bdata1@meta.data[colnames(bdata3),]$species_aln_short
bdata4$species_aln_short <- bdata1@meta.data[colnames(bdata4),]$species_aln_short

######### 
######### real mock data
######### 

adata_dataset_name <- "v9_wgs_mock5_pta"
adata1 <- readRDS(file.path(datasets_dir,adata_dataset_name,"cache_adata_align.RDS"))
adata2 <- readRDS(file.path(datasets_dir,adata_dataset_name,"cache_adata_kraken.RDS"))
adata3 <- readRDS(file.path(datasets_dir,adata_dataset_name,"cache_adata_infokmer.RDS"))
adata4 <- readRDS(file.path(datasets_dir,adata_dataset_name,"cache_adata_cs.RDS"))

#take cell names from adata1
adata2$species_aln_short <- adata1@meta.data[colnames(adata2),]$species_aln_short
adata3$species_aln_short <- adata1@meta.data[colnames(adata3),]$species_aln_short
adata4$species_aln_short <- adata1@meta.data[colnames(adata4),]$species_aln_short



do_label <- FALSE
a_p1 <- DimPlot(object = adata1, label = do_label, group.by = "species_aln_short") + 
  xlab("BWA1") + ylab("BWA2") + 
  ggtitle("") + 
  scale_color_manual(values = list_color_mock)

a_p2 <- DimPlot(object = adata2[,!is.na(adata2$species_aln_short)], label = do_label, group.by = "species_aln_short") + 
  xlab("KRAKEN1") + ylab("KRAKEN2") + 
  ggtitle("") + 
  scale_color_manual(values = list_color_mock)

a_p3 <- DimPlot(object = adata3[,!is.na(adata3$species_aln_short)], label = do_label, group.by = "species_aln_short") + #one cluster is a bit meah
  xlab("INFOKMER1") + ylab("INFOKMER2") + 
  ggtitle("") + 
  scale_color_manual(values = list_color_mock)

a_p4 <- DimPlot(object = adata4[,!is.na(adata4$species_aln_short)], label = do_label, group.by = "species_aln_short") + 
  xlab("CS1") + ylab("CS2") + 
  ggtitle("") + 
  scale_color_manual(values = list_color_mock)





b_p1 <- DimPlot(object = bdata1, label = do_label, group.by = "species_aln_short") + 
  xlab("BWA1") + ylab("BWA2") + 
  ggtitle("") + 
  scale_color_manual(values = list_color_mock)

b_p2 <- DimPlot(object = bdata2[,!is.na(bdata2$species_aln_short)], label = do_label, group.by = "species_aln_short") + 
  xlab("KRAKEN1") + ylab("KRAKEN2") + 
  ggtitle("") + 
  scale_color_manual(values = list_color_mock)

b_p3 <- DimPlot(object = bdata3[,!is.na(bdata3$species_aln_short)], label = do_label, group.by = "species_aln_short") + #one cluster is a bit meah
  xlab("INFOKMER1") + ylab("INFOKMER2") + 
  ggtitle("") + 
  scale_color_manual(values = list_color_mock)

b_p4 <- DimPlot(object = bdata4[,!is.na(bdata4$species_aln_short)], label = do_label, group.by = "species_aln_short") + 
  xlab("CS1") + ylab("CS2") + 
  ggtitle("") + 
  scale_color_manual(values = list_color_mock)


library(patchwork)

plot_list <- list(
  b_p1+theme(legend.position = "none")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+xlab("")+ylab(""),
  b_p2+theme(legend.position = "none")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+xlab("")+ylab(""),
  b_p3+theme(legend.position = "none")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+xlab("")+ylab(""),
  b_p4+theme(legend.position = "none")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+xlab("")+ylab(""),
  
  a_p1+theme(legend.position = "none")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+xlab("")+ylab(""),
  a_p2+theme(legend.position = "none")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+xlab("")+ylab(""),
  a_p3+theme(legend.position = "none")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+xlab("")+ylab(""),
  a_p4+theme(legend.position = "none")+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+xlab("")+ylab("")
)
ptot <- wrap_plots(plot_list, ncol = 4) &
  theme(
    plot.margin = margin(0, 0, 0, 0),
    axis.title  = element_blank(),
    axis.text   = element_blank(),
    axis.ticks  = element_blank()
  )
ptot
ggsave(plot = ptot, file.path(plotDirAll, "fig3_umaps.png"), width = 16*0.5, height = 8*0.5)




FeaturePlot(object = adata1, features = "nCount_species_cnt") + 
  xlab("BWA1") + ylab("BWA2") + 
  ggtitle("") 



################################################################################
########################## Fig 3 JL graph ######################################
################################################################################

PlotJohnsonLindenstraussMinDim(list_eps = c(0.1,0.2,0.3))
ggsave(file.path(plotDirAll, "fig3_JL.svg"), width = 4, height = 2)





################################################################################
########################## Fig 3 kmer hisogram #################################
################################################################################

bascetRoot <- file.path("/husky/henriksson/atrandi",adata_dataset_name)


p_minh <- file.path(bascetRoot,"minhash_hist.csv")
p_use_kmers <- file.path(bascetRoot,"use_kmers.txt")

kmerHist <- BascetReadMinhashHistogram(bascetRoot)

### KMER count histogram
kmerHist$rank <- 1:nrow(kmerHist)
ggplot(kmerHist[sample(1:nrow(kmerHist),min(30000, nrow(kmerHist))),], aes(rank, cnt)) +   
  rasterise(geom_point(), dpi = 300, dev = "ragg_png") +
  scale_x_log10() + 
  scale_y_log10() +
  theme_bw() +
  ylab("Count")+
  xlab("Rank")
ggsave(file.path(plotDirAll,"fig3_info_kmerHist.svg"), width = 3, height = 3)
ggsave(file.path(plotDirAll,"fig3_info_kmerHist.pdf"), width = 3, height = 3)



################################################################################
########################## Fig 3 kmer on umap ##################################
################################################################################

adata <- adata3 

most_common_kmer <- tail(rownames(adata)[order(rowSums(adata@assays$infokmer))], n=1)
most_common_kmer

### Most common KMER
FeaturePlot(adata, features = most_common_kmer) +   
  xlab("IK1") + ylab("IK2")
ggsave(file.path(plotDirAll,"fig3_info_kmer_ex.png"), width = 4, height = 4)


