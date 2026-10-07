library(Zorn)

#list_color_mock <- c(khroma::color("muted")(n = 9),"#DDDDDD")

library(Seurat)
library(Signac)
library(sqldf)
library(stringr)
library(ggplot2)
library(Matrix)
library(GenomicRanges)
library(future)
library(ggplot2)
library(ggalluvial)
library(ggrastr)


plotDirAll <- file.path("~/github/scmetag-revision/plots/all")

bascetRoot <- "/husky/henriksson/atrandi/v9_wgs_saliva2_mda"




################################################################################
##################### Load Quast data ##########################################
################################################################################


adata <- readRDS(file.path(bascetRoot,"cache_adata_kraken.RDS"))
quast_obs <- BascetReadObsFromH5(
  bascetRoot,
  "quast"
)
adata <- BascetAddMetaData(adata, quast_obs)





################################################################################
##################### Quast-based analysis #####################################
################################################################################

if(FALSE) {
  
  ggplot(adata@meta.data, aes(total_length, largest_contig, color=log10(log_cnt))) +  
    rasterise(geom_point(), dpi = 300, dev = "ragg_png") +
    scale_x_log10() + 
    scale_y_log10() +
    theme_bw() +
    xlab("Total contig length") + 
    ylab("Largest contig")

  ggplot(adata@meta.data, aes(total_length, contigs, color=log10(log_cnt))) + 
    rasterise(geom_point(), dpi = 300, dev = "ragg_png") +
    scale_x_log10() + 
    scale_y_log10() +
    theme_bw() +
    xlab("Total contig length") + 
    ylab("Number of contigs")

  ggplot(adata@meta.data, aes(n50, contigs, color=log10(log_cnt))) + 
    rasterise(geom_point(), dpi = 300, dev = "ragg_png") +
    scale_x_log10() + 
    scale_y_log10() +
    theme_bw() +
    xlab("N50") +  #N50 ~ weighted median statistic such that 50% of the entire assembly is contained in contigs or scaffolds equal to or larger than this value.
    ylab("Number of contigs")

}



ggplot(adata@meta.data[adata@meta.data$scDblFinder.score<0.68,], aes(n50, total_length, color=log10(log_cnt))) + 
  rasterise(geom_point(), dpi = 300, dev = "ragg_png") +
  scale_x_log10() + 
  scale_y_log10() +
  theme_bw() +
  xlab("N50") +
  ylab("Total contig length")
ggsave( file.path(plotDirAll,"saliva_quastQC.pdf"), width = 5, height = 5)

file.path(plotDirAll,"saliva_quastQC.pdf")

