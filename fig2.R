library(Zorn)


list_color_mock <- c(khroma::color("muted")(n = 9),"#DDDDDD")
list_color_mock_total <- c(list_color_mock, "#000000")

library(Seurat)
library(Signac)
library(sqldf)
library(stringr)
library(ggplot2)
library(Matrix)
library(GenomicRanges)
library(future)



plotDirAll <- file.path("~/github/scmetag-revision/plots/all")
if(!file.exists(plotDirAll)){
  dir.create(plotDirAll, recursive = TRUE)
  print("made plot dir all")
}


################################################################################
################# Fig 2 kneeplots for different lysis ########################## per species
################################################################################


adata_mda <- readRDS("/husky/henriksson/atrandi/v9_wgs_mock4_mda/cache_adata_align.RDS")
DefaultAssay(adata_mda) <- "species_cnt"

adata_pta <- readRDS("/husky/henriksson/atrandi/v9_wgs_mock5_pta/cache_adata_align.RDS")
DefaultAssay(adata_pta) <- "species_cnt"



## option of showTotal=TRUE
p1 <- KneeplotPerSpecies(adata_mda, showTotal = TRUE) + 
  scale_color_manual(values = list_color_mock_total) + 
  scale_y_log10(limits=c(1e3,1e7)) 
p2 <- KneeplotPerSpecies(adata_pta, showTotal = TRUE) + 
  scale_color_manual(values = list_color_mock_total) + 
  scale_y_log10(limits=c(1e3,1e7)) 

ptot <- egg::ggarrange(
  p1,p2,
  nrow=1)
ptot

ggsave(plot = ptot, file.path(plotDirAll, "fig2_strain_kneeplots.svg"), width = 9, height = 2)
#ggsave(plot = ptot, file.path(plotDirAll, "fig2_strain_kneeplots.svg"), width = 12, height = 4)



################################################################################
################# Fig 2 kneeplots for different lysis ########################## all species as a whole
################################################################################

#no longer used
if(FALSE) {
  plot_one_kneeplot_all <- function(adata){
    df <- data.frame(
      cnt=sort(colSums(adata@assays$species_cnt$counts), decreasing = TRUE)
    )
    df$index <- 1:nrow(df)
    ggplot(df, aes(index, cnt)) + geom_line() + scale_x_log10() + scale_y_log10() + theme_bw() + ylab("Read count") + xlab("Cell index")
  }
  
  p1 <- plot_one_kneeplot_all(adata_mda) + scale_y_log10(limits=c(1e3,1e7)) 
  p2 <- plot_one_kneeplot_all(adata_pta) + scale_y_log10(limits=c(1e3,1e7)) 
  
  ptot <- egg::ggarrange(
    p1,p2,
    nrow=1)
  ptot
  
  ggsave(plot = ptot, file.path(plotDirAll, "fig2_allread_kneeplots.svg"), width = 12, height = 4)
  
}



################################################################################
################# Fig 2 count per species ###################################### ################# can we define cutoffs differently?
################################################################################

#10k reads make sense?



plot_one_speciesfrac <- function(adata_mda){
  df <- data.frame(species=adata_mda@meta.data$species_aln[adata_mda$nCount_species_cnt>10000])
  df <- sqldf::sqldf("select species, count(*) as cnt from df group by species")
  df$frac <- df$cnt / sum(df$cnt) * 100
  df$species <- stringr::str_split_i(df$species," \\(",1)

  df$gram <- "pos"
  df$gram[stringr::str_detect(df$species,"coli")] <- "neg"
  df$gram[stringr::str_detect(df$species,"sphaer")] <- "neg"
  
  df$species <- factor(df$species, levels=rev(unique(df$species)))

  ggplot(df, aes(species, frac, fill=gram)) + 
    geom_bar(stat="identity") +  
    coord_flip() + 
    scale_fill_manual(values=c("#bb5566", "#004488")) +
    ylab("Fraction %") +
    xlab("Species") +
    theme_bw() 
  
}


ptot <- egg::ggarrange(
  plot_one_speciesfrac(adata_mda),
  plot_one_speciesfrac(adata_pta),
  nrow=1)
ptot
ggsave(plot = ptot, file.path(plotDirAll, "fig2_allread_species_count_barplot.svg"), width = 10, height = 2.5)





##########
########## Purity plot
##########



listp <- list()
for(i in 1:2) {
  if(i==1) {
    maxy <- 1100
    adata <- adata_mda
    adata_name <- "MDA"
  } else {
    maxy <- 3700
    adata <- adata_pta
    adata_name <- "PTA"
  }

  #this assumes no phix or xantham left
  bp <- data.frame(
    maxc = MatrixGenerics::colMaxs(adata@assays$species_cnt$counts),
    totc = colSums(adata@assays$species_cnt$counts)
  )
  bp$restc <- bp$totc - bp$maxc
  mean(bp$maxc/bp$totc>0.8)
  adata$barnyard_doublet_score <- bp$maxc/bp$totc < 0.8
  adata$by_isdoublet <- bp$maxc/bp$totc < 0.95
  
  bp$purity <- bp$maxc/bp$totc*100
  onep <- ggplot(bp, aes(x=purity)) + 
    geom_histogram(color = "black", fill = "steelblue") + 
    xlim(0,100) + 
    ylim(0,maxy) +
    theme_bw() + 
    ylab("Cell count") + 
    xlab("Purity %") + ggtitle(adata_name)
  onep
  
  print(mean(bp$purity>0.9*100))
  
  listp[[adata_name]] <- onep
}
ptot <- egg::ggarrange(plots=listp,nrow=1)
ggsave(plot=ptot, file.path(plotDirAll, "fig2_purity.svg"), width = 6, height = 2)


