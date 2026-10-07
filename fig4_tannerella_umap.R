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
##################### Compare samples with Kraken - load #######################
################################################################################


adata <- readRDS(file.path(bascetRoot,"cache_adata_kraken.RDS"))


#### 
tdata <- readRDS("/husky/henriksson/atrandi/v9_wgs_saliva2_mda/seurat_tannerella.RDS")
df_tan <- data.frame(
  row.names = colnames(tdata),
  tan_log_cnt = tdata$log_cnt
)
adata <- BascetAddMetaData(adata, df_tan)
adata$tan_log_cnt[is.na(adata$tan_log_cnt)] <- 0



################################################################################
##################### Tanerella counts on UMAPs ################################
################################################################################


tdata <- readRDS("/husky/henriksson/atrandi/v9_wgs_saliva2_mda/seurat_tannerella.RDS")
df_tan <- data.frame(
  row.names = colnames(tdata),
  tan_log_cnt = tdata$log_cnt
)


adata1 <- readRDS(file.path(bascetRoot,"cache_adata_kraken.RDS"))
adata2 <- readRDS(file.path(bascetRoot,"cache_adata_infokmer.RDS"))
adata3 <- readRDS(file.path(bascetRoot,"cache_adata_cs.RDS"))



plot_tan <- function(adata) {
  adata <- BascetAddMetaData(adata, df_tan)
  adata$tan_log_cnt[is.na(adata$tan_log_cnt)] <- 0
  FeaturePlot(adata, "tan_log_cnt", raster=TRUE, raster.dpi = c(200,200))
}

p1 <- plot_tan(adata1)
p2 <- plot_tan(adata2)
p3 <- plot_tan(adata3)

ptot <- egg::ggarrange(plots=list(p1,p2,p3),ncol=1)
ptot
ggsave(plot = ptot, file.path(plotDirAll,"tan_umap.pdf"), width=3, height=8)

