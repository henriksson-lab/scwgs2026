### Assemble quast output data for supplemental
library(Zorn)

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


################################################################################
##################### Load Quast data ##########################################
################################################################################

quast_out_dir <-file.path(plotDirAll, "sup_quast")
dir.create(quast_out_dir)


list_bascetRoot <- c(
  saliva="/husky/henriksson/atrandi/v9_wgs_saliva2_mda",
  MDA="/husky/henriksson/atrandi/v9_wgs_mock4_mda",
  PTA="/husky/henriksson/atrandi/v9_wgs_mock5_pta"
)

for(n in names(list_bascetRoot)) {
  print(n)
  bascetRoot <- list_bascetRoot[[n]]
  
  adata <- readRDS(file.path(bascetRoot,"cache_adata_kraken.RDS"))
  quast_obs <- BascetReadObsFromH5(
    bascetRoot,
    "quast"
  )
  adata <- BascetAddMetaData(adata, quast_obs)

  tow <- adata@meta.data
  tow <- tow[,-c(1:6)]
  
  write.csv(
    file=file.path(quast_out_dir, paste0(n,".csv")), 
    quote = FALSE,
    tow
  )  

}
