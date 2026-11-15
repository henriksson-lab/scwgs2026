#!/usr/bin/env Rscript
# Replot the existing Ling Kraken2 embedding; no analysis is recomputed.
# Optional arguments: dataset directory, output directory.
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)>=1L) args[1] else "/husky/henriksson/atrandi/v9_ling"
plotDir <- if(length(args)>=2L) args[2] else file.path(root,"plots")
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})
adata <- readRDS(file.path(root,"cache_adata_kraken.RDS"))
groups <- c("library","genus","sample","sequencing")
stopifnot(all(groups %in% colnames(adata[[]])),
          "kraken_umap" %in% Reductions(adata))
titles <- c("Library","Genus","Sample","Sequencing")
panels <- lapply(seq_along(groups),function(i) {
  group <- groups[i]
  p <- DimPlot(adata,reduction="kraken_umap",group.by=group,
               pt.size=1,raster=TRUE,raster.dpi=c(1000,1000),
               shuffle=TRUE,seed=42,combine=FALSE)[[1]]
  if(group=="genus") {
    counts <- table(adata$genus)
    genera <- names(counts)[counts>=50 &
      !names(counts) %in% c("","NA","Unclassified","unclassified")]
    if(length(genera)) p <- LabelClusters(p,id="genus",clusters=genera,
                                         repel=TRUE,size=3,seed=42)
    p <- p + NoLegend()
  }
  p + labs(title=titles[i],x="KRAKEN1",y="KRAKEN2") +
    theme(axis.text=element_blank(),axis.ticks=element_blank(),
          plot.title=element_text(face="bold"))
})
p <- wrap_plots(panels,ncol=2,nrow=2) + plot_annotation(tag_levels="A")
dir.create(plotDir,recursive=TRUE,showWarnings=FALSE)
for(ext in c("pdf","png")) {
  destination <- file.path(plotDir,paste0("umap_kraken_library_genus_sample_sequencing_2x2.",ext))
  ggsave(destination,plot=p,width=14,height=11,units="in",dpi=300,
         bg="white",limitsize=FALSE)
  message("Saved ",destination)
}
