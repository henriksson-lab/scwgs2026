
bascetRoot <- "/husky/henriksson/atrandi//v9_wgs_mock5_pta"
adata <- readRDS(file.path(bascetRoot,"cache_adata_kraken.RDS"))



p1 <- FeaturePlot(adata, "scDblFinder.score")
p2 <- DimPlot(adata, group.by = "scDblFinder.class")

plotDirAll <- file.path("~/github/scmetag-revision/plots/all")


ptot <- egg::ggarrange(p1,p2,nrow=1)
ptot
ggsave(plot = ptot, filename = file.path(plotDirAll,"supnote_doublet.png"), width = 10, height = 5)

file.path(plotDirAll,"supnote_doublet.png")
