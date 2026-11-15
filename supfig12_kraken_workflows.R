#!/usr/bin/env Rscript
# Supplementary Information Figure 12: Kraken UMAPs, colored by alignment species.
# Uses post_unified.R caches; only the RNA-style reduction needs computing.
# Usage: Rscript supfig12_kraken_workflows.R [dataset directory] [output directory]
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(patchwork)
})
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[1] else "/husky/henriksson/atrandi/v9_wgs_mock5_pta"
plot_dir <- if (length(args) >= 2L) args[2] else file.path("plots", "v9_wgs_mock5_pta")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

atac <- readRDS(file.path(root, "cache_adata_kraken.RDS"))
alignment <- readRDS(file.path(root, "cache_adata_align.RDS"))[[]]
stopifnot(DefaultAssay(atac) == "RNA", "kraken_umap" %in% Reductions(atac),
          "lsi" %in% Reductions(atac), "species_aln_short" %in% colnames(alignment))

# Keep the cached Kraken population (>1000 counts) for both reductions.
# Match alignment species by barcode only for plotting.
barcode <- function(x) sub("^BASCET_", "", x)
stopifnot(!anyDuplicated(barcode(rownames(alignment))))
species <- alignment$species_aln_short[
  match(barcode(colnames(atac)), barcode(rownames(alignment)))
]
keep <- !is.na(species) & nzchar(species)
stopifnot(any(keep))
message("Plotting ", sum(keep), " of ", ncol(atac), " Kraken cells with alignment species labels.")
colors <- setNames(c(khroma::color("muted")(n = 9), "#DDDDDD"),
                   sort(unique(na.omit(alignment$species_aln_short))))

# Cache just the new embedding locally; delete this file to recompute it.
rna_cache <- file.path(plot_dir, "supfig12_kraken_rna_umap.RDS")
if (file.exists(rna_cache)) {
  rna_umap <- readRDS(rna_cache)
  stopifnot(identical(rownames(rna_umap), colnames(atac)))
} else {
  # Start from raw counts: the cached data layer contains TF-IDF values.
  rna <- CreateSeuratObject(GetAssayData(atac, assay = "RNA", layer = "counts"),
                            assay = "RNA", min.cells = 0, min.features = 0)
  # Kraken RNA-seq branch of post_unified.R: 2000 variable features, 20 PCs.
  set.seed(42)
  rna <- NormalizeData(rna)
  rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
  rna <- ScaleData(rna, features = rownames(rna))
  rna <- RunPCA(rna, features = VariableFeatures(rna), seed.use = 42)
  rna <- RunUMAP(rna, dims = 1:20,
                 reduction.name = "kraken_umap", seed.use = 42)
  rna_umap <- Embeddings(rna, "kraken_umap")
  saveRDS(rna_umap, rna_cache)
  rm(rna)
}

# Standard Seurat UMAPs, sharing the same alignment-defined species colors.
atac$species_aln_short <- factor(species, levels = names(colors))
atac[["rna_umap"]] <- CreateDimReducObject(rna_umap, key = "RNAUMAP_", assay = "RNA")
plot_data <- atac[, keep]
p_rna <- DimPlot(plot_data, reduction = "rna_umap", group.by = "species_aln_short",
                 cols = colors, label = TRUE, repel = TRUE,
                 raster = TRUE, raster.dpi = c(1000, 1000), pt.size = 3) +
  NoLegend() + labs(title = "RNA-seq", x = "KRAKEN1", y = "KRAKEN2")
p_atac <- DimPlot(plot_data, reduction = "kraken_umap", group.by = "species_aln_short",
                  cols = colors, label = TRUE, repel = TRUE,
                  raster = TRUE, raster.dpi = c(1000, 1000), pt.size = 3) +
  NoLegend() + labs(title = "ATAC-seq", x = "KRAKEN1", y = "KRAKEN2")
p <- p_rna + p_atac
for (ext in c("pdf", "png")) {
  path <- file.path(plot_dir, paste0("supplementary_figure12_kraken_workflows.", ext))
  ggsave(path, plot = p, width = 12, height = 4.5, dpi = 300, bg = "white")
  message("Saved ", path)
}

# Feature weights for the first ten cached SVD/LSI components (Kraken taxids).
svd_weights <- VizDimLoadings(atac, reduction = "lsi", dims = 1:10,
                             nfeatures = 15, balanced = TRUE, ncol = 5) +
  plot_annotation(title = "Kraken SVD feature weights (TF-IDF / LSI)")
# Signac's companion diagnostic, also used in post_unified.R.
svd_depth <- DepthCor(atac, reduction = "lsi", n = 10) +
  labs(title = "Kraken SVD components: correlation with read depth",
       subtitle = "Assay: RNA; reduction: LSI")
for (ext in c("pdf", "png")) {
  ggsave(file.path(plot_dir, paste0("supplementary_figure12_kraken_svd_weights.", ext)),
         plot = svd_weights, width = 20, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(plot_dir, paste0("supplementary_figure12_kraken_svd_depthcor.", ext)),
         plot = svd_depth, width = 7, height = 4, dpi = 300, bg = "white")
}
