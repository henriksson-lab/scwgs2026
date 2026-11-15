#!/usr/bin/env Rscript
# Supplementary Information Figure 12: aligned-count UMAPs, colored by species.
# Uses post_unified.R caches; only the RNA-style reduction needs computing.
# Usage: Rscript supfig12_alignment_workflows.R [dataset directory] [output directory]
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

atac <- readRDS(file.path(root, "cache_adata_align.RDS"))
alignment <- atac[[]]
stopifnot(DefaultAssay(atac) == "chrom_cnt", "chrom_umap" %in% Reductions(atac),
          "lsi" %in% Reductions(atac), "species_aln_short" %in% colnames(alignment))

# Use the same cached aligned chromosome counts and cells for both workflows.
species <- alignment$species_aln_short
keep <- !is.na(species) & nzchar(species)
stopifnot(any(keep))
message("Plotting ", sum(keep), " of ", ncol(atac), " aligned cells with species labels.")
species_levels <- sort(unique(species[keep]))
colors <- setNames(c(khroma::color("muted")(n = 9), "#DDDDDD"),
                   sort(unique(na.omit(alignment$species_aln_short))))

# Cache just the new embedding locally; delete this file to recompute it.
rna_cache <- file.path(plot_dir, "supfig12_alignment_rna_umap.RDS")
if (file.exists(rna_cache)) {
  rna_umap <- readRDS(rna_cache)
  stopifnot(identical(rownames(rna_umap), colnames(atac)))
} else {
  # Start from raw counts: the cached data layer contains TF-IDF values.
  rna <- CreateSeuratObject(GetAssayData(atac, assay = "chrom_cnt", layer = "counts"),
                            assay = "chrom_cnt", min.cells = 0, min.features = 0)
  # Alignment RNA-seq branch of post_unified.R: all features, up to 50 PCs.
  set.seed(42)
  rna <- NormalizeData(rna)
  rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = nrow(rna))
  rna <- ScaleData(rna, features = rownames(rna))
  rna <- RunPCA(rna, features = VariableFeatures(rna),
                npcs = min(50, nrow(rna) - 1), seed.use = 42)
  rna <- RunUMAP(rna, dims = 1:min(50, nrow(rna) - 2),
                 reduction.name = "chrom_umap", seed.use = 42)
  rna_umap <- Embeddings(rna, "chrom_umap")
  saveRDS(rna_umap, rna_cache)
  rm(rna)
}

panel <- function(embedding, title) {
  df <- data.frame(UMAP1 = embedding[keep, 1], UMAP2 = embedding[keep, 2],
                   species = factor(species[keep], levels = species_levels))
  p <- ggplot(df, aes(UMAP1, UMAP2, color = species)) +
    geom_point(size = 0.3) +
    scale_color_manual(values = colors, drop = FALSE, name = "Species (alignment)") +
    coord_fixed() +
    labs(title = title, x = "BWA1", y = "BWA2") +
    theme_classic(base_size = 11) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          legend.position = "none")
  LabelClusters(p, id = "species", repel = TRUE, box = FALSE,
                size = 2.7, color = "black", seed = 42,
                max.overlaps = Inf, max.time = 2)
}
p <- (panel(rna_umap, "RNA-seq workflow\n(log normalization / PCA)") +
      panel(Embeddings(atac, "chrom_umap"), "ATAC-seq workflow\n(TF-IDF / LSI)")) +
  plot_layout(ncol = 2) +
  plot_annotation(title = "Supplementary Information Figure 12", tag_levels = "A")
for (ext in c("pdf", "png")) {
  path <- file.path(plot_dir, paste0("supplementary_figure12_alignment_workflows.", ext))
  ggsave(path, plot = p, width = 13, height = 5, dpi = 300, bg = "white")
  message("Saved ", path)
}

# Feature weights for the first ten cached SVD/LSI components (chromosomes).
svd_weights <- VizDimLoadings(atac, reduction = "lsi", dims = 1:10,
                             nfeatures = 15, balanced = TRUE, ncol = 5) +
  plot_annotation(title = "Aligned-count SVD feature weights (TF-IDF / LSI)")
# Signac's companion diagnostic, also used in post_unified.R.
svd_depth <- DepthCor(atac, reduction = "lsi", n = 10) +
  labs(title = "Aligned-count SVD components: correlation with read depth",
       subtitle = "Assay: chrom_cnt; reduction: LSI")
for (ext in c("pdf", "png")) {
  ggsave(file.path(plot_dir, paste0("supplementary_figure12_alignment_svd_weights.", ext)),
         plot = svd_weights, width = 20, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(plot_dir, paste0("supplementary_figure12_alignment_svd_depthcor.", ext)),
         plot = svd_depth, width = 7, height = 4, dpi = 300, bg = "white")
}
