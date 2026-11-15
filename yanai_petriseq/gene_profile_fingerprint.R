# Provenance for the optional variable-gene warning investigation.
YanaiGeneProfileFingerprint <- function(root) {
  mapping <- read.csv(file.path(root,"yanai_gene_count_libraries.csv"))
  paths <- normalizePath(file.path(root,mapping$counts),mustWork=TRUE)
  info <- file.info(paths)
  stopifnot(all(is.finite(info$size)))
  digest::digest(list(
    feature_selection="zorn_sorted_raw_gene_ids/vst2000/v2",
    gene_sort_locale=Sys.getlocale("LC_COLLATE"),
    selection=digest::digest(file=file.path(root,"umap","gene_cell_selection.csv"),
                             algo="sha256",serialize=FALSE),
    counts=data.frame(file=paths,size=info$size,mtime=as.numeric(info$mtime)),
    Seurat_version=as.character(utils::packageVersion("Seurat"))),algo="sha256")
}
