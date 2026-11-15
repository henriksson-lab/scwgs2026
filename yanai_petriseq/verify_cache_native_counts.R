# Compare every selected cached gene count with independently decoded native H5.
VerifyYanaiCacheNativeCounts <- function(root,object) {
  counts <- SeuratObject::GetAssayData(object,assay="RNA",layer="counts")
  metadata <- object@meta.data
  stopifnot(identical(colnames(counts),rownames(metadata)),"library" %in% names(metadata))
  mapping <- read.csv(file.path(root,"yanai_gene_count_libraries.csv"))
  stopifnot(setequal(unique(metadata$library),mapping$library_prefix))
  records <- vector("list",nrow(mapping))
  for(k in seq_len(nrow(mapping))) {
    path <- file.path(root,mapping$counts[k])
    genes <- gsub("_","-",as.character(rhdf5::h5read(path,"/var/_index")),fixed=TRUE)
    cells <- as.character(rhdf5::h5read(path,"/obs/_index"))
    stopifnot(!anyDuplicated(genes),!anyDuplicated(cells),all(rownames(counts) %in% genes))
    native <- methods::new("dgCMatrix",Dim=as.integer(c(length(genes),length(cells))),
      Dimnames=list(genes,cells),p=as.integer(rhdf5::h5read(path,"/X/indptr")),
      i=as.integer(rhdf5::h5read(path,"/X/indices")),
      x=as.numeric(rhdf5::h5read(path,"/X/data")))
    methods::validObject(native)
    selected <- colnames(counts)[metadata$library==mapping$library_prefix[k]]
    stopifnot(length(selected)>0L,all(selected %in% cells))
    native <- native[,selected,drop=FALSE]
    omitted <- setdiff(genes,rownames(counts))
    if(length(omitted)) stopifnot(Matrix::nnzero(native[omitted,,drop=FALSE])==0L)
    expected <- native[rownames(counts),,drop=FALSE]
    actual <- counts[,selected,drop=FALSE]
    delta <- actual-expected
    differences <- Matrix::nnzero(delta)
    if(differences>0L) {
      tab <- Matrix::summary(delta)
      tab <- head(tab[tab$x!=0,,drop=FALSE],1000L)
      write.csv(data.frame(gene=rownames(actual)[tab$i],cell_id=colnames(actual)[tab$j],
        cache_minus_native=tab$x),file.path(root,"umap","cache_native_count_mismatch.csv"),row.names=FALSE)
      stop("Cached gene counts differ from native H5 for ",mapping$library_prefix[k],
           ": ",differences," entries")
    }
    records[[k]] <- data.frame(library=mapping$library_prefix[k],cells=length(selected),
      native_genes=length(genes),cached_genes=nrow(counts),
      omitted_native_genes_zero=TRUE,exact_every_cell_gene_count=TRUE,
      selected_gene_read_assignments=sum(actual))
    message("Cached counts match every native gene entry: ",mapping$library_prefix[k])
  }
  table <- do.call(rbind,records)
  stopifnot(sum(table$cells)==ncol(counts))
  write.table(table,file.path(root,"umap","cache_native_counts.validation.tsv"),
              row.names=FALSE,quote=FALSE,sep="\t")
  report <- list(status="verified",scope="every cached cell/gene entry against native H5, including zero genes omitted from Seurat",
    libraries=nrow(mapping),cells=ncol(counts),cached_genes=nrow(counts),
    exact_every_cell_gene_count=TRUE,omitted_native_genes_zero=TRUE,
    checked_at=format(Sys.time(),tz="UTC",usetz=TRUE))
  jsonlite::write_json(report,file.path(root,"umap","cache_native_counts.validation.json"),
                       auto_unbox=TRUE,pretty=TRUE)
  invisible(report)
}
