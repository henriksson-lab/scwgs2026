#!/usr/bin/env Rscript
# Describe the full UMAP input without changing its cell or feature selection.
root <- Sys.getenv("YANAI_ROOT","/husky/henriksson/atrandi/v9_rnaseq_yanai")
suppressPackageStartupMessages(library(data.table))
tab <- fread(file.path(root,"umap","gene_cell_selection.csv"),
             select=c("library","gene_read_counts","genes_detected","retained"))
valid <- nrow(tab)
tab <- tab[retained==TRUE]
stopifnot(valid==31876328L,nrow(tab)==864570L,
          all(tab$gene_read_counts>1000),all(tab$genes_detected>=1))
manifest <- read.delim(file.path(root,"fastq_manifest.tsv"))
tab[,species:=manifest$scientific_name[match(library,manifest$library_prefix)]]
stopifnot(!anyNA(tab$species),uniqueN(tab$library)==88L)
summarize <- function(x) list(cells=nrow(x),
  gene_read_count_quantiles=as.list(quantile(x$gene_read_counts,c(0,.01,.05,.1,.25,.5,.75,.9,.95,1))),
  detected_gene_quantiles=as.list(quantile(x$genes_detected,c(0,.01,.05,.1,.25,.5,.75,.9,.95,1))),
  cells_at_most_1_gene=sum(x$genes_detected<=1),
  cells_at_most_2_genes=sum(x$genes_detected<=2),
  cells_at_most_5_genes=sum(x$genes_detected<=5),
  cells_at_most_10_genes=sum(x$genes_detected<=10),
  cells_at_most_50_genes=sum(x$genes_detected<=50))
report <- list(scope="every cell selected by the normal >1000 gene-read threshold",
  input_barcodes=valid,libraries=uniqueN(tab$library),all=summarize(tab),
  species=lapply(split(tab,tab$species),summarize),
  selection_modified=FALSE,checked_at=format(Sys.time(),tz="UTC",usetz=TRUE))
byLibrary <- tab[,.(cells=.N,median_genes=as.numeric(median(genes_detected)),
  cells_at_most_1_gene=sum(genes_detected<=1),cells_at_most_5_genes=sum(genes_detected<=5),
  cells_at_most_10_genes=sum(genes_detected<=10)),by=.(library,species)]
fwrite(byLibrary,file.path(root,"umap","gene_input_complexity.tsv"),sep="\t")
jsonlite::write_json(report,file.path(root,"umap","gene_input_complexity.json"),
                     auto_unbox=TRUE,pretty=TRUE)
print(report$all)
