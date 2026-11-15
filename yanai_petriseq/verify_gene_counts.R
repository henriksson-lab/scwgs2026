#!/usr/bin/env Rscript
# Audit every native gene matrix against all observed valid barcodes and GFF genes.
VerifyYanaiGeneCounts <- function(root,selection=NULL) {
  mapping <- read.csv(file.path(root,"yanai_gene_count_libraries.csv"))
  sharding <- readRDS(file.path(root,"yanai_sharding.rds"))
  stopifnot(!anyDuplicated(mapping$group),!anyDuplicated(mapping$library_prefix),
            setequal(mapping$library_prefix,sharding$group_prefix),
            all(mapping$count_unit=="primary_mapped_read_alignment"),
            all(mapping$strand=="unstranded"),
            all(mapping$overlap_assignment=="full_count_to_each_overlapping_gene"),
            all(mapping$gene_id_attribute=="locus_tag"))
  annotation <- unique(mapping$annotation)
  stopifnot(length(annotation)==1L)
  reference <- jsonlite::read_json(file.path(dirname(annotation),"reference.json"),
                                   simplifyVector=TRUE)
  samtools <- Sys.getenv("SAMTOOLS",Sys.which("samtools"))
  if(!nzchar(samtools)) samtools <- path.expand("~/miniconda3/bin/samtools")
  stopifnot(file.exists(samtools))
  gff <- read.delim(annotation,header=FALSE,comment.char="#",quote="",
                    colClasses="character")
  attributes <- gff[gff[[3]]=="gene",9]
  stopifnot(all(grepl("(^|;)locus_tag=",attributes)))
  genes <- sub(".*(?:^|;)locus_tag=([^;]+).*","\\1",attributes,perl=TRUE)
  stopifnot(length(genes)==11800L,!anyDuplicated(genes))
  report <- lapply(seq_len(nrow(mapping)),function(i) {
    bam <- file.path(root,mapping$alignment[i])
    stopifnot(file.exists(bam),file.exists(paste0(bam,".bai")))
    quick <- system2(samtools,c("quickcheck","-v",shQuote(bam)),stdout=TRUE,stderr=TRUE)
    if(!is.null(attr(quick,"status")) && attr(quick,"status")!=0L) {
      stop("BAM integrity check failed: ",bam,"\n",paste(quick,collapse="\n"))
    }
    stats <- system2(samtools,c("idxstats",shQuote(bam)),stdout=TRUE,stderr=TRUE)
    if(!is.null(attr(stats,"status")) && attr(stats,"status")!=0L) {
      stop("BAM index check failed: ",bam,"\n",paste(stats,collapse="\n"))
    }
    stats <- read.table(text=stats,sep="\t",col.names=c("id","length","mapped","unmapped"))
    contigs <- stats[stats$id!="*",]
    index <- match(reference$contigs$id,contigs$id)
    stopifnot(nrow(contigs)==nrow(reference$contigs),!anyDuplicated(contigs$id),
              !anyNA(index),all(contigs$length[index]==reference$contigs$length))
    group <- match(mapping$library_prefix[i],sharding$group_prefix)
    expected <- sharding$list_group_bc_count[[group]]$cellid
    picked <- sharding$list_picked_cells[[group]]
    stopifnot(!anyDuplicated(expected),setequal(expected,picked),
              length(expected)==mapping$retained_cells[i])
    seed <- readLines(file.path(root,paste0("gene_cells.",mapping$group[i],".txt")))
    stopifnot(!anyDuplicated(seed),setequal(seed,expected))
    path <- file.path(root,mapping$counts[i])
    ids <- as.character(rhdf5::h5read(path,"/obs/_index"))
    features <- as.character(rhdf5::h5read(path,"/var/_index"))
    attrs <- rhdf5::h5readAttributes(path,"/X")
    ptr <- as.numeric(rhdf5::h5read(path,"/X/indptr"))
    values <- as.numeric(rhdf5::h5read(path,"/X/data"))
    indices <- as.integer(rhdf5::h5read(path,"/X/indices"))
    stopifnot(!anyDuplicated(ids),setequal(ids,expected),
              all(startsWith(ids,paste0(mapping$library_prefix[i],"_"))),
              !anyDuplicated(features),setequal(features,genes),
              identical(attrs[["encoding-type"]],"csr_matrix"),
              all(attrs$shape==c(length(ids),length(genes))),
              length(ptr)==length(ids)+1L,ptr[1]==0,all(diff(ptr)>=0),
              tail(ptr,1)==length(values),length(indices)==length(values),
              all(indices>=0L & indices<length(genes)),
              all(is.finite(values)),all(values>=0),all(values==floor(values)))
    cumulative <- c(0,cumsum(values))
    counts <- cumulative[ptr[-1]+1]-cumulative[ptr[-length(ptr)]+1]
    if(!is.null(selection)) {
      rows <- which(selection$library==mapping$library_prefix[i])
      index <- match(ids,selection$cell_id[rows])
      stopifnot(length(rows)==length(ids),!anyNA(index),
                all(counts==selection$gene_read_counts[rows][index]),
                all(diff(ptr)==selection$genes_detected[rows][index]))
    }
    data.frame(library=mapping$library_prefix[i],barcodes=length(ids),
               genes=length(features),gene_read_counts=sum(counts),
               zero_gene_count_barcodes=sum(counts==0),
               exact_valid_barcode_coverage=TRUE,exact_reference_genes=TRUE,
               bam_integrity_valid=TRUE,bam_reference_contigs_valid=TRUE,
               sparse_counts_valid=TRUE,
               umap_input_counts_checked=!is.null(selection))
  })
  report <- do.call(rbind,report)
  write.table(report,file.path(root,"gene_input_validation.tsv"),
              sep="\t",quote=FALSE,row.names=FALSE)
  report
}
if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  root <- if(length(args)) args[1] else
    Sys.getenv("YANAI_ROOT","/husky/henriksson/atrandi/v9_rnaseq_yanai")
  selection <- read.csv(file.path(root,"umap","gene_cell_selection.csv"))
  print(VerifyYanaiGeneCounts(root,selection),row.names=FALSE)
}
