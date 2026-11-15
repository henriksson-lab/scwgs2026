# Independent whole-BAM oracle using R/Bioconductor interval overlaps.
# Compares EVERY cell/gene entry with the native Bascet read-count matrix.
# Intended for the isolated single-library preflight, not UMAP processing.
VerifyBamGeneCounts <- function(root,yieldSize=250000L) {
  stopifnot(requireNamespace("Matrix",quietly=TRUE))
  mapping <- read.csv(file.path(root,"yanai_gene_count_libraries.csv"))
  stopifnot(nrow(mapping)==1L)
  path <- file.path(root,mapping$counts)
  cells <- as.character(rhdf5::h5read(path,"/obs/_index"))
  features <- as.character(rhdf5::h5read(path,"/var/_index"))
  ptr <- as.integer(rhdf5::h5read(path,"/X/indptr"))
  values <- as.numeric(rhdf5::h5read(path,"/X/data"))
  indices <- as.integer(rhdf5::h5read(path,"/X/indices"))
  # CSR (cells x genes) becomes CSC (genes x cells) by transposing slots.
  actual <- methods::new("dgCMatrix",Dim=as.integer(c(length(features),length(cells))),
    Dimnames=list(features,cells),p=ptr,i=indices,x=values)
  methods::validObject(actual)
  oracle <- Matrix::sparseMatrix(i=integer(),j=integer(),x=numeric(),
    dims=dim(actual),dimnames=dimnames(actual))
  gff <- read.delim(mapping$annotation,header=FALSE,comment.char="#",quote="",colClasses="character")
  gff <- gff[gff[[3]]=="gene",]
  ids <- sub(".*(?:^|;)locus_tag=([^;]+).*","\\1",gff[[9]],perl=TRUE)
  geneIndex <- match(ids,features)
  stopifnot(!anyNA(geneIndex),!anyDuplicated(geneIndex))
  genes <- GenomicRanges::GRanges(seqnames=gff[[1]],
    ranges=IRanges::IRanges(start=as.integer(gff[[4]]),end=as.integer(gff[[5]])),strand="*")
  params <- Rsamtools::ScanBamParam(what=c("rname","pos","cigar"),tag="CB",
    flag=Rsamtools::scanBamFlag(isUnmappedQuery=FALSE,isSecondaryAlignment=FALSE,
                              isSupplementaryAlignment=FALSE))
  bam <- Rsamtools::BamFile(file.path(root,mapping$alignment),yieldSize=yieldSize)
  open(bam);on.exit(close(bam))
  processed <- 0
  repeat {
    chunk <- Rsamtools::scanBam(bam,param=params)[[1]]
    n <- length(chunk$pos)
    if(!n) {
      # A filtered yield can be empty before EOF; later primary alignments
      # must still be checked (and trailing unmapped blocks still read).
      if(Rsamtools::isIncomplete(bam)) next
      break
    }
    barcodes <- as.character(chunk$tag$CB)
    stopifnot(length(barcodes)==n,!anyNA(barcodes))
    cellIndex <- match(barcodes,cells)
    stopifnot(!anyNA(cellIndex))
    # Parse each distinct CIGAR once per chunk. M,D,N,=,X consume reference.
    # Insertions and soft/hard clips do not extend the reference interval.
    types <- unique(chunk$cigar)
    spans <- vapply(regmatches(types,gregexpr("[0-9]+[MDN=X]",types,perl=TRUE)),
      function(parts)sum(as.integer(sub(".$","",parts))),integer(1))
    widths <- spans[match(chunk$cigar,types)]
    stopifnot(all(widths>0L),all(chunk$pos>=1L))
    reads <- GenomicRanges::GRanges(seqnames=as.character(chunk$rname),
      ranges=IRanges::IRanges(start=chunk$pos,width=widths),strand="*")
    hits <- GenomicRanges::findOverlaps(reads,genes,ignore.strand=TRUE)
    block <- Matrix::sparseMatrix(i=geneIndex[S4Vectors::subjectHits(hits)],
      j=cellIndex[S4Vectors::queryHits(hits)],x=rep(1,length(hits)),dims=dim(actual))
    oracle <- oracle+block
    processed <- processed+n
    if(processed%%1000000L==0L) message("Independent BAM oracle: ",processed," primary mapped records")
  }
  difference <- oracle-actual
  if(Matrix::nnzero(difference)!=0L) {
    stop("Independent whole-BAM gene counts disagree with Bascet: ",
         Matrix::nnzero(difference)," cell/gene entries differ")
  }
  report <- list(status="verified",scope="every primary mapped BAM record and every cell/gene entry",
    library=mapping$library_prefix,cells=length(cells),genes=length(features),
    primary_mapped_records=processed,gene_read_assignments=sum(oracle),
    nonzero_entries=Matrix::nnzero(oracle),exact_all_cell_gene_counts=TRUE,
    includes_mapq_zero_and_duplicates=TRUE,umi_deduplication=FALSE,
    unstranded=TRUE,full_count_for_each_overlapping_gene=TRUE,
    counts_sha256=digest::digest(file=path,algo="sha256",serialize=FALSE),
    checked_at=format(Sys.time(),tz="UTC",usetz=TRUE))
  jsonlite::write_json(report,file.path(root,"bam_gene_counts.validation.json"),auto_unbox=TRUE,pretty=TRUE)
  message("Independent whole-BAM oracle matches every Bascet cell/gene count")
  invisible(report)
}
