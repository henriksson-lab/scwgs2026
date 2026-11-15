# Gene read counts after Bascet BWA-MEM2 alignment. Source from prep_yanai.R.
# This helper also permits rerunning the counting stage on existing alignments.
RunYanaiGeneCounts <- function(bascetRoot, debstat, annotation,
                              runner, bascetInstance) {
  if(!all(c("countMode", "cellFiles") %in% names(formals(BascetCountFeature)))) {
    stop("The local Zorn checkout needs countMode/cellFiles support for gene read counting")
  }
  # Bascet's historical countfeature default deduplicates UMIs. Explicit reads
  # mode counts primary mapped alignment records, with no UMI deduplication.
  # It is unstranded, includes MAPQ=0 primary alignments, excludes unmapped,
  # secondary and supplementary records, and gives a full count to EVERY
  # overlapping gene. Thus overlapping genes can share a read's count.
  # This applies to interval overlap (alignment reference span), not exon union.
  groups <- length(debstat$group_prefix)
  stopifnot(groups>0L, length(debstat$list_picked_cells)==groups)
  bamNames <- detectShardsForFile(bascetRoot, "aligned_pos")
  expected <- paste0("aligned_pos.",seq_len(groups),".bam")
  if(!identical(bamNames, expected)) {
    stop("Expected one aligned_pos.N.bam per independent library in group order")
  }
  stopifnot(all(file.exists(file.path(bascetRoot,paste0(bamNames,".bai")))))

  # A barcode with no gene-overlapping reads still belongs in the expression
  # matrix. Seed every library's rows from the retained debarcoding barcodes.
  cellFiles <- file.path(bascetRoot,paste0("gene_cells.",seq_len(groups),".txt"))
  for(i in seq_len(groups)) {
    cells <- debstat$list_picked_cells[[i]]
    stopifnot(length(cells)>0L, !anyDuplicated(cells))
    writeLines(cells, cellFiles[i])
  }
  job <- BascetCountFeature(
    bascetRoot, inputName="aligned_pos", outputName="genecount",
    gffFile=annotation, useFeature="gene",
    attrGeneId="locus_tag", attrGeneName="Name",
    countMode="reads", cellFiles=cellFiles,
    runner=runner, bascetInstance=bascetInstance
  )
  WaitForJob(job)
  mapping <- data.frame(
    group=seq_len(groups), library_prefix=debstat$group_prefix,
    alignment=bamNames, counts=paste0("genecount.",seq_len(groups),".h5"),
    retained_cells=lengths(debstat$list_picked_cells),
    count_unit="primary_mapped_read_alignment", strand="unstranded",
    overlap_assignment="full_count_to_each_overlapping_gene",
    gene_id_attribute="locus_tag", annotation=annotation
  )
  stopifnot(all(file.exists(file.path(bascetRoot,mapping$counts))))
  write.csv(mapping,file.path(bascetRoot,"yanai_gene_count_libraries.csv"),row.names=FALSE)
  message("Gene read counts ready for ",groups," independent libraries; zero-count cells retained")
  invisible(mapping)
}
