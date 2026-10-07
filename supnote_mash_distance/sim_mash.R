library(stringr)
library(Biostrings)


################################################################################
###### Single species, 100 cells, log-scale coverage 0.01x .. 10x ##############
################################################################################


insert_max_size <- 800
insert_min_size <- 400
num_cells <- 100
readlen <- 150
qread <- str_flatten(rep("F", readlen))

cov_min <- 0.01
cov_max <- 10
coverages <- 10 ** seq(log10(cov_min), log10(cov_max), length.out = num_cells)

outdir <- "./sim_mash"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

fa_rootdir <- "/husky/henriksson/atrandi/bwamem_ref/simulated_ref10/separate"
all_fa_files <- list.files(fa_rootdir)
list_genomes <- unique(str_split_i(all_fa_files, " ", 1))

#### Pick the first genome only
cur_genome <- list_genomes[1]
list_fa <- all_fa_files[stringr::str_starts(all_fa_files, cur_genome)]
print(list_fa)
allfa <- Biostrings::readDNAStringSet(file.path(fa_rootdir, list_fa))
allfa_str <- as.character(allfa)
allfa_len <- str_length(allfa)
genome_total_len <- sum(allfa_len)

cat("Using genome:", cur_genome, "\n")
cat("Total genome length:", genome_total_len, "\n")


#### For each cell, simulate reads at its target coverage
for (cell_i in 1:num_cells) {

  cur_cov <- coverages[cell_i]
  cell_name <- sprintf("cell%05d_cov%.4fx", cell_i, cur_cov)

  #### Paired-end coverage: num_pairs * 2 * readlen / genome_len = coverage
  num_reads <- max(1, round(cur_cov * genome_total_len / (2 * readlen)))

  cat(sprintf("cell %d/%d  cov=%.4fx  pairs=%d\n", cell_i, num_cells, cur_cov, num_reads))

  #### Sample chromosome/plasmid index for each read, weighted by length
  use_chrom <- sample(1:length(allfa_str), prob = allfa_len, replace = TRUE, size = num_reads)

  r1_list <- character(0)
  r2_list <- character(0)

  for (cur_chrom in unique(use_chrom)) {
    get_num_reads <- sum(use_chrom == cur_chrom)

    ch <- allfa_str[cur_chrom]
    totlen <- str_length(ch)
    start_pos <- round(runif(get_num_reads, min = 1, max = totlen - insert_max_size))
    end_pos <- start_pos + round(runif(get_num_reads, min = insert_min_size, max = insert_max_size))

    all_inserts <- str_sub(ch, start_pos, end_pos)

    insert_r1 <- str_sub(all_inserts, 1, readlen)
    insert_r2_rc <- str_sub(all_inserts, str_length(all_inserts) - (readlen - 1), str_length(all_inserts))
    insert_r2 <- as.character(reverseComplement(DNAStringSet(insert_r2_rc)))

    r1_list <- c(r1_list, insert_r1)
    r2_list <- c(r2_list, insert_r2)
  }

  #### Write FASTQ files for this cell
  read_ids <- sprintf("@%s_read%07d", cell_name, seq_along(r1_list))

  fq1_lines <- as.vector(rbind(paste0(read_ids, "/1"), r1_list, "+", qread))
  fq2_lines <- as.vector(rbind(paste0(read_ids, "/2"), r2_list, "+", qread))

  fq1_name <- file.path(outdir, paste0(cell_name, "_R1.fastq.gz"))
  fq2_name <- file.path(outdir, paste0(cell_name, "_R2.fastq.gz"))

  con1 <- gzfile(fq1_name, "w")
  writeLines(fq1_lines, con1)
  close(con1)

  con2 <- gzfile(fq2_name, "w")
  writeLines(fq2_lines, con2)
  close(con2)
}

cat("done\n")
