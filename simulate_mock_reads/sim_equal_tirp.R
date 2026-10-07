library(stringr)
library(Biostrings)


################################################################################
###### Make a mock community, with plasmids, roughly uniform distribution ######
################################################################################


insert_max_size <- 800
num_cells <- 50
num_shards <- 20
readlen <- 150
qread <- str_flatten(rep("F", readlen))

fa_rootdir <- "/husky/henriksson/atrandi/bwamem_ref/simulated_ref10/separate"
all_fa_files <- list.files(fa_rootdir)
list_genomes <- unique(str_split_i(all_fa_files, " ", 1))

#### Pre-load all fasta sequences, grouped by genome.
#### Each genome may have multiple fasta files (chromosomes, plasmids); a single
#### cell should draw reads uniformly across all bases of all fasta files that
#### belong to its genome.
genome_data <- list()
for (cur_genome in list_genomes) {
  list_fa <- all_fa_files[stringr::str_starts(all_fa_files, cur_genome)]
  print(list_fa)
  allfa <- Biostrings::readDNAStringSet(file.path(fa_rootdir, list_fa))
  genome_data[[cur_genome]] <- list(
    str = as.character(allfa),
    len = str_length(allfa)
  )
}


pos_index <- 1
for (cur_shard in 1:num_shards) {
  allreads_df <- list()
  reads_i <- 1
  print(paste("shard", cur_shard))

  for (cur_genome in list_genomes) {

    allfa_str <- genome_data[[cur_genome]]$str
    allfa_len <- genome_data[[cur_genome]]$len

    #### For each cell
    for (cell_i in 1:num_cells) {

      cell_name <- sprintf("cell%05d", pos_index)

      #num_reads <- 1000  #Randomize on log scale
      num_reads <- round(10**runif(1, 3, 4.5))

      #### Sample chromosome/plasmid index for each read, weighted by length
      #### (= uniform across all bases of all fasta files of the genome).
      use_chrom <- sample(1:length(allfa_str), prob = allfa_len, replace = TRUE, size = num_reads)

      #### For each chromosome/plasmid that got picked, generate its reads. All
      #### reads from this cell — regardless of which fasta they came from —
      #### share the same cell_name.
      for (cur_chrom in unique(use_chrom)) {
        get_num_reads <- sum(use_chrom == cur_chrom)

        ch <- allfa_str[cur_chrom]
        totlen <- str_length(ch)
        start_pos <- round(runif(get_num_reads, min = 1, max = totlen - insert_max_size))
        end_pos <- start_pos + round(runif(get_num_reads, min = 400, max = insert_max_size))

        all_inserts <- str_sub(ch, start_pos, end_pos)

        insert_r1 <- str_sub(all_inserts, 1, readlen)
        insert_r2_rc <- str_sub(all_inserts, str_length(all_inserts) - (readlen - 1), str_length(all_inserts))
        insert_r2 <- as.character(reverseComplement(DNAStringSet(insert_r2_rc)))

        allreads_df[[reads_i]] <- data.frame(
          name = cell_name,
          pos1 = 1,
          pos2 = 1,
          r1 = insert_r1,
          r2 = insert_r2,
          q1 = qread,
          q2 = qread,
          umi = ""
        )
        reads_i <- reads_i + 1
      }

      pos_index <- pos_index + 1
    }
  }

  ## Flatten list, save
  allreads_flat <- do.call(rbind, allreads_df)

  # 5 should be fixed!
  fname <- paste0("/husky/henriksson/atrandi/simulated5/filtered.", cur_shard, ".tirp")

  print("writing fastq")
  write.table(allreads_flat, fname, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

  print("compressing")
  system(paste0("bgzip --threads 15 -f ", fname))
  print("indexing")
  system(paste0("tabix -f -p bed ", fname, ".gz"))

}

