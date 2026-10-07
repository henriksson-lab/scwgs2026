bascetRoot <- getwd()



if(TRUE){
  setwd("~/github/zorn")
  source("R/job_general.R")
  source("R/job_local.R")
  source("R/job_slurm.R")
  source("R/jobscript.R")
  source("R/bascet_file.R")
  source("R/zorn.R")
  source("R/shell.R")
  source("R/zorn_aggr.R")
  source("R/count_kmer.R")
  source("R/countsketch.R")
  source("R/refgenome.R")
  source("R/kraken.R")
  source("R/container.R")
  source("R/ext_tools.R")
  source("R/debarcoding.R")
  source("R/sra_import.R")
} else {
  library(Zorn)
}



## Check where the script is running
is_beagle <- Sys.info()["nodename"]=="beagle"


setwd(bascetRoot)


if(is_beagle) {
  bascetInstance.default <- getBascetDevDir("/home/mahogny/github/bascet")
#  bascetRunner.default <- LocalRunner(direct = TRUE, showScript=FALSE, mem="50g")
  bascetRunner.default <- LocalRunner(direct = TRUE, showScript=FALSE, mem="40g")

  small_runner <- bascetRunner.default
  big_runner <- bascetRunner.default

} else {
  bascetInstance.default <- getBascetDevDir("/home/m/mahogny/github/bascet")
#  bascetInstance.default <- getBascetBinary()
#  bascetRunner.default <- SlurmRunner(account="hpc2n2026-186", ncpu="10", mem="50g", logTime=TRUE, verbose=TRUE, time="0-14:00:00")  #henlab
  bascetRunner.default <- SlurmRunner(account="hpc2n2026-186", ncpu="12", mem="60g", logTime=TRUE, verbose=TRUE, time="0-14:00:00")  #henlab
#  bascetRunner.default <- SlurmRunner(account="hpc2n2026-186", ncpu="12", mem="80g", logTime=TRUE, verbose=TRUE, time="0-14:00:00")  #henlab


  small_runner <- SlurmRunner(bascetRunner.default, ncpu="2")
  big_runner <- SlurmRunner(account="hpc2n2026-186", ncpu="12", mem="100g", logTime=TRUE, verbose=TRUE, time="0-14:00:00")  #henlab

}

# skylake nodes: ~5GB/core
# zen3: ~6.5GB/core
# zen4: ~13GB/core
# largemem nodes: ~4GB/core
# So depending on how much corehours you want to waste (if the code is serial only) you can figure it out from that list. It does NOT apply to GPU nodes though, only the CPU ones.


################################################################################
############################### Debarcode ######################################
################################################################################

### Only perform if there is raw data (skipped if simulated data)
if(file.exists("rawdata.txt")) {
  rawmeta_dir <- readLines("rawdata.txt")
  rawmeta <- DetectRawFileMeta(rawmeta_dir)

  #Decide on barcode error tolerance
  barcode_error <- NULL
  if(file.exists("barcode_error.txt")) {
    barcode_error <- readLines("barcode_error.txt") ################################### TODO note not used below!!
  }

  BascetDebarcode(
    bascetRoot,
    rawmeta,
    chemistry="atrandi-wgs"
  )

  ### Figure out which drops to keep (roughly)
  debstat <- PrepareSharding(
    bascetRoot,
    inputName="debarcoded",
    minQuantile=0.95
  )

  DebarcodedKneePlot(debstat, filename = "kneeplot.pdf")

  ### Shardify i.e. divide into multiple sets of files for parallel processing.
  BascetShardify(
    debstat,
    numOutputShards = 20   ############ was 40
  )

}





################################################################################
################## kraken workflow - unfiltered ################################
################################################################################


if(is_beagle) {
  kraken_db <- "/data/henlab/kraken/standard-8"
} else {
  kraken_db <- "/home/m/mahogny/mystore/atrandi/kraken_ref/standard-8"
}

### Run Kraken2 on each cell
BascetRunKraken(
  bascetRoot,
  inputName="filtered",
  outputName="humankraken",
  useKrakenDB = kraken_db
)



################################################################################
################## Alignment workflow ##########################################
################################################################################

# https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(25)00254-1
# Bowtie2 VSL T2T-CHM13  best

if(is_beagle) {
  bwaref <- "/husky/henriksson/atrandi/bwamem_ref/human_phi_xantham/all.fa"
} else {
  bwaref <- "/home/m/mahogny/mystore/atrandi/bwamem_ref/human_phi_xantham/all.fa"
}

### Perform alignment
BascetAlignToReference(
  bascetRoot,
  useReference = bwaref,
  runner = big_runner
)



################################################################################
################## Host filtering ##############################################
################################################################################


BascetFilterAlignment(
    bascetRoot,
    inputName="aligned_cell",
    outputName="aligned_cell_nohost", ##### makes bam; would be nice to get tirp right away
    keepMapped=FALSE
)



BascetMapTransform(
  bascetRoot,
  inputName="aligned_cell_nohost",
  outputName = "filtered_nohost",
  outFormat="tirp.gz"
)




################################################################################
################## Alignment-filter workflow ##########################################  #TODO map out human and other reads
################################################################################

#human reads were removed by mapping the high-quality reads against the human genome (hg19) using Bowtie2 (ref. 53) (v2.4.4). After removal of low-quality (quality score <20) and human reads
#nohuman: https://github.com/mbhall88/nohuman
#bowtie2: https://www.metagenomics.wiki/tools/short-read/remove-host-sequences  


#align using normal command. then another command to filter and keep? as before - keeps aligned files for troubleshooting.
#but don't use samtools?


# https://pmc.ncbi.nlm.nih.gov/articles/PMC11878760/


# !!! https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(25)00254-1 !!!
#  Bowtie2 VSL T2T-CHM13  best
# T2T-CHM13v2.0  need to use; combine xantham + phiX

# bowtie2-build toy_dataset_contig_for_mapping.fasta toy_dataset_contig_for_mapping.btindex
# bowtie2 -x toy_dataset_contig_for_mapping.btindex -f -U toy_dataset_reads_for_mapping.fasta -S toy_dataset_mapped_species1.sam

# -f 12  Extract only ( -f ) alignments with both reads unmapped
# a tutorial also suggests -F 256    # Do not (  -F  ) extract alignments which are: <not primary alignment>

################################################################################
################## Count sketch workflow #######################################
################################################################################


BascetRunCountsketch(
  bascetRoot,
  sketchSize=16384,
  inputName="filtered_nohost"
)


################################################################################
################## kraken workflow #############################################
################################################################################


if(is_beagle) {
  kraken_db <- "/data/henlab/kraken/standard-8"
} else {
  kraken_db <- "/home/m/mahogny/mystore/atrandi/kraken_ref/standard-8"
}

### Run Kraken2 on each cell
BascetRunKraken(
  bascetRoot,
  inputName="filtered_nohost",
  useKrakenDB = kraken_db
)





################################################################################
################## custom kraken workflow ###################################### abc_erysi database
################################################################################

if(FALSE) {

if(is_beagle) {
  kraken_db_chetan_abce <- "/husky/henriksson/atrandi/from_chetan/chetan_databases/abc_erysi"  #abc_erysi  refseq database
} else {
  #kraken_db_chetan_abce <- "/home/m/mahogny/mystore/atrandi/kraken_ref/standard-8"
}

### Run Kraken2 on each cell
BascetRunKraken(
  bascetRoot,
  inputName="filtered_nohost",
  outputRawName="abce_kraken_raw",
  outputName="abce_kraken_counts",
  useKrakenDB = kraken_db_chetan_abce
)

}


################################################################################
################## custom kraken workflow ###################################### MAG database
################################################################################

if(FALSE) {

if(is_beagle) {
  kraken_db_chetan <- "/husky/henriksson/atrandi/from_chetan/chetan_databases/kraken2_db_uhgg_v2.0.1"  #MAG database
} else {
  #kraken_db_chetan <- "/home/m/mahogny/mystore/atrandi/kraken_ref/standard-8"
}

### Run Kraken2 on each cell
BascetRunKraken(
  bascetRoot,
  inputName="filtered_nohost",
  outputRawName="custom_kraken_raw",
  outputName="custom_kraken_counts",
  useKrakenDB = kraken_db_chetan
)

}


################################################################################
################## To fastq, for metagraph etc  ################################
################################################################################


### Get reads in fastq format
BascetMapTransform(
  bascetRoot,
  inputName="filtered_nohost",
  outputName="asfq_nohost",
  outFormat="R1.fq.gz",
  runner=SlurmRunner(bascetRunner.default, ncpu="4")  #in and out processes
)


#  inputName="filtered_nohost",		#is SE???
#  outputName="asfq_nohost",		#is PE


################################################################################
################## Informative KMER ############################################
################################################################################

### Compute minhashes for each cell
BascetComputeMinhash(
  bascetRoot,
  inputName = "filtered_nohost",
)


### Gather minhashes into a single histogram
BascetMakeMinhashHistogram(
  bascetRoot
)

### Pick KMERs to detect
p_minh <- file.path(bascetRoot,"minhash_hist.csv")
p_use_kmers <- file.path(bascetRoot,"use_kmers.txt")
if(file.exists(p_minh)) {
  kmerHist <- BascetReadMinhashHistogram(bascetRoot)
  if(!file.exists(p_use_kmers)) {
    #print(dim(kmerHist))    
    set.seed(666)
    kmerHist$rand_index <- sample(1:nrow(kmerHist)) #add a tie breaker
    kmerHist <- kmerHist[order(kmerHist$cnt, kmerHist$rand_index, decreasing=TRUE),]

    #Pick the most common KMERs    
    picked_kmers <- kmerHist$kmer[1:100000]
    writeLines(picked_kmers, p_use_kmers)
  }
}


### Build count table by looking up selected KMERs in per-cell KMER databases
if(file.exists("use_kmers.txt")){
 useKMERs <- readLines("use_kmers.txt")

 BascetQueryFq(
   bascetRoot,
   useKMERs=useKMERs
 )
}

################################################################################
################## De novo Assembly ############################################
################################################################################



### Assemble all genomes
BascetMapCellSKESA(
  bascetRoot,
  inputName = "filtered_nohost",
  maxReadsPerCell = 5000000, #was 1M,
  runner = big_runner
)


### 
BascetMapCellGECCO(
  bascetRoot 
)

### FastQC
BascetMapCellFASTQC(
  bascetRoot,
  inputName = "filtered_nohost"  # todo make mapcell work on fastq as input too
)


### quast QC
BascetMapCellQUAST(
  bascetRoot,
  inputName  = "contigs",
  outputName = "quast"
)

### abricate QC
#BascetMapCell(
#  bascetRoot,
#  withfunction = "_abricate",
#  inputName = "skesa",
#  outputName = "abricate",
#  args = list(DATABASE_DIR="ncbi")
#)

