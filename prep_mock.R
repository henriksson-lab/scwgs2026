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
  bascetRunner.default <- LocalRunner(direct = TRUE, showScript=FALSE, mem="50g")

  small_runner <- bascetRunner.default

} else {
  bascetInstance.default <- getBascetDevDir("/home/m/mahogny/github/bascet")
#  bascetInstance.default <- getBascetBinary()
  bascetRunner.default <- SlurmRunner(account="hpc2n2026-186", ncpu="10", mem="40g", logTime=TRUE, verbose=TRUE, time="0-14:00:00")  #henlab

  small_runner <- SlurmRunner(bascetRunner.default, ncpu="2")
  #big_runner <- SlurmRunner(bascetRunner.default, mem="80gb", ncpu="15")
  big_runner <- SlurmRunner(bascetRunner.default, mem="60gb", ncpu="10")

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
    numOutputShards = 40
  )

}





################################################################################
################## Alignment workflow ##########################################
################################################################################

#if(is_beagle) {
#  bwaref <- "/husky/henriksson/atrandi/bwa_ref/atcc/ref10_xp/ref10_xp.fa"
#  small_runner <- bascetRunner.default
#} else {
#  bwaref <- "/home/m/mahogny/mystore/atrandi/bwa_ref/atcc/ref10_xp/ref10_xp.fa"
#  small_runner <- SlurmRunner(bascetRunner.default, ncpu="2")
#}


if(is_beagle) {
  bwaref <- "/husky/henriksson/atrandi/bwamem_ref/atcc/ref10_xp/ref10_xp.fa"
} else {
  bwaref <- "/home/m/mahogny/mystore/atrandi/bwamem_ref/atcc/ref10_xp/ref10_xp.fa"
}



### Perform alignment
BascetAlignToReference(
  bascetRoot,
  useReference = bwaref  
)


### Make a count-per-chromosome count table -- keep xantham in so we can QC
BascetCountChrom(
  bascetRoot,
  minMatching = 40,
  runner = small_runner
)

################################################################################
################## xantham filtering ###########################################
################################################################################



if(is_beagle) {
  bwaref <- "/husky/henriksson/atrandi/bwamem_ref/phi_xantham/all.fa"
} else {
  bwaref <- "/home/m/mahogny/mystore/atrandi/bwamem_ref/phi_xantham/all.fa"
}

### Perform alignment
BascetAlignToReference(
  bascetRoot,
  useReference = bwaref,
  outputName = "xp_aligned"
)


BascetFilterAlignment(
    bascetRoot,
    inputName="xp_aligned_cell",
    outputName="xp_aligned_cell_nohost", ##### makes bam; would be nice to get tirp right away
    keepMapped=FALSE
)

#Convert filtered data to TIRP
BascetMapTransform(
  bascetRoot,
  inputName="xp_aligned_cell_nohost",
  outputName = "filtered_nohost",
  outFormat="tirp.gz"
)



################################################################################
################## Count sketch workflow #######################################
################################################################################


BascetRunCountsketch(
  bascetRoot,
  inputName = "filtered_nohost"
)



#stop(666)

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
################## To fastq, for metagraph etc  ################################
################################################################################


if(FALSE){
### Get reads in fastq format
BascetMapTransform(
  bascetRoot,
  inputName="filtered_nohost",
  outputName="asfq",
  outFormat="R1.fq.gz",
  runner=SlurmRunner(bascetRunner.default, ncpu="4", mem="20g")  #in and out processes
)
}



################################################################################
################## Informative KMER ############################################
################################################################################

### Compute minhashes for each cell
BascetComputeMinhash(
  bascetRoot,
  inputName = "filtered_nohost"
)


### Gather minhashes into a single histogram
BascetMakeMinhashHistogram(
  bascetRoot
)

### Pick KMERs to detect
print("pick informative KMERs")
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

 BascetQueryFq( ########### better name needed!!  DetectKMERs? CountKMERs?
   bascetRoot,
   useKMERs=useKMERs
 )
}



################################################################################
################## QC ##########################################################
################################################################################


### quast QC
BascetMapCellQUAST(
#BascetMapCell(
  bascetRoot,
  inputName = "contigs",
)


################################################################################
################## De novo Assembly ############################################
################################################################################


### Assemble all genomes
BascetMapCellSKESA(
  bascetRoot,
  inputName = "filtered_nohost",
  maxReadsPerCell = 5000000,
  runner = big_runner
)

### 
BascetMapCellGECCO(
  bascetRoot
)

### FastQC
BascetMapCellFASTQC(
  bascetRoot,
  inputName = "filtered_nohost"
)




if(FALSE) {

### requires blast. not ready

### AMRfinder QC
BascetMapCellAMRfinder(
  bascetRoot,
  inputName = "contigs",
  db = "/husky/henriksson/atrandi/amrfinder_db"
)

}


###########################

stop(666)



print(2)
### abricate QC
BascetMapCellAbricate(
  bascetRoot,
  inputName = "contigs",
  args = list(DATABASE_DIR="ncbi")
)

print(3)

