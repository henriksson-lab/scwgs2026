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
} else {
  library(Zorn)
}

## Check where the script is running
is_beagle <- Sys.info()["nodename"]=="beagle"


setwd(bascetRoot)


if(is_beagle) {
  bascetInstance.default <- getBascetDevDir("/home/mahogny/github/bascet")
  bascetRunner.default <- LocalRunner(direct = TRUE, showScript=FALSE, mem="50g")
} else {
  bascetInstance.default <- getBascetDevDir("/home/m/mahogny/github/bascet")
#  bascetInstance.default <- getBascetBinary()
  bascetRunner.default <- SlurmRunner(account="hpc2n2026-186", ncpu="10", mem="35g", logTime=TRUE, verbose=TRUE, time="0-14:00:00")  #henlab
}

# skylake nodes: ~5GB/core
# zen3: ~6.5GB/core
# zen4: ~13GB/core
# largemem nodes: ~4GB/core
# So depending on how much corehours you want to waste (if the code is serial only) you can figure it out from that list. It does NOT apply to GPU nodes though, only the CPU ones.


################################################################################
################## Alignment workflow ##########################################
################################################################################

#10 species picked from NCBI. not the same as ATCC
if(is_beagle) {
  bwaref <- "/husky/henriksson/atrandi/bwamem_ref/simulated_ref10/all.fa"
} else {
  bwaref <- "/home/m/mahogny/mystore/atrandi/bwamem_ref/simulated_ref10/all.fa"
}

### Perform alignment
BascetAlignToReference(
  bascetRoot,
  aligner="BWAMEM2",
  useReference = bwaref  
)


### Make a count-per-chromosome count table
BascetCountChrom(
  bascetRoot,
  minMatching = 40
)

################################################################################
################## Count sketch workflow #######################################
################################################################################


BascetRunCountsketch(
  bascetRoot,
  inputName = "filtered"
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
  inputName="filtered",
  useKrakenDB = kraken_db
)


################################################################################
################## Informative KMER ############################################
################################################################################



### Compute minhashes for each cell
BascetComputeMinhash(
  bascetRoot
)


### Gather minhashes into a single histogram
BascetMakeMinhashHistogram(
  bascetRoot
)

### Pick KMERs to detect
print("pick min kmers")
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

