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
#  small_runner <- bascetRunner.default
} else {
  bascetInstance.default <- getBascetDevDir("/home/m/mahogny/github/bascet")
#  bascetInstance.default <- getBascetBinary()
  bascetRunner.default <- SlurmRunner(account="hpc2n2026-186", ncpu="10", mem="35g", logTime=TRUE, verbose=TRUE, time="0-14:00:00")  #henlab
#  small_runner <- SlurmRunner(bascetRunner.default, ncpu="2", mem="20g")
}

# skylake nodes: ~5GB/core
# zen3: ~6.5GB/core
# zen4: ~13GB/core
# largemem nodes: ~4GB/core
# So depending on how much corehours you want to waste (if the code is serial only) you can figure it out from that list. It does NOT apply to GPU nodes though, only the CPU ones.


################################################################################
################## Alignment to tannerella #####################################
################################################################################

# https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(25)00254-1
# Bowtie2 VSL T2T-CHM13  best

if(is_beagle) {
  bwaref <- "/husky/henriksson/atrandi/bwamem_ref/tanerella/all.fa"
} else {
  bwaref <- "/home/m/mahogny/mystore/atrandi/bwamem_ref/tanerella/all.fa"
}

### Perform alignment
BascetAlignToReference(
  bascetRoot,
  inputName="filtered_nohost",
  outputName="aligned_tannerella",
  useReference = bwaref
)




### Make a count-per-chromosome count table
BascetCountChrom(
  bascetRoot,
  inputName="aligned_tannerella_cell",
  outputName="countchrom_tannerella",
  removeMultimapper = FALSE,
  minMatching = 40
)



################################################################################
################## get tannerella reads filtering ##############################
################################################################################


BascetFilterAlignment(
    bascetRoot,
    inputName="aligned_tannerella_cell",
    outputName="aligned_justtannerella_cell",
    keepMapped=TRUE  
#,
#    runner=SlurmRunner(bascetRunner.default, ncpu="5", mem="10g")
)

#Convert filtered data to TIRP
BascetMapTransform(
  bascetRoot,
  inputName="aligned_justtannerella_cell",
  outputName = "aligned_justtannerella_tirp",
  outFormat="tirp.gz"
)

#Convert filtered data to FQ
BascetMapTransform(
  bascetRoot,
  inputName="aligned_justtannerella_cell",
  outputName = "aligned_justtannerella_fq",
  outFormat="R1.fq.gz"
)

