#!/usr/bin/env Rscript
# Standard R/Zorn DNA UMAP workflow; see post_dna.R and post_unified.R.
arg <- grep("^--file=",commandArgs(),value=TRUE)
file <- tryCatch(sys.frame(1)$ofile,error=function(e) NULL)
if(is.null(file)) file <- sub("^--file=","",arg[1])
scriptDir <- dirname(normalizePath(file,mustWork=TRUE))
Sys.setenv(DNA_DATASET="LING")
source(file.path(scriptDir,"post_dna.R"))
