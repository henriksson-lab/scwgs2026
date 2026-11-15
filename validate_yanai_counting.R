#!/usr/bin/env Rscript
# Count one COMPLETE, already aligned library in an isolated directory.
# Uses the same RunYanaiGeneCounts/Zorn/local Bascet path as prep_yanai.R.
# Results stay here; the main 88-library workflow is not changed or seeded.
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
scriptDir <- dirname(normalizePath(script,mustWork=TRUE))
bulk <- normalizePath(Sys.getenv("YANAI_ROOT","/husky/henriksson/atrandi/v9_rnaseq_yanai"),mustWork=TRUE)
library <- Sys.getenv("YANAI_COUNT_VALIDATION_LIBRARY","D3_lib1")
root <- Sys.getenv("YANAI_COUNT_VALIDATION_ROOT",file.path(bulk,paste0("validation_counts_",library)))
if(dir.exists(root) && length(list.files(root,all.files=TRUE,no..=TRUE))) {
  stop("Choose a fresh YANAI_COUNT_VALIDATION_ROOT; preserving existing validation: ",root)
}
dir.create(root,recursive=TRUE,showWarnings=FALSE)
root <- normalizePath(root,mustWork=TRUE)
mapping <- read.csv(file.path(bulk,"yanai_shard_libraries.csv"))
row <- mapping[mapping$library_prefix==library,,drop=FALSE]
stopifnot(nrow(row)==1L)
originalBam <- file.path(bulk,paste0("aligned_pos.",row$group,".bam"))
stopifnot(file.exists(originalBam), file.exists(paste0(originalBam, ".bai")))
samtools <- Sys.getenv("SAMTOOLS",path.expand("~/miniconda3/bin/samtools"))
stopifnot(system2(samtools,c("quickcheck","-v",shQuote(originalBam)))==0L)
for(suffix in c("",".bai")) {
  stopifnot(file.symlink(paste0(originalBam,suffix),file.path(root,paste0("aligned_pos.1.bam",suffix))))
}
s <- readRDS(file.path(bulk,"yanai_sharding.rds"))
group <- match(library,s$group_prefix)
stopifnot(!is.na(group),group==row$group)
debstat <- list(group_prefix=s$group_prefix[group],
  list_picked_cells=s$list_picked_cells[group],
  list_group_bc_count=s$list_group_bc_count[group])
stopifnot(setequal(debstat$list_picked_cells[[1]],debstat$list_group_bc_count[[1]]$cellid))
saveRDS(debstat,file.path(root,"yanai_sharding.rds"))
rm(s);gc()
zornDir <- path.expand(Sys.getenv("YANAI_ZORN_DIR","~/github/zorn"))
for(file in c("job_general.R","job_local.R","job_slurm.R","jobscript.R",
              "bascet_file.R","zorn.R","shell.R","zorn_aggr.R","refgenome.R",
              "container.R","debarcoding.R")) source(file.path(zornDir,"R",file))
source(file.path(scriptDir,"yanai_petriseq","count_genes.R"))
instance <- getBascetDevDir(path.expand(Sys.getenv("YANAI_BASCET_DIR","~/github/bascet")),
                          targetType=Sys.getenv("YANAI_BASCET_TARGET","release"))
runner <- LocalRunner(direct=TRUE,showScript=TRUE,ncpu=Sys.getenv("YANAI_THREADS","10"),
                      mem=Sys.getenv("YANAI_MEM","24g"))
annotation <- normalizePath(Sys.getenv("YANAI_GFF","/husky/henriksson/atrandi/bwamem_ref/ling/all.gff"),mustWork=TRUE)
setwd(root)
jsonlite::write_json(list(scope="one complete library; independent counting preflight",
  library=library,original_group=row$group,original_bam=originalBam,
  main_workflow_outputs_modified=FALSE,bascet_binary=instance@bin,
  bascet_binary_md5=unname(tools::md5sum(instance@bin))),
  "provenance.json",auto_unbox=TRUE,pretty=TRUE)
RunYanaiGeneCounts(root,debstat,annotation,runner,instance)
source(file.path(scriptDir,"yanai_petriseq","verify_gene_counts.R"))
print(VerifyYanaiGeneCounts(root),row.names=FALSE)
source(file.path(scriptDir,"yanai_petriseq","verify_bam_gene_counts.R"))
VerifyBamGeneCounts(root)
capture.output(sessionInfo(),file="sessionInfo.txt")
message("Full-library count preflight passed: ",root)
