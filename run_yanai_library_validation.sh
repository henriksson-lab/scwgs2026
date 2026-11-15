#!/usr/bin/env bash
# Full-library RNA validation while the main 88-library import continues.
# Reuses an already completed debarcoding output, preserving original cell IDs.
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
bulk_root="${YANAI_ROOT:-/husky/henriksson/atrandi/v9_rnaseq_yanai}"
library="${YANAI_VALIDATION_LIBRARY:-D1_lib1}"
[[ "$library" =~ ^D[0-9]+_lib[0-9]+$ ]] || { echo 'Invalid library name' >&2; exit 1; }
validation_root="${bulk_root}/validation_full_${library}"
mkdir -p "$validation_root"
exec 9>"${validation_root}/.workflow.lock"
flock -n 9 || { echo 'Library validation is already running' >&2; exit 1; }
Rscript - "$bulk_root" "$validation_root" "$library" <<'RSCRIPT'
args <- commandArgs(trailingOnly=TRUE)
bulk <- normalizePath(args[1]); root <- normalizePath(args[2]); library <- args[3]
meta <- read.csv(file.path(bulk,"debarcoded.meta"),check.names=FALSE)
selected <- which(meta$prefix==library)
stopifnot(length(selected)==1L)
source <- file.path(bulk,paste0("debarcoded.",selected,".tirp.gz"))
stopifnot(file.exists(source),file.exists(paste0(source,".hist")))
rawdir <- file.path(root,"raw_bgzf")
dir.create(rawdir,showWarnings=FALSE)
link <- function(from,to) {
  if(!file.exists(to)) stopifnot(file.link(from,to))
  stopifnot(file.info(from)$size==file.info(to)$size)
}
for(mate in c("r1","r2")) {
  filename <- meta[[mate]][selected]
  link(file.path(meta$dir[selected],filename),file.path(rawdir,filename))
}
link(source,file.path(root,"debarcoded.1.tirp.gz"))
link(paste0(source,".hist"),file.path(root,"debarcoded.1.tirp.gz.hist"))
stored <- meta[selected,,drop=FALSE]
stored$dir <- rawdir
stored$shard <- 1L
write.csv(stored,file.path(root,"debarcoded.meta"),row.names=FALSE)
manifest <- read.delim(file.path(bulk,"fastq_manifest.tsv"),check.names=FALSE)
manifest <- manifest[manifest$library_prefix==library,,drop=FALSE]
stopifnot(nrow(manifest)==1L)
write.table(manifest,file.path(root,"fastq_manifest.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(rawdir,file.path(root,"rawdata.txt"))
RSCRIPT
export YANAI_ROOT="$validation_root"
export YANAI_RUNNER=local
export YANAI_THREADS="${YANAI_THREADS:-10}"
export YANAI_MEM="${YANAI_MEM:-24g}"
export YANAI_ARRAY_WORKERS=1
cd "$repo_dir"
Rscript prep_yanai.R
exec 8>"${bulk_root}/../.umap-analysis.lock"
flock 8
Rscript post_yanai.R
Rscript - "$validation_root" <<'RSCRIPT'
suppressPackageStartupMessages(library(Seurat))
root <- commandArgs(trailingOnly=TRUE)[1]
mapping <- read.csv(file.path(root,"yanai_gene_count_libraries.csv"))
selection <- read.csv(file.path(root,"umap/gene_cell_selection.csv"))
coordinates <- read.csv(file.path(root,"umap/gene_umap.csv"))
adata <- readRDS(file.path(root,"cache_adata_gene.RDS"))
embedding <- Embeddings(adata,"gene_umap")
stopifnot(nrow(mapping)==1L,nrow(selection)==mapping$retained_cells,
          setequal(selection$cell_id[selection$retained],colnames(adata)),
          identical(coordinates$cell_id,rownames(embedding)),
          all(is.finite(embedding)),
          isTRUE(all.equal(as.matrix(coordinates[,c("UMAP1","UMAP2")]),
                           embedding,tolerance=1e-10,check.attributes=FALSE)))
jsonlite::write_json(list(status="verified",scope="one complete library",
                         library=mapping$library_prefix,valid_barcodes=nrow(selection),
                         umap_cells=ncol(adata)),file.path(root,"umap/validation.json"),
                     auto_unbox=TRUE,pretty=TRUE)
RSCRIPT
