#!/usr/bin/env Rscript
# Recreate the pinned reference without changing matching existing files.
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
here <- dirname(normalizePath(script,mustWork=TRUE))
args <- commandArgs(trailingOnly=TRUE)
root <- "/husky/henriksson/atrandi/bwamem_ref/ling"
if(length(args)) {
  if(length(args)!=2L || args[1]!="--root") stop("Usage: prepare.R [--root directory]")
  root <- args[2]
}
dir.create(file.path(root,"source"),recursive=TRUE,showWarnings=FALSE)
root <- normalizePath(root,mustWork=TRUE)
sha <- function(path) digest::digest(file=path,algo="sha256",serialize=FALSE)
publish <- function(temp,path) {
  if(file.exists(path)) {
    if(sha(temp)!=sha(path)) stop("Refusing to replace different reference: ",path)
  } else if(!file.link(temp,path)) stop("Cannot publish ",path)
  unlink(temp)
}
refs <- jsonlite::read_json(file.path(here,"assemblies.json"))
manifest <- read.delim(file.path(here,"source_manifest.tsv"),colClasses="character")
for(i in seq_len(nrow(manifest))) {
  row <- manifest[i,]
  path <- file.path(root,"source",row$filename)
  valid <- function(p) file.exists(p) && file.info(p)$size==as.numeric(row$bytes) &&
    unname(tools::md5sum(p))==row$md5
  if(!file.exists(path)) {
    part <- paste0(path,".part")
    status <- system2("curl",shQuote(c("--fail","--location","--silent","--show-error",
      "--retry","5","--continue-at","-","--output",part,row$url)))
    if(status!=0L || !valid(part)) stop("Download/checksum failure: ",part)
    if(!file.link(part,path)) stop("Cannot publish ",path)
    unlink(part)
  }
  if(!valid(path)) stop("NCBI size/MD5 mismatch: ",path)
  message("Verified ",row$filename)
}
sequences <- numeric()
species <- character()
fastaTemp <- tempfile("reference-fasta-",tmpdir=root)
gffTemp <- tempfile("reference-gff-",tmpdir=root)
on.exit.cleanup <- c(fastaTemp,gffTemp)
tryCatch({
  fasta <- file(fastaTemp,"wt")
  tryCatch(for(ref in refs) {
    con <- gzfile(file.path(root,"source",paste0(ref$stem,"_genomic.fna.gz")),"rt")
    lines <- tryCatch(readLines(con),finally=close(con))
    headers <- which(startsWith(lines,">"))
    stopifnot(length(headers)>0L,headers[1]==1L)
    ids <- sub(" .*","",substring(lines[headers],2L))
    stopifnot(!anyDuplicated(ids),!any(ids %in% names(sequences)))
    ends <- c(headers[-1]-1L,length(lines))
    lengths <- vapply(seq_along(headers),function(i) {
      stopifnot(ends[i]>headers[i])
      sum(nchar(trimws(lines[seq.int(headers[i]+1L,ends[i])])))
    },numeric(1))
    sequences <- c(sequences,setNames(lengths,ids))
    species <- c(species,setNames(rep(ref$key,length(ids)),ids))
    writeLines(lines,fasta,useBytes=TRUE)
  },finally=close(fasta))
  allFeatures <- character()
  geneIds <- character()
  gffIds <- character()
  geneNames <- 0L
  gff <- file(gffTemp,"wt")
  tryCatch({
    writeLines("##gff-version 3",gff)
    for(ref in refs) {
      writeLines(paste("# source_assembly",ref$accession,ref$species),gff)
      con <- gzfile(file.path(root,"source",paste0(ref$stem,"_genomic.gff.gz")),"rt")
      lines <- tryCatch(readLines(con),finally=close(con))
      stopifnot(!any(startsWith(lines,"##FASTA")))
      lines <- lines[!startsWith(lines,"##gff-version")]
      records <- strsplit(lines[!startsWith(lines,"#")],"\t",fixed=TRUE)
      stopifnot(all(lengths(records)==9L))
      records <- do.call(rbind,records)
      start <- as.numeric(records[,4]);end <- as.numeric(records[,5])
      stopifnot(all(records[,1] %in% names(sequences)),all(is.finite(start)),
        all(is.finite(end)),all(start>=1 & start<=end),
        all(end<=sequences[records[,1]]))
      allFeatures <- c(allFeatures,records[,3])
      attrs <- strsplit(records[records[,3]=="gene",9],";",fixed=TRUE)
      attribute <- function(parts,key) {
        found <- parts[startsWith(parts,paste0(key,"="))]
        if(length(found)!=1L) return("")
        utils::URLdecode(substring(found,nchar(key)+2L))
      }
      geneIds <- c(geneIds,vapply(attrs,attribute,character(1),key="locus_tag"))
      gffIds <- c(gffIds,vapply(attrs,attribute,character(1),key="ID"))
      geneNames <- geneNames+sum(nzchar(vapply(attrs,attribute,character(1),key="Name")))
      writeLines(lines,gff,useBytes=TRUE)
    }
  },finally=close(gff))
  stopifnot(all(nzchar(geneIds)),!anyDuplicated(geneIds),
            all(nzchar(gffIds)),!anyDuplicated(gffIds),length(geneIds)==11800L,
            geneNames==11800L,sum(sequences)==11774727)
  publish(fastaTemp,file.path(root,"all.fa"))
  publish(gffTemp,file.path(root,"all.gff"))
  for(name in c("assemblies.json","source_manifest.tsv")) {
    temp <- tempfile(tmpdir=root)
    stopifnot(file.copy(file.path(here,name),temp))
    publish(temp,file.path(root,name))
  }
  report <- list(paper_doi="10.1038/s41586-023-06974-w",assemblies=refs,
    contigs=lapply(names(sequences),function(id) list(id=id,length=unname(sequences[id]),
                                                    species=unname(species[id]))),
    total_bases=sum(sequences),feature_counts=as.list(table(allFeatures)),
    unique_gene_locus_tags=length(geneIds),genes_with_Name=geneNames,
    bascet_countfeature=list(useFeature="gene",attrGeneId="locus_tag",attrGeneName="Name"),
    all_fa_sha256=sha(file.path(root,"all.fa")),all_gff_sha256=sha(file.path(root,"all.gff")))
  path <- file.path(root,"reference.json")
  if(file.exists(path)) {
    old <- jsonlite::read_json(path)
    stopifnot(old$all_fa_sha256==report$all_fa_sha256,
              old$all_gff_sha256==report$all_gff_sha256,
              old$unique_gene_locus_tags==report$unique_gene_locus_tags)
  } else jsonlite::write_json(report,path,auto_unbox=TRUE,pretty=TRUE)
  message("Ready: ",length(sequences)," contigs; ",sum(sequences)," bases; ",
          length(geneIds)," unique genes. Existing matching reference preserved.")
},finally=unlink(on.exit.cleanup))
