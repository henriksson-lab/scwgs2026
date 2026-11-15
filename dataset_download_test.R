#!/usr/bin/env Rscript
# Real shell/R integration: byte preservation, resume and checksum rejection.
script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
repo <- dirname(normalizePath(script,mustWork=TRUE))
testRoot <- tempfile("download R fixture ",tmpdir="/tmp")
dir.create(testRoot)
root <- file.path(testRoot,"yanai")
dir.create(file.path(root,"raw_gzip"),recursive=TRUE)
streams <- list(charToRaw("@pair1\nACGTN\n+\nIIIII\n@pair2\nTGCAA\n+\nJJJJJ\n"),charToRaw("@pair1\nGGGCC\n+\nHHHHH\n@pair2\nAACCT\n+\nFFFFF\n"))
inputs <- file.path(root,"raw_gzip",paste0("SRR123_",1:2,".fastq.gz"))
for(i in 1:2) {con<-gzfile(inputs[i],"wb");writeBin(streams[[i]],con);close(con)}
manifest <- file.path(testRoot,"fastq_manifest.tsv")
row <- data.frame(run_accession="SRR123",library_prefix="D1_lib1",fastq_ftp="fixture/SRR123_1.fastq.gz;fixture/SRR123_2.fastq.gz",fastq_md5=paste(unname(tools::md5sum(inputs)),collapse=";"),fastq_bytes=paste(file.info(inputs)$size,collapse=";"))
write.table(row,manifest,sep="\t",row.names=FALSE,quote=FALSE)
run <- function(script,args,expected=0L) {
  status <- suppressWarnings(system2("bash",shQuote(c(file.path(repo,script),args)),stdout=TRUE,stderr=TRUE))
  exit <- attr(status,"status");if(is.null(exit)) exit<-0L
  if(exit!=expected) stop(paste(status,collapse="\n"))
  cat(paste(status,collapse="\n"),"\n")
}
args <- c("--root",root,"--manifest",manifest,"--workers","2","--threads","1")
run("run_yanai_download.sh",args)
outputs <- file.path(root,"raw_bgzf",paste0("D1_lib1_S1_L001_R",1:2,"_001.fastq.gz"))
initial <- tools::md5sum(outputs)
for(i in 1:2) {con<-gzfile(outputs[i],"rb");actual<-readBin(con,"raw",n=10000);close(con);stopifnot(identical(actual,streams[[i]]))}
run("run_yanai_download.sh",args)
stopifnot(identical(initial,tools::md5sum(outputs)))
locked <- suppressWarnings(system2("flock",shQuote(c("-n",file.path(root,".download.lock"),
  "bash",file.path(repo,"run_yanai_download.sh"),args)),stdout=TRUE,stderr=TRUE))
stopifnot(identical(attr(locked,"status"),1L))
sourceBytes <- readBin(inputs[1],"raw",n=file.info(inputs[1])$size)
con<-file(inputs[1],"ab");writeBin(as.raw(0),con);close(con)
run("run_yanai_download.sh",args,1L)
writeBin(sourceBytes,inputs[1])
packed <- readBin(outputs[1],"raw",n=file.info(outputs[1])$size)
badPacked <- packed
badPacked[length(badPacked)-5L] <- as.raw(bitwXor(as.integer(badPacked[length(badPacked)-5L]),1L))
writeBin(badPacked,outputs[1])
run("run_yanai_download.sh",args,1L)
writeBin(packed,outputs[1])
lingRoot <- file.path(testRoot,"ling")
lingManifest <- file.path(testRoot,"manifest.tsv")
r <- data.frame(record="1",library="feces_deep",barcode_a="ACGT",role="barcode",filename="reads.fastq.gz",bytes=file.info(inputs[2])$size,md5=unname(tools::md5sum(inputs[2])),url=paste0("file://",inputs[2]))
write.table(r,lingManifest,sep="\t",row.names=FALSE,quote=FALSE)
# CURL's file URL requires escaped spaces; the shell still receives one argument.
r$url <- paste0("file://",gsub(" ","%20",inputs[2],fixed=TRUE))
write.table(r,lingManifest,sep="\t",row.names=FALSE,quote=FALSE)
args <- c("--root",lingRoot,"--manifest",lingManifest,"--workers","2")
run("run_ling_download.sh",args)
run("run_ling_download.sh",args)
stopifnot(unname(tools::md5sum(file.path(lingRoot,"source/feces_deep/reads.fastq.gz")))==r$md5)
cat("PASS: parallel mates preserve every inflated byte; reuse preserves compressed bytes; concurrent invocation, source corruption and BGZF CRC corruption rejected; checksum-verified curl download and reuse pass.\nFixture:",testRoot,"\n")
