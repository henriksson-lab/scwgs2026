################## Fig2: distribution of reads across species and gram pos/neg
library(Zorn)
library(Seurat)
library(Signac)
library(sqldf)
library(stringr)
library(ggplot2)
library(Matrix)
library(GenomicRanges)
library(future)
library(ggplot2)
library(ggalluvial)
library(ggrastr)



list_color_mock <- c(khroma::color("muted")(n = 9),"#DDDDDD")


##############################
### Novaseq pools

load_frac <- function(fname,poolname) {
  adata_mock1 <- readRDS(file.path(fname, "cache_adata_align.RDS"))
  DefaultAssay(adata_mock1) <- "species_cnt"
  
  onedf <- rowSums(adata_mock1@assays$species_cnt$counts)
  onefrac <- data.frame(
    species=names(onedf),
    frac=onedf/sum(onedf),
    pool=poolname
  )  
  onefrac
}



pool1 <- load_frac(
  "/husky/henriksson/atrandi/v9_wgs_mock1_mda",
  "mock1_MDA"
)

pool3 <- load_frac(
  "/husky/henriksson/atrandi/v9_wgs_mock3_mda",
  "mock3_MDA"
)

pool4 <- load_frac(
  "/husky/henriksson/atrandi/v9_wgs_mock4_mda",
  "mock4_MDA"
)

pool5 <- load_frac(
  "/husky/henriksson/atrandi/v9_wgs_mock5_pta",
  "mock5_PTA"
)


##############################
### Bulk MiSeq sample: species fractions from samtools idxstats
mapSeq2strain <- read.csv("/husky/henriksson/atrandi/bwa_ref/atcc/map_seq2strain.csv")
mapSeq2strain <- mapSeq2strain[mapSeq2strain$strain!="phix" & mapSeq2strain$strain!="Xanthomonas campestris",]

load_frac_idxstats <- function(fname,poolname) {
  idxstats <- read.table(
    text=system(paste("samtools idxstats", fname), intern=TRUE),
    sep="\t", col.names=c("id","len","mapped","unmapped")
  )
  idxstats <- merge(idxstats, mapSeq2strain)
  onedf <- tapply(idxstats$mapped, idxstats$strain, sum)
  data.frame(
    species=names(onedf),
    frac=as.double(onedf/sum(onedf)),
    pool=poolname
  )
}

pool2_miseq <- load_frac_idxstats(
  "/husky/henriksson/atrandi/v9_wgs_mock2_mda_miseq/mock2_miseq.sorted.bam",
  "mock2_MDA_miseq"
)


alldf <- rbind(
  pool1,
  pool3,
  pool4,
  pool5,
  pool2_miseq
)



alldf$gram <- "pos"
alldf$gram[stringr::str_detect(alldf$species,"coli")] <- "neg"
alldf$gram[stringr::str_detect(alldf$species,"sphaer")] <- "neg"



##############################
############################## Supplemental note: Fraction for all (selected subset in fig2)
##############################

ggplot(alldf, aes(pool, frac*100, fill=species)) + 
  geom_bar(stat="identity") +
  theme_bw() + 
  scale_fill_manual(values = list_color_mock) + 
  scale_x_discrete(limits=rev)+
  ylab("Fraction (%)")+
  xlab("Sample")+
  coord_flip() 
ggsave(file.path(plotDirAll,"supnote_fractions_all.pdf"), width = 8, height = 2)


##############################
############################## Supplemental note: Gram pos/neg for all (selected subset in fig2)
##############################


### Mean fraction per species within each gram class, normalized to sum to 1 per pool
gramdf <- sqldf("select pool, gram, avg(frac) as frac from alldf group by pool, gram")
gramdf$frac <- gramdf$frac / ave(gramdf$frac, gramdf$pool, FUN=sum)

ggplot(gramdf, aes(pool, frac*100, fill=gram)) +
  geom_bar(stat="identity") +
  theme_bw() + 
  scale_x_discrete(limits=rev)+
  ylab("Fraction (%)")+
  xlab("Sample")+
  coord_flip() +
  scale_fill_manual(values=c("#bb5566", "#004488")) 
ggsave(file.path(plotDirAll,"supnote_fractions_gram.pdf"), width = 6, height = 2)

