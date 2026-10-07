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


plotDirAll <- file.path("~/github/scmetag-revision/plots/all")

################################################################################
##################### Human% sc vs bulk ########################################
################################################################################


get_frac_human <- function(fname){

  bami <- str_split_fixed(system(paste("samtools idxstats",fname), intern = TRUE),"\t",4)
  bami <- as.data.frame(bami)
  colnames(bami) <- c("chr","len","mapped","unmapped")
  bami$len <- as.integer(bami$len)
  bami$mapped <- as.integer(bami$mapped)
  bami$unmapped <- as.integer(bami$unmapped)
  unmapped <- bami[bami$chr=="*",]$unmapped
  m <- sum(bami[str_starts(bami$chr, "chr"),]$mapped)
  
  data.frame(
    unmapped=unmapped,
    mapped=m
  )
}

#Get human aligned reads before filtering
human_frac_bulk <- get_frac_human("/husky/henriksson/atrandi/rawdata/251128_jyoti_mock_bulk/bulk/after_filter/bulk_sorted.bam")

## Gather all single-cell aligned reads before host filtering
human_frac_sc <- NULL
for(i in c(1:20)) {
  print(i)
  human_frac_sc <- rbind(
    human_frac_sc,
    get_frac_human(paste0("/husky/henriksson/atrandi/v9_wgs_saliva2_mda/aligned_pos.",i,".bam"))
    )
}

get_human_perc <- function(df) {
  sum(as.double(df$unmapped)) / ( sum(as.double(df$mapped)) + sum(as.double(df$unmapped)))
}


toplot_human_perc <- data.frame(
  frac=c(get_human_perc(human_frac_bulk),  get_human_perc(human_frac_sc)),
  src=c("Bulk","scMetaG")
)
toplot_human_perc

ggplot(toplot_human_perc, aes(y=frac*100, x=src)) + geom_bar(stat="identity",fill="darkblue") +
  theme_bw()+xlab("")+ylab("Fraction (%)")
ggsave(filename = file.path(plotDirAll,"humanfrac.svg"), width = 1, height = 3)


file.path(plotDirAll,"humanfrac.svg")


################################################################################
##################### Compare samples with Kraken - load #######################
################################################################################

### Read matrix, rename
mat_saliva2 <- ReadBascetCountMatrix("/husky/henriksson/atrandi/v9_wgs_saliva2_mda","kraken_counts", verbose=FALSE)      
adata <- CreateSeuratObject(mat_saliva2, project = "proj", min.cells = 0, min.features = 0)

## Add KRAKEN consensus taxonomy to metadata
kraken_taxid <- KrakenFindConsensusTaxonomy(mat_saliva2)
adata <- BascetAddMetaData(adata, kraken_taxid)



################################################################################
############ Read bulk data and compare ########################################
################################################################################

bulk_taxid <- arrow::read_feather("fig4_bulk_kraken/nohost_bulk.feather")
colnames(bulk_taxid) <- c("tax_id","cnt")


bulk_taxid <- as.data.frame(bulk_taxid)
bulk_taxid <- bulk_taxid[order(bulk_taxid$cnt, decreasing = TRUE),]

############### 
############### Bulk: make cnt
############### 

get_freq_bulk <- function(taxlevel="phylum") {
  
  taxonomizr::prepareDatabase(getAccessions=FALSE)
  df_taxid_bulk <- data.frame(taxonomizr::getTaxonomy(
    bulk_taxid$tax_id,
    desiredTaxa = c(
      "phylum", "class", "order", "family", "genus","species")
  ))
  df_taxid_bulk$taxid <- bulk_taxid$tax_id
  df_taxid_bulk$cnt <- bulk_taxid$cnt
  df_taxid_bulk <- df_taxid_bulk[!is.na(df_taxid_bulk$species),] #keep entries with known species
  
  
  ############### Sum frequency on chosen level
  bulk_freq <- sqldf(paste0("select sum(cnt) as cnt, `",taxlevel,"` as taxgroup from df_taxid_bulk group by taxgroup"))
  bulk_freq <- bulk_freq[!is.na(bulk_freq$taxgroup),]
  bulk_freq$freq <- bulk_freq$cnt / sum(bulk_freq$cnt)
  bulk_freq$protocol <- "Bulk"
  bulk_freq
}



############### 
############### Single cell, per cell: make cnt
############### 


get_freq_sc_percell <- function(taxlevel="phylum") {
  df_taxid_sc <- adata@meta.data
  
  sc_freq_percell <- sqldf(paste0("select count(`",taxlevel,"`) as cnt, `",taxlevel,"` as taxgroup from df_taxid_sc group by taxgroup"))
  sc_freq_percell <- sc_freq_percell[!is.na(sc_freq_percell$taxgroup),]
  sc_freq_percell$freq <- sc_freq_percell$cnt / sum(sc_freq_percell$cnt)
  sc_freq_percell$protocol <- "SC: by cell"
  sc_freq_percell
}



############### 
############### Single cell, per read
############### 

# "phylum", "class", "order", "family", "genus","species"

get_freq_sc_perread <- function(taxlevel="phylum") {
  
  #### Get counts per tax-id across cells
  M <- as(GetAssayData(object = adata, assay = "RNA", layer = "counts"), "TsparseMatrix")
  gather_taxid <- data.frame(
    taxid_index = M@j + 1,  #Note that indices are zero-based within Matrix package
    cnt = M@x
  )
  col_taxid <- stringr::str_remove(rownames(M),"taxid-")
  gather_taxid$tax_id <- col_taxid[gather_taxid$taxid_index]
  gather_taxid <- sqldf("select sum(cnt) as cnt, tax_id from gather_taxid group by tax_id order by cnt desc")
  
  #### Merge taxid with taxonomy table
  df_taxid_sc_perread <- data.frame(taxonomizr::getTaxonomy(
    gather_taxid$tax_id,
    desiredTaxa = c(
      "phylum", "class", "order", "family", "genus","species")
  ))
  df_taxid_sc_perread$taxid <- gather_taxid$tax_id
  df_taxid_sc_perread$cnt <- gather_taxid$cnt
  df_taxid_sc_perread
  
  ############### Compute frequency on chosen level
  onlevel <- sqldf(paste0("select sum(cnt) as cnt, `",taxlevel,"` as taxgroup from df_taxid_sc_perread group by taxgroup order by cnt desc"))
  onlevel <- onlevel[!is.na(onlevel$taxgroup),]
  onlevel$freq <- onlevel$cnt / sum(onlevel$cnt)
  onlevel$protocol <- "SC: by read"
  
  onlevel  
}  




################################################################################
############### Alluvial plot comparing SC and bulk ############################
################################################################################


bunch_other <- function(df, cutoff=0.03) {
  #Bunch together "other" for readability
  totfreq <- sqldf("select sum(freq) as sumfreq, taxgroup from df group by taxgroup")
  remove_phylum <- totfreq$taxgroup[totfreq$sumfreq<cutoff]
  df_with_other <- df
  df_with_other$taxgroup[df_with_other$taxgroup %in% remove_phylum] <- "Other"
  df_with_other <- sqldf("select taxgroup, protocol, sum(freq) as freq from df_with_other group by taxgroup, protocol")
  df_with_other
}

make_al_plot <- function(df_with_other) {
  ggplot(df_with_other,
         aes(x = protocol, stratum = taxgroup, alluvium = interaction(taxgroup, taxgroup),
             y = freq, fill = taxgroup, label = taxgroup)) +
    geom_flow(stat = "alluvium", aes(fill = taxgroup), alpha = 0.5) +
    geom_stratum(alpha = 0.8) +
    geom_text(stat = "stratum", size = 3, color = "black", fontface = "bold") +
    scale_x_discrete(expand = c(0.1, 0.1)) +
    labs(
      x     = "",
      y     = "Frequency",
      fill  = "taxgroup"
    ) +
    theme_minimal(base_size = 13) +
    #theme(plot.margin = unit(c(1,3,3,1), "cm")) +
    theme(legend.position = "none") 
}



taxlevel <- "phylum"
if(FALSE) { ############ other levels we can plot at
  taxlevel <- "class"
  taxlevel <- "order" 
  taxlevel <- "genus"
}


sc_freq_percell <- get_freq_sc_percell(taxlevel)
sc_freq_perread <- get_freq_sc_perread(taxlevel)
bulk_freq <- get_freq_bulk(taxlevel)

df <- rbind(
  sc_freq_percell,
  bulk_freq,
  sc_freq_perread
)
df$protocol <- factor(df$protocol, levels=c(
  "SC: by cell",
  "Bulk",
  "SC: by read"
))

df_with_other <- bunch_other(df, 0.05)
#df_with_other <- bunch_other(df, 0.1)

#make_al_plot(df) + ggtitle(taxlevel)
#make_al_plot(df_with_other) + ggtitle(taxlevel)
#ggsave(file.path(plotDirAll,"bulk_sc_phylum.svg"), width=3, height=6)

  
make_al_plot(df_with_other) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_fill_manual(values =  khroma::color("muted")(9))
ggsave(file.path(plotDirAll,"bulk_sc_phylum.svg"), width=3, height=6)





################################################################################
########### Supplemental fig: Tanerella abundance on bulk and SC ###############
################################################################################


taxlevel <- "species"
sc_freq_percell <- get_freq_sc_percell(taxlevel)
bulk_freq <- get_freq_bulk(taxlevel)

compare_sc_bulk <- merge(
  all=TRUE,
  data.frame(
    cnt_sc = sc_freq_percell$cnt,
    tacgroup = sc_freq_percell$taxgroup
  ),
  data.frame(
    cnt_bulk = bulk_freq$cnt,
    tacgroup = bulk_freq$taxgroup
  )
)
compare_sc_bulk[is.na(compare_sc_bulk)] <- 0


#Interesting species
compare_sc_bulk[compare_sc_bulk$cnt_bulk==0 & compare_sc_bulk$cnt_sc>100,]

compare_sc_bulk[compare_sc_bulk$tacgroup=="Tannerella sp. oral taxon 808",]
compare_sc_bulk$is_tannerella <- compare_sc_bulk$tacgroup=="Tannerella sp. oral taxon 808"
compare_sc_bulk$lab <- ""
compare_sc_bulk$lab[compare_sc_bulk$is_tannerella] <- "Tannerella"



### Plot tanerella in bulk vs single cell abundance
compare_sc_bulk <- compare_sc_bulk[order(compare_sc_bulk$is_tannerella, decreasing = FALSE),]
ggplot(compare_sc_bulk, aes(cnt_sc+1, cnt_bulk+1, label=lab, color=is_tannerella)) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10() +
  xlab("SC counts + 1")+
  ylab("Bulk counts + 1")+
  ggrepel::geom_text_repel(
    max.overlaps = 100,
    min.segment.length = 100,
    data=compare_sc_bulk[compare_sc_bulk$is_tannerella,], 
    color="black"
  ) + theme_bw() +
  theme(legend.position = "none")
ggsave(file.path(plotDirAll,"tanerella_abundance.png"), width = 5, height = 5)



