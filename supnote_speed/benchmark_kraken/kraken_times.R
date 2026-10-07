library(ggplot2)


flag_redo <- FALSE

###########
########### get times -- standard-8 
###########

num_cells <- 2742 #from hist file
for_cells <- 10000

if(flag_redo) {
  time_db <- system.time(system(
    "kraken2 --db /data/henlab/kraken/standard-8 --thread 10 --output out --gzip-compressed --paired empty.R1.fq.gz empty.R2.fq.gz"
  ))[3]
} else {
  time_db <- 42.118 
}

if(flag_redo) {
  time_db_reads <- system.time(system(
    "kraken2 --db /data/henlab/kraken/standard-8 --thread 10 --output out --gzip-compressed --paired fastp.1.R1.fq.gz fastp.1.R2.fq.gz"
  ))[3] 
} else {
  time_db_reads <- 160.794 
}


time_reads <- time_db_reads - time_db

#num_reads <- system("zcat fastp.1.R1.fq.gz | wc -l", intern=TRUE)
#39 844 353  #40M reads


time_naive_kraken  <- time_db * for_cells + time_reads/num_cells * for_cells
time_bascet_kraken <- time_db + time_reads/num_cells * for_cells

### How much slower is naive Kraken vs Bascet Kraken?
time_naive_kraken/time_bascet_kraken
#887x




ktimes <- data.frame(
  method=c(
    "Standard db, Naive KRAKEN2",
    "Standard db, Bascet KRAKEN2"
  ),
  time=c(time_naive_kraken, time_db_reads)
)

#method,time
#Naive KRAKEN2,standard-8,10967.16
#Bascet KRAKEN2,standard-8,38.801

###########
########### Plotting
###########

# standard-8 db: 7.5 gb

ktimes <- read.table("kraken_times.txt",sep=",", header = TRUE)
ktimes$method <- paste(ktimes$db, ktimes$method)
ktimes <- ktimes[ktimes$db=="Standard-8 db",]  ##### skip the big database
ktimes$method <- c("Naive KRAKEN2", "Bascet KRAKEN2") 
ktimes$time_h <- ktimes$time/3600

ktimes

ggplot(ktimes, aes(method, time/3600)) + 
  geom_bar(stat="identity") +
  theme_bw() + ylab("Time (hours), for 10k cells") + xlab("Method") +
  coord_flip()
ggsave("kraken_days_standard8.svg", width = 5, height=1.5)


#write.csv(ktimes, "/data/henriksson/kraken_for_julian.txt")


