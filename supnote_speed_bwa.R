 

numcell <- 10e3
startup <- 11.1
time_per_cell <- 0.61


total_time_naive <- (startup+time_per_cell)*numcell
total_time_batch <- startup+time_per_cell*numcell


df <- data.frame(
  samp=c("Naive","Batched"),
  time=c(total_time_naive,total_time_batch)
)

ggplot(df, aes(samp,time/3600)) + 
  geom_bar(stat="identity") + 
  coord_flip()+
  xlab("")+
  ylab("Time(h)") +
  theme_bw()


plotDirAll <- file.path("~/github/scmetag-revision/plots/all")
ggsave(filename = file.path(plotDirAll,"plot_perf_bwamem2.pdf"), width = 6, height = 4)