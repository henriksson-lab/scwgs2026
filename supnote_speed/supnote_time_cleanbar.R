 
plotDirAll <- file.path("~/github/scmetag-revision/plots/all")


numreads <- 9371176 

time_bascet <- 128.65
time_cleanbar <- 1187.32

time_cleanbar/time_bascet #9.2x

df <- data.frame(
  time=c(time_bascet, time_cleanbar),
  tool=c("Bascet","CleanBar")
)

ggplot(df, aes(tool,time)) + 
  geom_bar(stat="identity", fill="darkblue") +
  theme_bw() +
  xlab("Software") + ylab("Wall time (s)")+
  coord_flip()
ggsave(filename = file.path(plotDirAll,"time_cleanbar.pdf"), width = 4, height = 1)


#how long for cleanbar to debarcode 800M reads => 28h
time_cleanbar/9371176*800e6/3600
