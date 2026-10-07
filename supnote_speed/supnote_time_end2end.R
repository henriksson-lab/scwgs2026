library(ggplot2)
library(khroma)


log_file <- "timelogs/v9_wgs_mock5_pta_BENCH/logsummary.txt"
#log_file <- "alltimes.v9_wgs_mock5_pta.log"
raw <- readLines(log_file)


if(FALSE) {
  ###### For running on beagle with local runner

  # Keep only the per-shard totals (each invocation ends with "Success! elapsed=Xs")
  success_lines <- grep("Success!", raw, fixed = TRUE, value = TRUE)
  
  pattern <- "^logs/([^.]+)\\.([^.]+)\\.log:.*elapsed=([0-9.]+)s"
  m <- regmatches(success_lines, regexec(pattern, success_lines))
  keep <- vapply(m, length, integer(1)) == 4L
  m <- m[keep]
  
  m <- success_lines
  
  dat <- data.frame(
    activity = vapply(m, `[`, character(1), 2),
    shard    = vapply(m, `[`, character(1), 3),
    duration = as.numeric(vapply(m, `[`, character(1), 4)),
    stringsAsFactors = FALSE
  )
  

  ###### Rename activities
  activity_map <- c(
    Zgetraw             = "Debarcode",
    Ztransform          = "Shardify",
    Zfastqc             = "FastQC",
    Zaln                = "Align (BWA)",
    ZfilterAln          = "Filter Align",
    Zcountchrom         = "CountChrom",
    Z_kraken_fq         = "KRAKEN2",
    Zskesa              = "SKESA"
  )
  dat <- dat[dat$activity %in% names(activity_map), ]
  dat$activity <- activity_map[dat$activity]
  
  
} else {
  #For running on HPC2N with slurm runner
  
  s_proc <- raw[seq(from=1, to=length(raw), by=2)]
  s_elapsed <- raw[seq(from=2, to=length(raw), by=2)]
  s_proc_name <- str_split_i(s_proc,"=",3)
  s_time <- str_split_i(s_elapsed,"=",2)
  s_time_sec <- as.double(str_remove(s_time,"s"))

  s_proc_start <- str_split_i(s_proc," ",1)
  
  dat <- data.frame(
    starttime=s_proc_start,
    activity=s_proc_name,
    duration=s_time_sec
  )
  
  dat <- dat[order(dat$starttime),] 
  
#  s_time_sec <- as.double(str_remove(s_time,"s"))
#  s_proc
  
  unique(dat$activity)
  sort(unique(dat$activity))
  
  ###### Rename activities
  activity_map <- c(
    Debarcode      = "Debarcode",
    Shardify       = "Shardify",
    #Zfastqc             = "FastQC", #### we did not run here
    Align          = "Align (BWA-MEM2)",  ### note .. there are two alignment steps! how to handle? guess we can report the latter for simplicity
    Filterbam      = "Filter BAM",
    Countchrom     = "CountChrom",
    Countsketch    = "Countsketch",
    Kraken         = "Kraken2",
    MinhashFq      = "MinhashFq",
    MinhashHist    = "MinhashHist",
    DetectKmerFq   = "DetectKmerFq",
    Skesa          = "SKESA"
  )
  length(activity_map)==length(unique(dat$activity))  #check that all activities are mapped
  
  dat <- dat[dat$activity %in% names(activity_map), ]
  dat$activity <- activity_map[dat$activity]
  
}


# Plot order follows the order of activity_map
activity_levels <- unname(activity_map)

#bright <- color("bright")
#bright <- color("muted")
bright <- khroma::color("sunset")(11)
fill_scale <- scale_fill_manual(values = setNames(bright,
                                                  activity_levels))

# Build the bar with labels floating above each segment at 45°
plot_time <- function(dat_agg, ylab_text) {
  ordered <- rev(activity_levels)
  d <- dat_agg[match(ordered, dat_agg$activity), ]
  d <- d[!is.na(d$activity), ]
  d$activity <- factor(d$activity, levels = ordered)
  # position_stack() puts the LAST factor level at y=0 and stacks upward,
  # so midpoints must be computed on the reversed duration vector.
  dur_rev <- rev(d$duration)
  ymid_rev <- cumsum(dur_rev) - dur_rev / 2
  d$ymid_h <- rev(ymid_rev) / 3600

  ggplot(d, aes(1, duration/3600, fill = activity)) +
    geom_bar(stat = "identity", position = position_stack(), width = 0.4) +
    geom_text(aes(x = 1.25, y = ymid_h, label = activity),
              angle = 45, hjust = 0, vjust = 0) +
    fill_scale +
    coord_flip(clip = "off") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.y     = element_blank(),
          axis.ticks.y    = element_blank(),
          axis.title.y    = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.margin = margin(60, 20, 5, 5)) +
    ylab(ylab_text) + xlab("") +
    expand_limits(x = c(0.5, 2))
}

###########
########### Wall time
###########

fmt_hm <- function(x) {
  out <- rep(NA_character_, length(x))
  ok <- !is.na(x)
  total_min <- round(x[ok] * 60)
  h <- total_min %/% 60
  m <- total_min %% 60
  vals <- character(length(h))
  vals[h == 0]         <- sprintf("%dmin", m[h == 0])
  vals[h > 0 & m == 0] <- sprintf("%dh",   h[h > 0 & m == 0])
  vals[h > 0 & m != 0] <- sprintf("%dh%dmin", h[h > 0 & m != 0], m[h > 0 & m != 0])
  out[ok] <- vals
  out
}

dat_max <- sqldf::sqldf("select activity, median(duration) as duration from dat group by activity")
p1 <- plot_time(dat_max, "Wall time") +
  scale_y_continuous(breaks = seq(0, 24, by = 1/3), labels = fmt_hm)
p1

###########
########### CPU time
###########


dat_sum <- sqldf::sqldf("select activity, sum(duration) as duration from dat group by activity")
p2 <- plot_time(dat_sum, "CPU time (Hours)")


plotDirAll <- file.path("~/github/scmetag-revision/plots/all")

ptot <- egg::ggarrange(p1,p2, ncol=1)
ggsave(plot=ptot, filename = file.path(plotDirAll,"supnote_time_end2end.svg"), width = 5, height = 4)





###########
########### Time for the different UMAP methods
###########

annotate_per_umap <- function(dat_max) {
  dat_max$umap_type <- ""
  dat_max$umap_type[dat_max$activity=="Align (BWA-MEM2)"] <- "Count per chromosome"
  dat_max$umap_type[dat_max$activity=="CountChrom"]       <- "Count per chromosome"
  
  dat_max$umap_type[dat_max$activity=="Countsketch"] <- "Countsketch"
  
  dat_max$umap_type[dat_max$activity=="DetectKmerFq"] <- "Informative KMER"
  dat_max$umap_type[dat_max$activity=="MinhashFq"]    <- "Informative KMER"
  
  dat_max$umap_type[dat_max$activity=="Kraken2"] <- "KMER database"
  
  dat_max <- dat_max[dat_max$umap_type != "",]
  dat_max
}




plot_time_umap <- function(dat_agg, ylab_text) {
  ordered <- rev(activity_levels)
  d <- dat_agg[match(ordered, dat_agg$activity), ]
  d <- d[!is.na(d$activity), ]
  d$activity <- factor(d$activity, levels = ordered)
  # position_stack() puts the LAST factor level at y=0 and stacks upward,
  # so midpoints must be computed on the reversed duration vector.
  dur_rev <- rev(d$duration)
  ymid_rev <- cumsum(dur_rev) - dur_rev / 2
  d$ymid_h <- rev(ymid_rev) / 60
  
  ggplot(d, aes(umap_type, duration/60, fill = activity)) +
    geom_bar(stat = "identity", position = position_stack(), width = 0.4) +
#    geom_text(aes(x = 1.25, y = ymid_h, label = activity),
#              angle = 45, hjust = 0, vjust = 0) +
    fill_scale +
    coord_flip(clip = "off") +
    theme_bw() +
    theme(legend.position = "none",
          #axis.text.y     = element_blank(),
          axis.ticks.y    = element_blank(),
          axis.title.y    = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()
    #      plot.margin = margin(60, 20, 5, 5)
    ) +
    ylab(ylab_text) 
    #xlab("") #+
#    expand_limits(x = c(0.5, 2))
}


dat_per_umap <- annotate_per_umap(dat_max)
plot_time_umap(dat_per_umap, "Wall time (minutes)") 
ggsave(filename = file.path(plotDirAll,"supnote_time_umaps_walltime.svg"), width = 5, height = 2)



dat_per_umap <- annotate_per_umap(dat_sum)
plot_time_umap(dat_per_umap, "CPU time (minutes)") 
ggsave(filename = file.path(plotDirAll,"supnote_time_umaps_cputime.svg"), width = 5, height = 2)


