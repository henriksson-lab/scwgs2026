
log_file <- "speed_zheng2022.log"
raw <- readLines(log_file)

# Keep only the per-shard totals (each invocation ends with "Success! elapsed=Xs")
success_lines <- grep("Success!", raw, fixed = TRUE, value = TRUE)

pattern <- "^(?:logs/)?([^.]+)\\.([^.]+)\\.log:.*elapsed=([0-9.]+)s"
m <- regmatches(success_lines, regexec(pattern, success_lines))
keep <- vapply(m, length, integer(1)) == 4L
m <- m[keep]

dat <- data.frame(
  activity = vapply(m, `[`, character(1), 2),
  shard    = vapply(m, `[`, character(1), 3),
  duration = as.numeric(vapply(m, `[`, character(1), 4)),
  stringsAsFactors = FALSE
)

activity_map <- c(
  ZimportSRA          = "Import SRA",
  Zgetraw             = "Debarcode",
  Ztransform          = "Shardify",
  Zfastqc             = "FastQC",
  Zaln                = "Align (BWA)",
  ZfilterAln          = "Filter Align",
  Zcountchrom         = "CountChrom",
  Z_kraken_fq         = "KRAKEN2",
  Zskesa              = "SKESA",
  Zcs                 = "Count Sketch"
#  Zmap__minhash_fq    = "MinHash Map",
#  bascet_minhash_hist = "MinHash Hist"
)
dat <- dat[dat$activity %in% names(activity_map), ]
activity_levels <- unname(activity_map[names(activity_map) %in% unique(dat$activity)])
dat$activity <- activity_map[dat$activity]

# Plot order follows the order of activity_map

library(ggplot2)
library(khroma)
#bright <- color("bright")
bright <- color("muted")
fill_scale <- scale_fill_manual(values = setNames(bright(length(activity_levels)),
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
    ylab(ylab_text) + xlab("")
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
  scale_y_continuous(labels = fmt_hm)

###########
########### CPU time
###########

#dat_sum <- sqldf::sqldf("select activity, sum(duration) as duration from dat group by activity")
p2 <- plot_time(dat_sum, "Wall time (Hours)")  #CPU time really, but here download dominated so more honest

#ptot <- egg::ggarrange(p1,p2, ncol=1)
p2
ggsave(plot=p2, filename = "plots/time_zheng2022.svg", width = 5, height = 2) 
