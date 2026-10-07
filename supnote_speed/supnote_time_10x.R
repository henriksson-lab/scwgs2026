#walltime_seconds=8515.30
#walltime_hms=2:21:55

library(ggplot2)
library(khroma)



log_file <- "time_10x_bascet.txt"
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
  Zdebarcode          = "Debarcode",
  Zshardify           = "Shardify",
  Zaln                = "Align",
  Zcountf             = "Count Features"
)

dat$activity <- ifelse(
  dat$activity %in% names(activity_map),
  activity_map[dat$activity],
  dat$activity
)

# Plot order follows the known workflow, then appends any future categories.
activity_levels <- c(unname(activity_map), setdiff(unique(dat$activity), unname(activity_map)))
cellranger_time_hours <- 8515.30/3600

#bright <- color("bright")
fill_scale <- scale_fill_manual(values = c(setNames(bright(length(activity_levels)),
                                           activity_levels),
                                      CellRanger = "grey65"))

# Build the bar with labels floating above each segment at 45°
bright <- color("muted")
plot_time <- function(dat_agg, ylab_text, show_text = TRUE) {
  ordered <- rev(activity_levels)
  d <- dat_agg[match(ordered, dat_agg$activity), ]
  d <- d[!is.na(d$activity), ]
  d$activity <- factor(d$activity, levels = ordered)
  d$workflow <- "BASCET"
  # position_stack() puts the LAST factor level at y=0 and stacks upward,
  # so midpoints must be computed on the reversed duration vector.
  dur_rev <- rev(d$duration)
  ymid_rev <- cumsum(dur_rev) - dur_rev / 2
  d$ymid_h <- rev(ymid_rev) / 3600

  d_compare <- data.frame(
    activity = factor("CellRanger", levels = c(ordered, "CellRanger")),
    workflow = "CellRanger",
    duration = cellranger_time_hours * 3600,
    ymid_h = cellranger_time_hours / 2,
    stringsAsFactors = FALSE
  )
  d$activity <- factor(d$activity, levels = c(ordered, "CellRanger"))
  d_plot <- rbind(d, d_compare)
  d_plot$workflow <- factor(d_plot$workflow, levels = c("CellRanger", "BASCET"))

  p <- ggplot(d_plot, aes(workflow, duration/3600, fill = activity)) +
    geom_bar(stat = "identity", position = position_stack(), width = 0.4) +
    #fill_scale + #########################################################################
    scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0)) +
    coord_flip(clip = "off") +
    theme_bw() +
    theme(legend.position = "none",
          axis.ticks.y    = element_blank(),
          axis.title.y    = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.margin = margin(60, 20, 5, 5)) +
    ylab(ylab_text) + xlab("")

  if (show_text) {
    p <- p +
      geom_text(data = d, aes(x = "BASCET", y = ymid_h, label = activity),
                angle = 45, hjust = 0, vjust = 0)
  }

  p
}

###########
########### Wall time
###########

dat_max <- sqldf::sqldf("select activity, median(duration) as duration from dat group by activity")
p1 <- plot_time(dat_max, "Wall time (Hours)", show_text = FALSE)
p1

time_bascet <- sum(dat_max$duration)/3600
time_CR <- cellranger_time_hours
time_CR/time_bascet


###########
########### CPU time
###########

dat_sum <- sqldf::sqldf("select activity, sum(duration) as duration from dat group by activity")
p2 <- plot_time(dat_sum, "CPU time (Hours)", show_text = TRUE)

ptot <- egg::ggarrange(p2,p1, ncol=1)
ggsave(plot=ptot, filename = "plots/time_end2end.svg", width = 5, height = 4)

#/data/henriksson/github/scmetag-revision/supnote_speed/plots/time_end2end.svg




