library(ggplot2)

# Usage: Rscript supnote_speed_performance_shard.R <timings.tsv|run_directory> [plot_directory]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || length(args) > 2) {
  stop("Usage: Rscript supnote_speed_performance_shard.R <timings.tsv|run_directory> [plot_directory]")
}
timing_file <- if (dir.exists(args[1])) file.path(args[1], "timings.tsv") else args[1]
plot_dir <- if (length(args) == 2) args[2] else file.path(dirname(timing_file), "plots")
df <- read.delim(timing_file, check.names = FALSE)
required <- c("cells", "read_pairs", "repeat", "step", "wall_seconds", "exit_code")
if (!all(required %in% names(df))) stop("Timing table is missing required columns")
if (any(!df$step %in% c("debarcode", "shardify"))) stop("Unknown benchmark step")
if (any(!is.finite(df$exit_code))) stop("Invalid exit code in timing table")
if (any(df$exit_code != 0)) warning("Excluding failed commands from the plot")
df <- df[df$exit_code == 0, , drop = FALSE]
if (!nrow(df)) stop("No successful timings to plot")
if (any(!is.finite(df$cells) | df$cells < 1 | df$cells != floor(df$cells)) ||
    any(!is.finite(df$wall_seconds) | df$wall_seconds < 0)) {
  stop("Invalid cell counts or wall times")
}
if (anyDuplicated(df[c("cells", "repeat", "step")])) stop("Duplicate repetition/step rows")

summary_df <- do.call(rbind, lapply(split(df, list(df$cells, df$step), drop = TRUE), function(d) {
  data.frame(cells = d$cells[1], step = d$step[1], repetitions = nrow(d),
             median_seconds = median(d$wall_seconds),
             min_seconds = min(d$wall_seconds), max_seconds = max(d$wall_seconds))
}))
summary_df <- summary_df[order(summary_df$step, summary_df$cells), ]
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
write.table(summary_df, file.path(plot_dir, "performance_shard_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

step_levels <- c("debarcode", "shardify")
step_labels <- c("Debarcoding and sorting", "Sharding (10 shards)")
df$step <- factor(df$step, levels = step_levels, labels = step_labels)
summary_df$step <- factor(summary_df$step, levels = step_levels, labels = step_labels)
# A one-cell-count smoke run has points and ranges but no connecting line.
line_df <- summary_df[ave(summary_df$cells, summary_df$step, FUN = length) > 1, ]
cell_breaks <- sort(unique(df$cells))
cell_labels <- format(cell_breaks, big.mark = ",", scientific = FALSE, trim = TRUE)
df$wall_hours <- df$wall_seconds / 3600
summary_df$median_hours <- summary_df$median_seconds / 3600
summary_df$min_hours <- summary_df$min_seconds / 3600
summary_df$max_hours <- summary_df$max_seconds / 3600
line_df$median_hours <- line_df$median_seconds / 3600

p <- ggplot(df, aes(cells, wall_hours, color = step)) +
  geom_point(alpha = 0.5, size = 1.8) +
  geom_linerange(data = summary_df,
                 aes(x = cells, ymin = min_hours, ymax = max_hours, color = step),
                 inherit.aes = FALSE, linewidth = 0.6) +
  geom_line(data = line_df, aes(y = median_hours, group = step), linewidth = 0.8) +
  geom_point(data = summary_df, aes(y = median_hours), shape = 18, size = 2.8) +
  facet_wrap(~step, ncol = 1, scales = "free_y", drop = TRUE) +
  scale_color_manual(values = c("#0072B2", "#D55E00"), drop = FALSE) +
  labs(y = "Wall time (hours)", title = "Bascet runtime",
       caption = "Points: repetitions. Diamonds and lines: median. Ranges: min-max. Separate y scales.") +
  theme_bw() +
  theme(legend.position = "none", panel.grid.minor = element_blank())

# Keep the existing log-x view and add the requested log-log comparison.
linear_y <- scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08)))
plots <- list(
  linear = p + linear_y + scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)) +
    xlab("Number of cells"),
  logx = p + linear_y + scale_x_log10(breaks = cell_breaks, labels = cell_labels) +
    xlab("Number of cells (log scale)")
)
if(all(df$wall_seconds > 0)) {
  log_breaks <- 10^(0:floor(log10(max(df$cells))))
  log_breaks <- log_breaks[log_breaks >= min(df$cells)]
  plots$loglog <- p +
    scale_x_log10(breaks = log_breaks,
                  labels = format(log_breaks, big.mark = ",", scientific = FALSE, trim = TRUE)) +
    scale_y_log10() +
    labs(x = "Number of cells (log scale)", y = "Wall time (hours, log scale)")
} else {
  warning("Skipping log-log plot: some recorded durations are zero")
}
for (name in names(plots)) {
  prefix <- file.path(plot_dir, paste0("plot_perf_bascet_", name))
  ggsave(paste0(prefix, ".pdf"), plots[[name]], width = 7, height = 6)
  ggsave(paste0(prefix, ".png"), plots[[name]], width = 7, height = 6, dpi = 300)
}
message("Saved runtime plots and summary to ", plot_dir)
