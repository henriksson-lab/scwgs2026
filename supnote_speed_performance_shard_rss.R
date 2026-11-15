library(ggplot2)

# Usage: Rscript supnote_speed_performance_shard_rss.R <run_directory> [plot_directory]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || length(args) > 2) {
  stop("Usage: Rscript supnote_speed_performance_shard_rss.R <run_directory> [plot_directory]")
}
run_dir <- args[1]
plot_dir <- if(length(args) == 2) args[2] else file.path(run_dir, "plots")
rss <- read.delim(file.path(run_dir, "rss.tsv"), check.names = FALSE)
timings <- read.delim(file.path(run_dir, "timings.tsv"), check.names = FALSE)
stopifnot(all(c("cells", "repeat", "step", "elapsed_seconds", "pid", "rss_bytes") %in% names(rss)),
          all(c("cells", "repeat", "step", "peak_rss_bytes", "exit_code") %in% names(timings)))
timings <- timings[timings$exit_code == 0 & is.finite(timings$peak_rss_bytes), , drop = FALSE]
if (!nrow(timings)) stop("No successful memory measurements")
if (any(!is.finite(rss$elapsed_seconds) | rss$elapsed_seconds < 0 | !is.finite(rss$rss_bytes) | rss$rss_bytes < 0)) {
  stop("Invalid RSS measurements")
}
counts <- sort(unique(timings$cells))
step_names <- c("debarcode", "shardify")
step_labels <- c("Debarcoding and sorting", "Sharding (10 shards)")
count_labels <- paste(format(counts, big.mark = ",", scientific = FALSE, trim = TRUE), ifelse(counts == 1, "cell", "cells"))
label_data <- function(d) {
  d$step <- factor(d$step, levels = step_names, labels = step_labels)
  d$cell_label <- factor(d$cells, levels = counts, labels = count_labels)
  d$repetition <- factor(d[["repeat"]])
  d
}
rss <- label_data(rss)
timings <- label_data(timings)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

rss_plot <- droplevels(rss[rss$cells >= 10000, , drop = FALSE])
trace_timings <- droplevels(timings[timings$cells >= 10000, , drop = FALSE])
trace_counts <- sort(unique(trace_timings$cells))
if(length(trace_counts)) {
present <- unique(rss_plot[c("step", "cell_label")])
panels <- unique(trace_timings[c("step", "cell_label")])
panel_key <- function(d) paste(d$step, d$cell_label)
missing <- panels[!panel_key(panels) %in% panel_key(present), , drop = FALSE]
missing$x_min <- missing$y_min <- rep(0, nrow(missing))
missing$x_max <- missing$y_max <- missing$x <- missing$y <- rep(0, nrow(missing))
missing$label <- rep("", nrow(missing))
for(i in seq_len(nrow(missing))) {
  d <- trace_timings[trace_timings$step == missing$step[i] & trace_timings$cell_label == missing$cell_label[i], ]
  duration <- if("wrapper_wall_seconds" %in% names(d)) d$wrapper_wall_seconds else d$wall_seconds
  missing$x_max[i] <- max(duration)
  missing$y_max[i] <- max(d$peak_rss_bytes) / 1024^3 * 1.1
  missing$x[i] <- missing$x_max[i] / 2
  missing$y[i] <- missing$y_max[i] / 2
  missing$label[i] <- sprintf("No periodic samples\nPeak RSS: %.3f-%.3f GiB",
                             min(d$peak_rss_bytes) / 1024^3, max(d$peak_rss_bytes) / 1024^3)
}
interval <- ""
if(requireNamespace("jsonlite", quietly = TRUE) && file.exists(file.path(run_dir, "metadata.json"))) {
  metadata <- jsonlite::fromJSON(file.path(run_dir, "metadata.json"))
  if(!is.null(metadata$rss_interval_seconds)) interval <- paste0("Requested interval: ", metadata$rss_interval_seconds, " s. ")
}
p <- ggplot(rss_plot, aes(elapsed_seconds, rss_bytes / 1024^3, color = repetition,
                     group = interaction(repetition, pid))) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_blank(data = missing, aes(x = x_min, y = y_min), inherit.aes = FALSE) +
  geom_blank(data = missing, aes(x = x_max, y = y_max), inherit.aes = FALSE) +
  geom_text(data = missing, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3) +
  facet_wrap(~step + cell_label, ncol = length(trace_counts), scales = "free", drop = FALSE) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.04))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Elapsed time since step submission (seconds)", y = "Bascet RSS (GiB)",
       color = "Repetition", title = "Bascet resident memory over time",
       caption = paste0(interval, "Only the Bascet process is sampled; R workers are excluded.\nPanels use separate axes. Peak RSS is measured independently by GNU time.")) +
  theme_bw() +
  theme(legend.position = "top", panel.grid.minor = element_blank())
width <- max(7, 2.7 * length(trace_counts))
ggsave(file.path(plot_dir, "plot_perf_bascet_rss.pdf"), p, width = width, height = 6)
ggsave(file.path(plot_dir, "plot_perf_bascet_rss.png"), p, width = width, height = 6, dpi = 200)
} else {
  message("No inputs with at least 10,000 cells; skipping the RSS-over-time plot")
}

peaks <- do.call(rbind, lapply(split(timings, list(timings$cells, timings$step), drop = TRUE), function(d) {
  data.frame(cells = d$cells[1], step = d$step[1], median_peak = median(d$peak_rss_bytes))
}))
lines <- peaks[ave(peaks$cells, peaks$step, FUN = length) > 1, ]
p_peak <- ggplot(timings, aes(cells, peak_rss_bytes / 1024^3, color = step)) +
  geom_point(alpha = 0.6) +
  geom_line(data = lines, aes(y = median_peak / 1024^3, group = step), linewidth = 0.8) +
  facet_wrap(~step, scales = "free_y", ncol = 1) +
  scale_x_log10(breaks = counts, labels = format(counts, big.mark = ",", scientific = FALSE, trim = TRUE)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Number of cells (log scale)", y = "Peak Bascet RSS (GiB)",
       title = "Peak resident memory", caption = "GNU time peak RSS. Points: repetitions. Lines: medians. Separate y scales.") +
  theme_bw() + theme(legend.position = "none", panel.grid.minor = element_blank())
ggsave(file.path(plot_dir, "plot_perf_bascet_peak_rss.pdf"), p_peak, width = 7, height = 6)
ggsave(file.path(plot_dir, "plot_perf_bascet_peak_rss.png"), p_peak, width = 7, height = 6, dpi = 300)
message("Saved RSS traces and peak RSS plots to ", plot_dir)
