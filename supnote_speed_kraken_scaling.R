library(ggplot2)

# Standard-8 timing model from supnote_speed/benchmark_kraken/kraken_times.R.
# Use the explicit benchmark timings; the old bar plot reads a different
# set of saved values from kraken_times.txt. Keep the original code intact.
startup <- 42.118 # seconds, measured with empty input
benchmark_total <- 160.794 # seconds, including startup
benchmark_cells <- 2742
time_per_cell <- (benchmark_total - startup) / benchmark_cells
num_cells <- unique(round(seq(1, 1000000, length.out = 10001)))

# Cell-by-cell processing loads the database for every cell; batching loads
# it once. Both models assume the same processing time per cell.
df <- rbind(
  data.frame(cells = num_cells, method = "Naive Kraken2",
             seconds = num_cells * (startup + time_per_cell)),
  data.frame(cells = num_cells, method = "Bascet-Kraken2",
             seconds = startup + num_cells * time_per_cell)
)
df$method <- factor(df$method, levels = c("Naive Kraken2", "Bascet-Kraken2"))

p <- ggplot(df, aes(cells, seconds / 3600, color = method)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("#666666", "#0072B2")) +
  scale_x_continuous(breaks = c(1, 250000, 500000, 750000, 1000000),
                     labels = c("1", "250,000", "500,000", "750,000", "1,000,000")) +
  labs(x = "Number of cells", y = "Estimated total runtime (hours)",
       color = NULL, title = "Kraken2 (Standard-8 database)") +
  theme_bw() +
  theme(legend.position = "top", panel.grid.minor = element_blank())

plot_dir <- file.path("plots", "all")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(plot_dir, "plot_perf_kraken2_scaling.pdf"), p, width = 6, height = 4)
ggsave(file.path(plot_dir, "plot_perf_kraken2_scaling.png"), p, width = 6, height = 4, dpi = 300)
