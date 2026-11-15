library(ggplot2)

# Same timing assumptions as supnote_speed_bwa.R; keep the old plot intact.
startup <- 11.1 # seconds
time_per_cell <- 0.61 # seconds
num_cells <- unique(round(seq(1, 1000000, length.out = 10001)))

# Cell-by-cell processing pays startup for every cell; batching pays it once.
df <- rbind(
  data.frame(cells = num_cells, method = "Naive",
             seconds = num_cells * (startup + time_per_cell)),
  data.frame(cells = num_cells, method = "Batched",
             seconds = startup + num_cells * time_per_cell)
)
df$method <- factor(df$method, levels = c("Naive", "Batched"))

p <- ggplot(df, aes(cells, seconds / 3600, color = method)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("#666666", "#0072B2")) +
  scale_x_continuous(breaks = c(1, 250000, 500000, 750000, 1000000),
                     labels = c("1", "250,000", "500,000", "750,000", "1,000,000")) +
  labs(x = "Number of cells", y = "Estimated total runtime (hours)",
       color = NULL, title = "BWA-MEM2") +
  theme_bw() +
  theme(legend.position = "top", panel.grid.minor = element_blank())

plot_dir <- file.path("plots", "all")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(plot_dir, "plot_perf_bwamem2_scaling.pdf"), p, width = 6, height = 4)
ggsave(file.path(plot_dir, "plot_perf_bwamem2_scaling.png"), p, width = 6, height = 4, dpi = 300)
