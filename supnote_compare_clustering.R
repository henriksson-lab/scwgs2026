library(ggplot2)

datasets <- c("simulated5", "v8_wgs_mock5_pta")

methods <- c("alignment", "kraken", "infokmer", "countsketch")

grid_methods <- c("alignment", "kraken", "infokmer", "countsketch")

method_names <- c(
  alignment = "Alignment",
  kraken = "Kraken",
  infokmer = "InfoKmer",
  countsketch = "CountSketch"
)

dataset_names <- c(
  simulated5 = "Simulated",
  v8_wgs_mock5_pta = "Mock PTA"
)

root_dir <- "supnote_compare_clustering_sim"
plot_dir <- file.path("plots", "all")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

format_res <- function(x) {
  format(x, trim = TRUE, scientific = FALSE)
}

strip_panel_prefix <- function(x) {
  sub("^.*___", "", x)
}

read_summary <- function() {
  rows <- lapply(datasets, function(dataset) {
    path <- file.path(root_dir, dataset, "compare_clustering_sim_summary.csv")
    if (!file.exists(path)) {
      stop("Missing summary file: ", path)
    }
    df <- read.csv(path, check.names = FALSE)
    if (!"dataset" %in% colnames(df)) {
      df <- cbind(dataset = dataset, df)
    }
    df
  })
  all_cols <- unique(unlist(lapply(rows, colnames)))
  rows <- lapply(rows, function(df) {
    missing_cols <- setdiff(all_cols, colnames(df))
    for (col in missing_cols) {
      df[[col]] <- NA
    }
    df[, all_cols, drop = FALSE]
  })
  summary_df <- do.call(rbind, rows)
  summary_df$missing_alignment_species[is.na(summary_df$missing_alignment_species)] <- ""
  summary_df$dataset_label <- dataset_names[summary_df$dataset]
  summary_df$method_label <- method_names[summary_df$method]
  summary_df$param_label <- paste0(
    "k=", summary_df$selected_k_param,
    "\nres=", format_res(summary_df$selected_resolution),
    "\nclusters=", summary_df$n_clusters
  )
  summary_df$dataset_label <- factor(summary_df$dataset_label, levels = dataset_names[datasets])
  summary_df$method <- factor(summary_df$method, levels = methods)
  summary_df$method_label <- factor(summary_df$method_label, levels = method_names[methods])
  summary_df
}

read_confusion <- function(summary_df) {
  rows <- list()
  for (dataset in datasets) {
    for (method in methods) {
      path <- file.path(root_dir, dataset, paste0("compare_clustering_", method, "_confusion.csv"))
      if (!file.exists(path)) {
        stop("Missing confusion file: ", path)
      }

      df <- read.csv(path, check.names = FALSE)
      cluster_cols <- setdiff(colnames(df), c("method", "resolution", "species_aln_short"))
      long <- do.call(rbind, lapply(cluster_cols, function(cluster) {
        data.frame(
          dataset = dataset,
          method = method,
          species_aln_short = df$species_aln_short,
          cluster = cluster,
          n = df[[cluster]],
          check.names = FALSE
        )
      }))
      rows[[paste(dataset, method, sep = "_")]] <- long
    }
  }

  confusion <- do.call(rbind, rows)
  confusion <- merge(
    confusion,
    summary_df[, c("dataset", "method", "dataset_label", "method_label")],
    by = c("dataset", "method"),
    all.x = TRUE
  )
  confusion$method <- factor(confusion$method, levels = methods)
  confusion$dataset_label <- factor(confusion$dataset_label, levels = dataset_names[datasets])
  confusion$method_label <- factor(confusion$method_label, levels = method_names[methods])

  panel_orders <- lapply(split(confusion, list(confusion$dataset, confusion$method), drop = TRUE), function(panel_df) {
    dataset <- panel_df$dataset[1]
    method <- panel_df$method[1]

    mat <- xtabs(n ~ species_aln_short + cluster, data = panel_df)
    cluster_sizes <- colSums(mat)
    cluster_order <- names(sort(cluster_sizes, decreasing = TRUE))

    dominant_cluster <- apply(mat, 1, function(x) names(which.max(x))[1])
    dominant_n <- apply(mat, 1, max)
    species_sizes <- rowSums(mat)
    dominant_rank <- match(dominant_cluster, cluster_order)
    species_order <- names(sort(
      dominant_rank * 1e12 - dominant_n * 1e6 - species_sizes,
      decreasing = FALSE
    ))

    list(
      species_panel = paste(dataset, method, rev(species_order), sep = "___"),
      cluster_panel = paste(dataset, method, cluster_order, sep = "___")
    )
  })
  species_levels <- unlist(lapply(panel_orders, function(x) x$species_panel), use.names = FALSE)
  cluster_levels <- unlist(lapply(panel_orders, function(x) x$cluster_panel), use.names = FALSE)

  confusion$species_panel <- paste(
    confusion$dataset,
    confusion$method,
    confusion$species_aln_short,
    sep = "___"
  )
  confusion$cluster_panel <- paste(
    confusion$dataset,
    confusion$method,
    confusion$cluster,
    sep = "___"
  )
  confusion$species_panel <- factor(confusion$species_panel, levels = species_levels)
  confusion$cluster_panel <- factor(confusion$cluster_panel, levels = cluster_levels)
  confusion
}

read_grid_scan <- function(summary_df) {
  rows <- list()
  for (dataset in datasets) {
    for (method in grid_methods) {
      path <- file.path(root_dir, dataset, paste0("compare_clustering_", method, "_resolution_scan.csv"))
      if (!file.exists(path)) {
        stop("Missing resolution scan file: ", path)
      }

      df <- read.csv(path, check.names = FALSE)
      df$dataset <- dataset
      df$method <- method
      rows[[paste(dataset, method, sep = "_")]] <- df
    }
  }

  all_cols <- unique(unlist(lapply(rows, colnames)))
  rows <- lapply(rows, function(df) {
    missing_cols <- setdiff(all_cols, colnames(df))
    for (col in missing_cols) {
      df[[col]] <- NA
    }
    df[, all_cols, drop = FALSE]
  })
  scan <- do.call(rbind, rows)
  scan$resolution_label <- format_res(scan$resolution)
  scan <- merge(
    scan,
    summary_df[, c("dataset", "method", "dataset_label", "method_label", "selected_k_param", "selected_resolution")],
    by = c("dataset", "method"),
    all.x = TRUE
  )
  scan$selected <- scan$k.param == scan$selected_k_param & scan$resolution == scan$selected_resolution
  scan$method <- factor(scan$method, levels = grid_methods)
  scan$dataset_label <- factor(scan$dataset_label, levels = dataset_names[datasets])
  scan$method_label <- factor(scan$method_label, levels = method_names[grid_methods])
  scan$k_label <- factor(scan$k.param, levels = sort(unique(scan$k.param)))
  scan$resolution_label <- factor(
    scan$resolution_label,
    levels = format_res(sort(unique(scan$resolution)))
  )
  scan
}

summary_df <- read_summary()
panel_labels <- unique(summary_df[, c("dataset_label", "method_label", "param_label")])
grid_panel_labels <- unique(summary_df[summary_df$method %in% grid_methods, c("dataset_label", "method_label", "param_label")])
write.csv(
  summary_df[, intersect(c(
    "dataset",
    "method",
    "n_cells",
    "selected_k_param",
    "selected_resolution",
    "n_clusters",
    "overall_purity",
    "largest_cluster_fraction",
    "min_cluster_size",
    "n_species_ground_truth",
    "missing_alignment_species"
  ), colnames(summary_df))],
  file.path(root_dir, "compare_clustering_sim_summary_all.csv"),
  row.names = FALSE,
  na = ""
)

confusion_df <- read_confusion(summary_df)
confusion_plot <- ggplot(confusion_df, aes(cluster_panel, species_panel, fill = n)) +
  geom_tile(color = "white", linewidth = 0.15) +
  geom_label(
    data = panel_labels,
    aes(x = Inf, y = Inf, label = param_label),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.1,
    size = 2.3,
    linewidth = 0.15,
    fill = "white"
  ) +
  scale_fill_gradient(low = "white", high = "#2166ac", name = "Cells") +
  scale_x_discrete(labels = strip_panel_prefix) +
  scale_y_discrete(labels = strip_panel_prefix) +
  facet_wrap(dataset_label ~ method_label, ncol = 4, scales = "free") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 9),
    panel.grid = element_blank(),
    strip.text = element_text(size = 8),
    legend.position = "right"
  ) +
  xlab("Seurat cluster") +
  ylab("Alignment species ground truth")

ggsave(
  file.path(plot_dir, "compare_clustering_confusion_3x2.svg"),
  plot = confusion_plot,
  width = 13,
  height = 6,
  limitsize = FALSE
)

ggsave(
  file.path(plot_dir, "compare_clustering_confusion_4x2.svg"),
  plot = confusion_plot,
  width = 13,
  height = 6,
  limitsize = FALSE
)

scan_df <- read_grid_scan(summary_df)
grid_scan_plot <- ggplot(scan_df, aes(resolution_label, k_label, fill = n_clusters)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = n_clusters), size = 2.2) +
  geom_label(
    data = grid_panel_labels,
    aes(x = Inf, y = Inf, label = param_label),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.1,
    size = 2.3,
    linewidth = 0.15,
    fill = "white"
  ) +
  geom_point(
    data = scan_df[scan_df$selected, , drop = FALSE],
    shape = 21,
    size = 4,
    stroke = 1,
    fill = NA,
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "#b2182b",
    mid = "white",
    high = "#2166ac",
    midpoint = 10,
    name = "Clusters"
  ) +
  facet_grid(dataset_label ~ method_label, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 9),
    panel.grid = element_blank(),
    strip.text = element_text(size = 8),
    legend.position = "right"
  ) +
  xlab("Resolution") +
  ylab("k")

ggsave(
  file.path(plot_dir, "compare_clustering_grid_scan_3x2.svg"),
  plot = grid_scan_plot,
  width = 13,
  height = 5,
  limitsize = FALSE
)

ggsave(
  file.path(plot_dir, "compare_clustering_grid_scan_4x2.svg"),
  plot = grid_scan_plot,
  width = 13,
  height = 5,
  limitsize = FALSE
)
