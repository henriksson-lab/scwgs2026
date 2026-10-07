library(Seurat)
library(ggplot2)
library(stringr)

dataset_configs <- list(
  simulated5 = list(
    dataset_name = "simulated5",
    enabled = FALSE,
    expected_n_species = 10,
    expected_species_map = "map_seq2strain.sim.csv"
  ),
  v8_wgs_mock5_pta = list(
    dataset_name = "v8_wgs_mock5_pta",
    enabled = TRUE,
    expected_n_species = 10,
    expected_species_map = "map_seq2strain.sim.csv"
  )
)


########## Helper functions
canonical_cell_names <- function(x) {
  stringr::str_remove(x, "^BASCET_")
}

load_expected_species <- function(map_path) {
  if (!file.exists(map_path)) {
    warning("Cannot find expected strain map: ", map_path)
    return(character(0))
  }

  mapSeq2strain <- read.csv(map_path)
  sort(unique(stringr::str_split_i(mapSeq2strain$strain, stringr::fixed(" ("), 1)))
}

read_cached_object <- function(bascetRoot, filename) {
  path <- file.path(bascetRoot, filename)
  if (!file.exists(path)) {
    stop("Missing cached object: ", path)
  }
  readRDS(path)
}

confusion_to_df <- function(tab, method, resolution) {
  df <- as.data.frame.matrix(tab)
  df <- data.frame(species_aln_short = rownames(df), df, check.names = FALSE)
  df$method <- method
  df$resolution <- resolution
  df[, c("method", "resolution", "species_aln_short", setdiff(colnames(df), c("method", "resolution", "species_aln_short")))]
}


########## 
plot_confusion <- function(tab, method, resolution) {
  df <- as.data.frame(tab)
  colnames(df) <- c("species_aln_short", "cluster", "n")
  ggplot(df, aes(cluster, species_aln_short, fill = n)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    ) +
    xlab("Seurat cluster") +
    ylab("Alignment species ground truth") +
    ggtitle(sprintf("%s, resolution %.4g", method, resolution))
}

add_alignment_truth <- function(adata, alignment_meta, method) {
  cell_index <- match(canonical_cell_names(colnames(adata)), rownames(alignment_meta))
  keep_cells <- !is.na(cell_index)
  if (sum(keep_cells) == 0) {
    stop("No overlapping cells between ", method, " and alignment ground truth")
  }

  adata <- adata[, keep_cells]
  alignment_index <- cell_index[keep_cells]
  adata$species_aln <- alignment_meta[alignment_index, "species_aln"]
  adata$species_aln_short <- alignment_meta[alignment_index, "species_aln_short"]
  adata <- adata[, !is.na(adata$species_aln_short)]
  adata
}

cluster_metrics <- function(tab) {
  n <- sum(tab)
  cluster_sizes <- colSums(tab)
  nonempty_cluster_sizes <- cluster_sizes[cluster_sizes > 0]
  dominant_n <- apply(tab, 2, max)

  data.frame(
    n_clusters = length(nonempty_cluster_sizes),
    overall_purity = sum(dominant_n) / n,
    largest_cluster_fraction = max(nonempty_cluster_sizes) / n,
    min_cluster_size = min(nonempty_cluster_sizes),
    max_cluster_size = max(nonempty_cluster_sizes),
    median_cluster_size = median(nonempty_cluster_sizes),
    n_singleton_clusters = sum(nonempty_cluster_sizes == 1)
  )
}

select_resolution <- function(
    adata,
    method,
    reduction,
    dims,
    expected_n_clusters,
    max_resolution = NULL,
    max_k_param = NULL,
    resolutions = c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2),
    k.params = c(10, 20, 50, 100, 200, 500, 1000),
    min_acceptable_purity = 0.7,
    max_acceptable_largest_cluster_fraction = 0.9,
    stop_row_at_n_clusters = {
      value <- Sys.getenv("STOP_ROW_AT_N_CLUSTERS", "")
      if (nzchar(value)) as.integer(value) else NA_integer_
    }) {

  if (nzchar(Sys.getenv("RESOLUTIONS", ""))) {
    resolutions <- as.numeric(trimws(strsplit(Sys.getenv("RESOLUTIONS"), ",", fixed = TRUE)[[1]]))
  }
  if (nzchar(Sys.getenv("K_PARAMS", ""))) {
    k.params <- as.integer(trimws(strsplit(Sys.getenv("K_PARAMS"), ",", fixed = TRUE)[[1]]))
  }
  if (!is.null(max_resolution)) {
    resolutions <- resolutions[resolutions <= max_resolution]
  }
  if (!is.null(max_k_param)) {
    k.params <- k.params[k.params <= max_k_param]
  }

  k.params <- k.params[k.params < ncol(adata)]
  scan_rows <- list()
  cluster_labels <- list()

  for (k.param in k.params) {
    message("  ", method, ": scanning k=", k.param)
    nn_graph <- paste0(method, "_nn_k", k.param)
    neighbor_graph <- paste0(method, "_snn_k", k.param)
    adata <- FindNeighbors(
      adata,
      reduction = reduction,
      dims = dims,
      k.param = k.param,
      graph.name = c(nn_graph, neighbor_graph),
      verbose = FALSE
    )

    for (res in resolutions) {
      message("    ", method, ": k=", k.param, ", res=", res)
      tmp <- FindClusters(
        adata,
        graph.name = neighbor_graph,
        resolution = res,
        algorithm = 1,
        verbose = FALSE
      )
      cluster_col <- paste0(neighbor_graph, "_res.", res)
      if (!cluster_col %in% colnames(tmp@meta.data)) {
        cluster_col <- "seurat_clusters"
      }
      labels <- tmp@meta.data[[cluster_col]]
      tab <- table(
        species_aln_short = tmp$species_aln_short,
        cluster = labels
      )
      metrics <- cluster_metrics(tab)
      key <- paste(k.param, res, sep = "_")
      scan_rows[[key]] <- data.frame(
        k.param = k.param,
        resolution = res,
        metrics
      )
      cluster_labels[[key]] <- labels
      if (!is.na(stop_row_at_n_clusters) && metrics$n_clusters > stop_row_at_n_clusters) {
        message(
          "    ",
          method,
          ": stopping k=",
          k.param,
          " row after res=",
          res,
          " produced ",
          metrics$n_clusters,
          " clusters"
        )
        break
      }
    }
  }

  scan <- do.call(rbind, scan_rows)
  rownames(scan) <- NULL
  scan$distance <- abs(scan$n_clusters - expected_n_clusters)
  scan$acceptable <- scan$overall_purity >= min_acceptable_purity &
    scan$largest_cluster_fraction <= max_acceptable_largest_cluster_fraction

  exact <- scan[scan$n_clusters == expected_n_clusters, , drop = FALSE]
  exact_acceptable <- exact[exact$acceptable, , drop = FALSE]
  if (nrow(exact_acceptable) > 0) {
    selected_row <- exact_acceptable[
      order(
        -exact_acceptable$overall_purity,
        exact_acceptable$largest_cluster_fraction,
        -exact_acceptable$min_cluster_size,
        exact_acceptable$k.param,
        exact_acceptable$resolution
      ),
      ,
      drop = FALSE
    ][1, ]
  } else {
    candidates <- scan[scan$acceptable, , drop = FALSE]
    if (nrow(candidates) == 0) {
      candidates <- scan
    } else if (nrow(exact) > 0) {
      warning(
        method,
        " produced exactly ",
        expected_n_clusters,
        " clusters only in collapsed/low-purity solutions; selecting the best acceptable nearby clustering"
      )
    }
    selected_row <- candidates[
      order(
        candidates$distance,
        -candidates$overall_purity,
        candidates$largest_cluster_fraction,
        -candidates$min_cluster_size,
        candidates$k.param,
        candidates$resolution
      ),
      ,
      drop = FALSE
    ][1, ]
    warning(
      method,
      " did not produce exactly ",
      expected_n_clusters,
      " clusters; selected closest k.param ",
      selected_row$k.param,
      " and resolution ",
      selected_row$resolution,
      " with ",
      selected_row$n_clusters,
      " clusters"
    )
  }

  selected_k <- selected_row$k.param
  selected_resolution <- selected_row$resolution
  selected_key <- paste(selected_k, selected_resolution, sep = "_")
  adata$cluster_selected <- cluster_labels[[selected_key]]

  if (selected_row$overall_purity < 0.8 || selected_row$largest_cluster_fraction > 0.9) {
    warning(
      method,
      " selected a weak clustering: purity=",
      round(selected_row$overall_purity, 3),
      ", largest_cluster_fraction=",
      round(selected_row$largest_cluster_fraction, 3)
    )
  }

  list(
    object = adata,
    k.param = selected_k,
    resolution = selected_resolution,
    n_clusters = selected_row$n_clusters,
    overall_purity = selected_row$overall_purity,
    largest_cluster_fraction = selected_row$largest_cluster_fraction,
    min_cluster_size = selected_row$min_cluster_size,
    scan = scan
  )
}

select_reduction <- function(adata, method, preferred_reductions) {
  reduction <- preferred_reductions[preferred_reductions %in% names(adata@reductions)][1]
  if (is.na(reduction)) {
    stop(
      "None of the preferred reductions for ",
      method,
      " were found. Available reductions: ",
      paste(names(adata@reductions), collapse = ", ")
    )
  }
  reduction
}

run_method <- function(
    adata,
    method,
    preferred_reductions,
    outDir,
    plotDir,
    expected_n_species,
    max_resolution = NULL,
    max_k_param = NULL,
    max_dims = 30) {
  reduction <- select_reduction(adata, method, preferred_reductions)
  n_dims <- ncol(Embeddings(adata, reduction = reduction))
  dims <- seq_len(min(max_dims, n_dims))

  result <- select_resolution(
    adata = adata,
    method = method,
    reduction = reduction,
    dims = dims,
    expected_n_clusters = expected_n_species,
    max_resolution = max_resolution,
    max_k_param = max_k_param
  )

  tab <- table(
    species_aln_short = result$object$species_aln_short,
    cluster = result$object$cluster_selected
  )

  write.csv(
    result$scan,
    file.path(outDir, paste0("compare_clustering_", method, "_resolution_scan.csv")),
    row.names = FALSE
  )
  write.csv(
    confusion_to_df(tab, method, result$resolution),
    file.path(outDir, paste0("compare_clustering_", method, "_confusion.csv")),
    row.names = FALSE
  )
  ggsave(
    file.path(plotDir, paste0("compare_clustering_", method, "_confusion.pdf")),
    plot_confusion(tab, method, result$resolution),
    width = 7,
    height = 5,
    limitsize = FALSE
  )

  list(
    method = method,
    k.param = result$k.param,
    resolution = result$resolution,
    n_clusters = result$n_clusters,
    overall_purity = result$overall_purity,
    largest_cluster_fraction = result$largest_cluster_fraction,
    min_cluster_size = result$min_cluster_size,
    n_cells = ncol(result$object),
    confusion = tab
  )
}

run_dataset <- function(config) {
  dataset_name <- config$dataset_name
  expected_n_species <- as.integer(config$expected_n_species)
  bascetRoot <- file.path("/husky/henriksson/atrandi", dataset_name)
  plotDir <- file.path(getwd(), "plots", dataset_name)
  outDir <- file.path("compare_clustering_sim", dataset_name)

  if (!dir.exists(plotDir)) {
    dir.create(plotDir, recursive = TRUE)
  }
  if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = TRUE)
  }

  message("Running clustering comparison for ", dataset_name)

  adata_align <- read_cached_object(bascetRoot, "cache_adata_align.RDS")

  if (!all(c("species_aln", "species_aln_short") %in% colnames(adata_align@meta.data))) {
    stop("Alignment object must contain species_aln and species_aln_short metadata")
  }

  alignment_meta <- adata_align@meta.data[, c("species_aln", "species_aln_short"), drop = FALSE]
  rownames(alignment_meta) <- canonical_cell_names(rownames(alignment_meta))

  expected_species <- load_expected_species(config$expected_species_map)
  observed_species <- sort(unique(alignment_meta$species_aln_short))
  missing_alignment_species <- setdiff(expected_species, observed_species)
  n_truth_species <- length(unique(alignment_meta$species_aln_short))
  if (n_truth_species != expected_n_species) {
    warning(
      dataset_name,
      " alignment ground truth contains ",
      n_truth_species,
      " species, not the expected ",
      expected_n_species
    )
  }
  if (length(missing_alignment_species) > 0) {
    warning(
      dataset_name,
      " species missing from alignment ground truth: ",
      paste(missing_alignment_species, collapse = ", ")
    )
  }

  write.csv(
    data.frame(species = missing_alignment_species),
    file.path(outDir, "compare_clustering_missing_alignment_species.csv"),
    row.names = FALSE
  )

  requested_methods <- Sys.getenv("METHODS", "")
  method_names <- c("alignment", "kraken", "infokmer", "countsketch")
  if (nzchar(requested_methods)) {
    selected_methods <- trimws(strsplit(requested_methods, ",", fixed = TRUE)[[1]])
    unknown_methods <- setdiff(selected_methods, method_names)
    if (length(unknown_methods) > 0) {
      stop("Unknown METHODS entries: ", paste(unknown_methods, collapse = ", "))
    }
  } else {
    selected_methods <- method_names
  }

  method_defs <- list(
    alignment = function() run_method(
      adata_align,
      "alignment",
      c("lsi", "pca"),
      outDir = outDir,
      plotDir = plotDir,
      expected_n_species = expected_n_species,
      max_resolution = 0.2,
      max_k_param = 200,
      max_dims = 30
    ),
    kraken = function() {
      adata_kraken <- read_cached_object(bascetRoot, "cache_adata_kraken.RDS")
      adata_kraken <- add_alignment_truth(adata_kraken, alignment_meta, "kraken")
      run_method(
        adata_kraken,
        "kraken",
        c("lsi", "pca"),
        outDir = outDir,
        plotDir = plotDir,
        expected_n_species = expected_n_species,
        max_resolution = 0.5,
        max_k_param = 200,
        max_dims = 30
      )
    },
    infokmer = function() {
      adata_infokmer <- read_cached_object(bascetRoot, "cache_adata_infokmer.RDS")
      adata_infokmer <- add_alignment_truth(adata_infokmer, alignment_meta, "infokmer")
      run_method(
        adata_infokmer,
        "infokmer",
        c("lsi", "pca"),
        outDir = outDir,
        plotDir = plotDir,
        expected_n_species = expected_n_species,
        max_resolution = 0.5,
        max_k_param = 200,
        max_dims = 30
      )
    },
    countsketch = function() {
      adata_cs <- read_cached_object(bascetRoot, "cache_adata_cs.RDS")
      adata_cs <- add_alignment_truth(adata_cs, alignment_meta, "countsketch")
      run_method(
        adata_cs,
        "countsketch",
        c("kmersketch", "pca"),
        outDir = outDir,
        plotDir = plotDir,
        expected_n_species = expected_n_species,
        max_resolution = 2,
        max_k_param = 500,
        max_dims = 50
      )
    }
  )

  results <- lapply(selected_methods, function(method) method_defs[[method]]())
  names(results) <- selected_methods

  summary_df <- do.call(
    rbind,
    lapply(results, function(x) {
      data.frame(
        dataset = dataset_name,
        method = x$method,
        n_cells = x$n_cells,
        selected_k_param = x$k.param,
        selected_resolution = x$resolution,
        n_clusters = x$n_clusters,
        overall_purity = x$overall_purity,
        largest_cluster_fraction = x$largest_cluster_fraction,
        min_cluster_size = x$min_cluster_size,
        n_species_ground_truth = n_truth_species,
        missing_alignment_species = paste(missing_alignment_species, collapse = "; ")
      )
    })
  )

  summary_path <- file.path(outDir, "compare_clustering_sim_summary.csv")
  if (length(selected_methods) < length(method_names) && file.exists(summary_path)) {
    old_summary <- read.csv(summary_path, check.names = FALSE, na.strings = character())
    if (!"dataset" %in% colnames(old_summary)) {
      old_summary$dataset <- dataset_name
    }
    old_summary$dataset[is.na(old_summary$dataset) | old_summary$dataset == ""] <- dataset_name
    old_summary$missing_alignment_species[is.na(old_summary$missing_alignment_species)] <- ""
    old_summary <- old_summary[!old_summary$method %in% selected_methods, , drop = FALSE]
    old_missing <- setdiff(colnames(summary_df), colnames(old_summary))
    for (col in old_missing) {
      old_summary[[col]] <- NA
    }
    new_missing <- setdiff(colnames(old_summary), colnames(summary_df))
    for (col in new_missing) {
      summary_df[[col]] <- NA
    }
    summary_df <- rbind(
      old_summary[, colnames(summary_df), drop = FALSE],
      summary_df
    )
    summary_df <- summary_df[match(method_names, summary_df$method), , drop = FALSE]
  }

  write.csv(
    summary_df,
    summary_path,
    row.names = FALSE
  )

  print(summary_df)
  summary_df
}

requested_dataset <- Sys.getenv("DATASET_NAME", "")
if (nzchar(requested_dataset)) {
  if (!requested_dataset %in% names(dataset_configs)) {
    stop("Unknown DATASET_NAME: ", requested_dataset)
  }
  dataset_configs[[requested_dataset]]$enabled <- TRUE
  dataset_configs <- dataset_configs[requested_dataset]
}

enabled_configs <- Filter(function(x) isTRUE(x$enabled), dataset_configs)
if (length(enabled_configs) == 0) {
  stop("No datasets are enabled")
}

all_summaries <- do.call(rbind, lapply(enabled_configs, run_dataset))
write.csv(
  all_summaries,
  file.path("compare_clustering_sim", "compare_clustering_sim_summary_all.csv"),
  row.names = FALSE
)
