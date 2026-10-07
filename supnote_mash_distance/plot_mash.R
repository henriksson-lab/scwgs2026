library(ggplot2)
library(stringr)
library(tidyr)
library(dplyr)


mat_path <- "mash_out/mash_distance_matrix.tsv"
cat("Reading:", mat_path, "\n")

dm <- read.table(mat_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
dm <- as.matrix(dm)

#### Strip path and extension from labels; keep cell + coverage tag.
clean_label <- function(x) {
  x <- basename(x)
  x <- str_remove(x, "_R[12]\\.fastq(\\.gz)?$")
  x <- str_remove(x, "\\.fastq(\\.gz)?$")
  x <- str_remove(x, "\\.msh$")
  x
}
rownames(dm) <- clean_label(rownames(dm))
colnames(dm) <- clean_label(colnames(dm))

#### Extract coverage from "cell00001_cov0.0100x" → 0.01 for ordering.
extract_cov <- function(x) {
  m <- str_match(x, "cov([0-9.]+)x")[, 2]
  as.numeric(m)
}
cov_row <- extract_cov(rownames(dm))
cov_col <- extract_cov(colnames(dm))

ord_r <- order(cov_row)
ord_c <- order(cov_col)
dm <- dm[ord_r, ord_c]

cat("Matrix:", nrow(dm), "x", ncol(dm), "\n")
cat("Distance range:", range(dm, na.rm = TRUE), "\n")

#### Long format for ggplot
df <- as.data.frame(dm)
df$row <- rownames(df)
df_long <- pivot_longer(df, -row, names_to = "col", values_to = "dist")

df_long$row <- factor(df_long$row, levels = rownames(dm))
df_long$col <- factor(df_long$col, levels = colnames(dm))

p <- ggplot(df_long, aes(x = col, y = row, fill = dist)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Mash dist") +
  scale_x_discrete(breaks = levels(df_long$col)[seq(1, ncol(dm), by = 10)]) +
  scale_y_discrete(breaks = levels(df_long$row)[seq(1, nrow(dm), by = 10)]) +
  coord_fixed() +
  labs(title = "Mash distance — single species, log-scale coverage",
       x = "cell (sorted by coverage)",
       y = "cell (sorted by coverage)") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("plots/mash_heatmap.svg", p, width = 8, height = 7)
ggsave("plots/mash_heatmap.png", p, width = 8, height = 7, dpi = 150)
cat("wrote plots/mash_heatmap.{svg,png}\n")


#### Companion plot: distance vs coverage of the lower-coverage cell in the pair,
#### compared against the highest-coverage cell (best self-reference).
ref_idx <- which.max(cov_col)
ref_name <- colnames(dm)[ref_idx]
dist_to_ref <- dm[, ref_idx]
df_ref <- data.frame(
  cell = rownames(dm),
  coverage = cov_row[ord_r],
  dist = dist_to_ref
)

p2 <- ggplot(df_ref, aes(x = coverage, y = dist)) +
  geom_point(size = 1.2) +
  geom_line(alpha = 0.4) +
  scale_x_log10() +
  labs(title = paste0("Mash distance vs coverage (reference = ", ref_name, ")"),
       x = "cell coverage (x, log10)",
       y = "Mash distance to highest-coverage cell") +
  theme_minimal(base_size = 11)

ggsave("plots/mash_dist_vs_cov.svg", p2, width = 7, height = 5)
ggsave("plots/mash_dist_vs_cov.png", p2, width = 7, height = 5, dpi = 150)
cat("wrote plots/mash_dist_vs_cov.{svg,png}\n")
