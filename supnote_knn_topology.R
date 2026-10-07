library(ggplot2)
library(patchwork)


################################################################################
###### Illustrative figure: all-to-all vs kNN graph on a 2D manifold ###########
################################################################################

set.seed(42)

#### Two arcs that share a 2D embedding — topology = "two curves".
n_per <- 30
t1 <- seq(0, pi, length.out = n_per)
t2 <- seq(0, pi, length.out = n_per)

arc1 <- data.frame(x = cos(t1),         y = sin(t1),         group = "A")
arc2 <- data.frame(x = cos(t2) + 1.6,   y = -sin(t2) + 0.3,  group = "B")

pts <- rbind(arc1, arc2)
pts$x <- pts$x + rnorm(nrow(pts), sd = 0.03)
pts$y <- pts$y + rnorm(nrow(pts), sd = 0.03)
pts$id <- seq_len(nrow(pts))
N <- nrow(pts)

D <- as.matrix(dist(pts[, c("x", "y")]))


#### All-to-all edges
all_edges <- expand.grid(i = 1:N, j = 1:N)
all_edges <- all_edges[all_edges$i < all_edges$j, ]
all_edges$dist <- D[cbind(all_edges$i, all_edges$j)]
all_edges$x    <- pts$x[all_edges$i]
all_edges$y    <- pts$y[all_edges$i]
all_edges$xend <- pts$x[all_edges$j]
all_edges$yend <- pts$y[all_edges$j]


#### kNN edges (k nearest neighbours per node, mutual-or-not)
k <- 4
knn_pairs <- list()
for (i in seq_len(N)) {
  nn <- order(D[i, ])[2:(k + 1)]   # skip self
  for (j in nn) {
    a <- min(i, j); b <- max(i, j)
    knn_pairs[[length(knn_pairs) + 1]] <- c(a, b)
  }
}
knn_df <- unique(do.call(rbind, knn_pairs))
knn_df <- data.frame(i = knn_df[, 1], j = knn_df[, 2])
knn_df$dist <- D[cbind(knn_df$i, knn_df$j)]
knn_df$x    <- pts$x[knn_df$i]
knn_df$y    <- pts$y[knn_df$i]
knn_df$xend <- pts$x[knn_df$j]
knn_df$yend <- pts$y[knn_df$j]


cat("Points:", N, "\n")
cat("All-to-all edges:", nrow(all_edges), "\n")
cat("kNN edges (k=", k, "): ", nrow(knn_df), "\n", sep = "")


#### Common aesthetic
base_theme <- theme_void(base_size = 11) +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
        legend.position = "none",
        plot.margin = margin(6, 6, 6, 6))

dist_range <- range(c(all_edges$dist, knn_df$dist))

p1 <- ggplot() +
  geom_segment(data = all_edges,
               aes(x = x, y = y, xend = xend, yend = yend, color = dist, alpha = 1 / (dist + 0.1)),
               linewidth = 0.2) +
  geom_point(data = pts, aes(x = x, y = y, fill = group),
             shape = 21, size = 3, stroke = 0.4, color = "white") +
  scale_color_viridis_c(option = "magma", direction = -1, limits = dist_range) +
  scale_alpha(range = c(0.05, 0.6)) +
  scale_fill_manual(values = c(A = "#2c7fb8", B = "#d95f0e")) +
  coord_equal() +
  labs(title = "All-to-all distances",
       subtitle = sprintf("%d edges  (N·(N-1)/2)", nrow(all_edges))) +
  base_theme

p2 <- ggplot() +
  geom_segment(data = knn_df,
               aes(x = x, y = y, xend = xend, yend = yend, color = dist),
               linewidth = 0.6) +
  geom_point(data = pts, aes(x = x, y = y, fill = group),
             shape = 21, size = 3, stroke = 0.4, color = "white") +
  scale_color_viridis_c(option = "magma", direction = -1, limits = dist_range) +
  scale_fill_manual(values = c(A = "#2c7fb8", B = "#d95f0e")) +
  coord_equal() +
  labs(title = sprintf("k-NN graph (k = %d)", k),
       subtitle = sprintf("%d edges", nrow(knn_df))) +  #  —  topology preserved
  base_theme


fig <- p1 + p2 + plot_layout(ncol = 2)

dir.create("plots", showWarnings = FALSE)
#out_svg <- "plots/knn_vs_full.svg"
out_png <- "plots/knn_vs_full.png"
#ggsave(out_svg, fig, width = 9, height = 4.5)
ggsave(out_png, fig, width = 9, height = 4.5, dpi = 150)


############


fig <- p1 + p2 + plot_layout(ncol = 1)

dir.create("plots", showWarnings = FALSE)
out_png <- "plots/knn_vs_full.png"
ggsave(out_png, fig, width = 4.5, height = 9, dpi = 150)