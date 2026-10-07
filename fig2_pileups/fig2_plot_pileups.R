library(duckdb)
library(DBI)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridisLite)
library(patchwork)
library(ggh4x)
library(farver)
library(ggtext)
library(ggrastr)
library(tibble)
library(ggpattern)

source("style.R")

if(FALSE) {
  source("make_pileups.R")
}



# ── Constants ─────────────────────────────────────────────────────────────────

SPECIES_MAP <- c(
  "ATCC_10987"   = "Bacillus pacificus",
  "ATCC_12228"   = "Staphylococcus epidermidis",
  "ATCC_15703"   = "Bifidobacterium adolescentis",
  "ATCC_17029"   = "Cereibacter sphaeroides",
  "ATCC_33323"   = "Lactobacillus gasseri",
  "ATCC_35702"   = "Clostridium beijerinckii",
  "ATCC_47077"   = "Enterococcus faecalis",
  "ATCC_700610"  = "Streptococcus mutans",
  "ATCC_700926"  = "Escherichia coli",
  "ATCC_BAA-816" = "Deinococcus radiodurans",
  "NZ_CDNB"      = "Xanthomonas campestris",
  "phix"         = "PhiX"
)

if(FALSE) {
  SPECIES_TOPLOT <- setdiff(unname(SPECIES_MAP), c("PhiX", "Xanthomonas campestris"))
} else {
  SPECIES_TOPLOT <- c("Enterococcus faecalis")
}
PLOT_DIR <- "plots"

dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)


# ── Annotate ──────────────────────────────────────────────────────────────────

pileup <- pileup.raw %>%
  dplyr::rename(weighted.depth = weighted_depth) %>%
  mutate(species = SPECIES_MAP[organism_id]) %>%
  filter(species != "PhiX") %>%
  filter(species %in% SPECIES_TOPLOT)

contigs <- contigs %>%
  left_join(genome.sizes, by = "organism_id") %>%
  mutate(
    species  = SPECIES_MAP[organism_id],
    bp.start = contig_offset,
    bp.end   = contig_offset + sizeof_reference
  )

# ── Cell ordering ─────────────────────────────────────────────────────────────

pileup <- pileup %>%
  group_by(group, id, species) %>%
  mutate(
    mean.depth  = mean(weighted.depth),
    total.depth = sum(weighted.depth)
  ) %>%
  ungroup()

pileup.ranks <- pileup %>%
  distinct(group, id, species, mean.depth, total.depth) %>%
  arrange(group, species, total.depth) %>%
  group_by(group, species) %>%
  mutate(rank = row_number()) %>%
  select(group, id, species, rank) %>%
  ungroup()

pileup <- pileup %>%
  left_join(pileup.ranks, by = join_by(group, id, species))

# ── Rescale to per-group per-species ──────────────────────────────────────────

pileup <- pileup %>%
  group_by(group, species, position) %>%
  mutate(
    cell.log.depth = log10(weighted.depth + 1),
    col_sum        = sum(cell.log.depth),
    depth_colnorm  = cell.log.depth / col_sum
  ) %>%
  ungroup() %>%
  group_by(group, species) %>%
  mutate(
    log.depth       = log10(col_sum + 1),
    norm.depth      = log.depth / max(log.depth),
    cell.norm.depth = cell.log.depth / max(cell.log.depth)
  ) %>%
  ungroup()

# ── Species order + palette ───────────────────────────────────────────────────

species.order  <- sort(setdiff(unname(SPECIES_MAP), c("PhiX", "Xanthomonas campestris")))
muted.colors   <- c(as.vector(khroma::colour("muted")(length(species.order) - 1)), "#DDDDDD")
species.colors <- setNames(muted.colors, species.order)
species.oklab  <- convert_colour(t(col2rgb(muted.colors)), from = "rgb", to = "oklab")
rownames(species.oklab) <- species.order

white.oklab <- matrix(
  convert_colour(matrix(c(255, 255, 255), nrow = 1), from = "rgb", to = "oklab"),
  nrow = length(species.order), ncol = 3, byrow = TRUE,
  dimnames = list(species.order, c("l", "a", "b"))
)
oklab.stops <- list(white.oklab, species.oklab, species.oklab)
#oklab.stops <- list(species.oklab, species.oklab, species.oklab) ####################### disabling oklab (does not work)

scale_fill_species   <- function(...) scale_fill_manual(values = species.colors, ...)
scale_colour_species <- function(...) scale_colour_manual(values = species.colors, ...)

# ── Heatmap fill via OKLab interpolation ─────────────────────────────────────

pileup <- pileup %>%
  mutate(
    fill = {
      idx       <- match(species, species.order)
      n.segs    <- length(oklab.stops) - 1L
      seg.width <- 1 / n.segs
      
      seg   <- pmin(floor(cell.norm.depth / seg.width), n.segs - 1L)
      t.seg <- (cell.norm.depth - seg * seg.width) / seg.width
      
      new.l <- new.a <- new.b <- numeric(length(idx))
      
      for (s in seq_len(n.segs)) {
        mask <- seg == (s - 1L)
        if (!any(mask)) next
        lo <- oklab.stops[[s]]
        hi <- oklab.stops[[s + 1L]]
        t  <- t.seg[mask]
        i  <- idx[mask]
        new.l[mask] <- lo[i, "l"] + t * (hi[i, "l"] - lo[i, "l"])
        new.a[mask] <- lo[i, "a"] + t * (hi[i, "a"] - lo[i, "a"])
        new.b[mask] <- lo[i, "b"] + t * (hi[i, "b"] - lo[i, "b"])
      }
      
      new.rgb <- convert_colour(cbind(l = new.l, a = new.a, b = new.b), from = "oklab", to = "rgb")
      new.rgb <- pmin(pmax(new.rgb, 0), 255)
      rgb(new.rgb[, 1], new.rgb[, 2], new.rgb[, 3], maxColorValue = 255)
    }
  )


# ── Axis limits ───────────────────────────────────────────────────────────────

bottom.limits <- pileup %>%
  group_by(group, species, position) %>%
  summarise(total.depth = sum(weighted.depth), .groups = "drop") %>%
  mutate(total.depth.t = log10(total.depth + 1e-9)) %>%
  group_by(species) %>%
  summarise(ymax = max(total.depth.t)) %>%
  deframe()

side.limits <- pileup %>%
  distinct(group, id, species, total.depth) %>%
  mutate(total.depth.t = log10(total.depth + 1)) %>%
  group_by(species) %>%
  summarise(xmax = max(total.depth.t)) %>%
  deframe()

# ── Plot functions ────────────────────────────────────────────────────────────

pileup.species.heatmap <- function(g, sp) {
  df.heat <- pileup %>%
    filter(species == sp & group == g)
  
  x.breaks.max <- max(df.heat$position)
  x.breaks <- c(
    scales::breaks_width(500000)(range(df.heat$position)),
    x.breaks.max
  ) %>%
    sort() %>%
    unique() %>%
    .[. <= x.breaks.max]
  
  ggplot(df.heat, aes(x = position, y = rank, fill = fill)) +
    rasterise(geom_raster(), dpi = 300, dev = "ragg_png") +
    scale_fill_identity() +
    scale_x_continuous(
      breaks = x.breaks,
      labels = function(x) ifelse(x == x.breaks.max, "", paste0(x / 1e6, " Mbp")),
      expand = c(0, 0)
    ) +
    scale_y_continuous(breaks = c(), expand = c(0, 0)) +
    labs(x = "", y = "") +
    theme_forpub(font_scale = 2) +
    theme(
      axis.line          = element_blank(),
      axis.line.x.bottom = element_line(colour = "#000000", linewidth = 0.25)
    )
}

pileup.species.bottom.hist <- function(g, sp) {
  n.cells <- pileup %>%
    filter(species == sp & group == g) %>%
    distinct(id) %>%
    nrow()
  
  df.bhist <- pileup %>%
    filter(species == sp & group == g) %>%
    group_by(position) %>%
    summarise(total.depth = sum(norm.depth), .groups = "drop") %>%
    mutate(total.depth.t = log10(total.depth + 1))
  
  ggplot(df.bhist, aes(x = position, y = total.depth.t)) +
    rasterise(geom_area(fill = species.colors[sp], position = "identity"), dpi = 300, dev = "ragg_png") +
    rasterise(geom_line(color = species.colors[sp], linewidth = 0.3), dpi = 300, dev = "ragg_png") +
    scale_x_continuous(
      breaks = scales::breaks_width(500000),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      position = "right",
      expand   = c(0, 0),
      limits   = c(0, bottom.limits[sp]),
      breaks   = c(0, 2, 4)#,
#      labels   = function(x) parse(text = paste0("10^", x))
    ) +
#    annotation_logticks(
#      sides = "r", outside = FALSE,
#      short = unit(4, "pt"), mid = unit(8, "pt"), long = unit(12, "pt"),
#      linewidth = 0.4
#    ) +
    coord_cartesian(clip = "off") +
    labs(title = paste0(g, " (n = ", n.cells, ")"), x = "", y = "") +
    theme_forpub(font_scale = 2) +
    theme(
      plot.title   = element_text(size = 24),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line    = element_blank()
    )
}

pileup.species.side.hist <- function(g, sp) {
  df.shist <- pileup %>%
    filter(species == sp & group == g) %>%
    distinct(id, rank, total.depth) %>%
    mutate(total.depth.t = log10(total.depth + 1))
  
  ggplot(df.shist, aes(x = total.depth.t, y = rank)) +
    rasterise(geom_area(fill = species.colors[sp], orientation = "y"), dpi = 300, dev = "ragg_png") +
    rasterise(geom_line(color = species.colors[sp], linewidth = 0.3), dpi = 300, dev = "ragg_png") +
    scale_x_continuous(
      expand   = c(0, 0),
      position = "top",
      limits   = c(0, side.limits[sp]),
      breaks   = c(0, 2, 4)#,
#      labels   = function(x) parse(text = paste0("10^", x))
    ) +
    scale_y_continuous(breaks = c(), expand = c(0, 0)) +
#    annotation_logticks(
#      sides = "t", outside = FALSE,
#      short = unit(4, "pt"), mid = unit(8, "pt"), long = unit(12, "pt"),
#      linewidth = 0.4
#    ) +
    coord_cartesian(clip = "off") +
#    labs(x = expression(log[10](depth)), y = "") +
    labs(x = expression(depth), y = "") +
    theme_forpub(font_scale = 2) +
    theme(
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line    = element_blank()
    )
}

# ── Render ────────────────────────────────────────────────────────────────────

pileup.species.block <- function(g, sp) {
  p.bhist <- pileup.species.bottom.hist(g, sp) + theme(plot.margin = margin(0, 0, 0, 0))
  p.heat  <- pileup.species.heatmap(g, sp)     + theme(plot.margin = margin(0, 0, 0, 0))
  p.shist <- pileup.species.side.hist(g, sp)   + theme(plot.margin = margin(0, 0, 0, 0))
  
  design <- "
    A##
    ###
    C#D
  "
  
  wrap_plots(
    A = p.bhist,
    C = p.heat,
    D = p.shist,
    design  = design,
    widths  = c(3, -2/3 + 0.05, 1),
    heights = c(1, -2/3 + 0.05, 3)
  )
}


############# Make individual pileups
for (sp in SPECIES_TOPLOT) {
  if (pileup %>% filter(species == sp, group %in% c("MDA", "PTA")) %>% nrow() == 0) {
    message("Skipping ", sp, " no data"); next
  }
  
  p_heatmap <- (pileup.species.block("MDA", sp) | pileup.species.block("PTA", sp)) +
    theme(panel.border = element_blank())
  
  sp.filename <- gsub(" ", "_", tolower(sp))
  ggsave(file.path(PLOT_DIR, paste0("pileup_", sp.filename, ".pdf")), p_heatmap, width = 24, height = 12)
  ggsave(file.path(PLOT_DIR, paste0("pileup_", sp.filename, ".svg")), p_heatmap, width = 24, height = 12)
}

p_heatmap #666

###################### 
###################### 
###################### 
###################### 


pileup.species.heatmap.all <- function(g, sp, minrank=5000) {
  df.heat <- pileup %>%
    filter(species == sp & group == g & rank > minrank) # & rank<maxrank
  
  x.breaks.max <- max(df.heat$position)
  x.breaks <- c(
    scales::breaks_width(500000)(range(df.heat$position)),
    x.breaks.max
  ) %>%
    sort() %>%
    unique() %>%
    .[. <= x.breaks.max]
  
  ggplot(df.heat, aes(x = position, y = rank, fill = fill)) +
    rasterise(geom_raster(), dpi = 300, dev = "ragg_png") +
    scale_fill_identity() +
    scale_x_continuous(
      breaks = x.breaks,
      labels = function(x) ifelse(x == x.breaks.max, "", paste0(x / 1e6, " Mbp")),
      expand = c(0, 0)
    ) +
    scale_y_continuous(breaks = c(), expand = c(0, 0)) +
    labs(x = "", y = "") +
    theme_forpub(font_scale = 2) +
    theme(
      axis.line          = element_blank(),
      axis.line.x.bottom = element_line(colour = "#000000", linewidth = 0.25)
    )
}

pileup.species.block.all <- function(g, sp) {
  #pileup.species.heatmap.all(g, sp)     + theme(plot.margin = margin(0, 0, 0, 0))
  
  
  
  p.heat  <- pileup.species.heatmap(g, sp)     + theme(plot.margin = margin(0, 0, 0, 0))
  p.shist <- pileup.species.side.hist(g, sp)   + theme(plot.margin = margin(0, 0, 0, 0))
  
  design <- "
    C#D
  "
  wrap_plots(
    C = p.heat,
    D = p.shist,
    design  = design,
    widths  = c(3, -2/3 + 0.05, 1),
    heights = c(3)
  )  
  
  
  
  
  
  
  
}

###################### Make all pileups in one figure
if(FALSE) {

  all_heatmap <- list()
  for (sp in SPECIES_TOPLOT) {
    cat(cp)
    if (pileup %>% filter(species == sp, group %in% c("MDA", "PTA")) %>% nrow() == 0) {
      next
    }

    p_heatmap <- (pileup.species.block.all("MDA", sp) | pileup.species.block.all("PTA", sp)) +
      theme(panel.border = element_blank())
    p_heatmap

    sp.filename <- gsub(" ", "_", tolower(sp))
    ggsave(file.path(PLOT_DIR, paste0("new_pileup_", sp.filename, ".pdf")), p_heatmap, width = 24, height = 4)
    #ggsave(file.path(PLOT_DIR, paste0("pileup_", sp.filename, ".svg")), p_heatmap, width = 24, height = 12)
  }

}


# ── Single cell pileups ───────────────────────────────────────────────────────

pileup.single.cells <- pileup %>%
  distinct(group, species, id, total.depth) %>%
  group_by(group, species) %>%
  mutate(
    ntile = case_when(
      total.depth > quantile(total.depth, 0.99) ~ "top_01",
      #total.depth > quantile(total.depth, 0.75) ~ "top_75",
      TRUE ~ NA
    )) %>%
  filter(!is.na(ntile)) %>%
  group_by(group, species, ntile) %>%
  arrange(group, species, total.depth) %>%
  dplyr::slice(ceiling(n() / 2)) %>%
  select(group, species, id, total.depth, ntile) %>%
  ungroup()

pileup.single <- function(cells) {
  sp <- unique(cells$species)
  
  df <- pileup %>%
    inner_join(cells %>% select(group, species, id), by = c("group", "species", "id")) %>%
    complete(
      nesting(group, species, id),
      position = seq(min(position), max(position), by = unique(bin_width)[1])
    ) %>%
    mutate(weighted.depth = replace_na(weighted.depth, 0)) %>%
    group_by(group, species, id) %>%
    fill(total.depth, .direction = "downup") %>%
    mutate(cell_label = paste0(group, " (n = ", round(total.depth, digits = 0), ")")) %>%
    arrange(group, position)
  print(df)
  
  ggplot(df, aes(x = position, y = log10(weighted.depth + 1))) +
    geom_area(fill = species.colors[sp]) +
    geom_line(color = species.colors[sp], linewidth = 0.3) +
    facet_grid(. ~ cell_label, axis.labels = "all") +
    scale_x_continuous(
      breaks = scales::breaks_width(500000),
      labels = function(x) paste0(x / 1e6, "Mbp"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand   = c(0, NA),
      position = "right"#,
      #labels   = function(x) parse(text = paste0("10^", x))
    ) +
#    annotation_logticks(
#      sides = "r", outside = FALSE,
#      short = unit(3, "pt"), mid = unit(8, "pt"), long = unit(10, "pt"),
#      linewidth = 0.4
#    ) +
    labs(title = cells$ntile[1], x = "", y = "Read count") +
    theme_forpub(base_size = 24, font_scale = 1) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.ticks.y = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )
}

for (sp in SPECIES_TOPLOT) {
  for (nt in unique(pileup.single.cells$ntile)) {
    cells <- pileup.single.cells %>% filter(species == sp, ntile == nt)
    if (nrow(cells) == 0) { message("Skipping ", sp, " ntile ", nt); next }
    
    plot_sc <- pileup.single(cells) + ggtitle("")
    plot_sc
    sp.filename <- gsub(" ", "_", tolower(sp))
    ggsave(file.path(PLOT_DIR, paste0("pileup_", sp.filename, "_single_", nt, ".pdf")), plot_sc, width = 12, height = 3)
    ggsave(file.path(PLOT_DIR, paste0("pileup_", sp.filename, "_single_", nt, ".svg")), plot_sc, width = 12, height = 3)
  }
}




############## Special montage for Fig 2 - single cell pileup
if(TRUE) {
  nt <- "top_01"
  sp <- "Enterococcus faecalis"
  
  plot_sc <- pileup.single(cells) + ggtitle("")
  plot_sc
  sp.filename <- gsub(" ", "_", tolower(sp))
  ggsave(file.path(PLOT_DIR, paste0("pileup_", sp.filename, "_single_", nt, ".pdf")), plot_sc, width = 12, height = 3)
  ggsave(file.path(PLOT_DIR, paste0("pileup_", sp.filename, "_single_", nt, ".svg")), plot_sc, width = 12, height = 3)
}



##############
############## Global read distribution
############## 
# 
# if(FALSE) {
#   df.dist <- pileup %>%
#     distinct(group, id, species, total.depth) %>%
#     group_by(group, species) %>%
#     summarise(total.depth = sum(total.depth), .groups = "drop") %>%
#     mutate(
#       total.depth.t = log10(total.depth),
#       species       = factor(species, levels = species.order),
#       group         = factor(group, levels = c("PTA", "MDA"))
#     )
#   
#   p.dist <- ggplot(df.dist, aes(x = total.depth.t, y = species, fill = species)) +
#     rasterise(geom_col(
#       aes(group = interaction(species, group)),
#       position = position_dodge(width = 0.8),
#       width    = 0.7,
#     ), dev = "ragg_png", dpi = 300) +
#     scale_fill_species(name = NULL) +
#     scale_y_discrete(name = NULL, limits = rev) +
#     coord_cartesian(clip = "off") +
#     theme_forpub(base_size = 24, font_scale = 1) +
#     theme(
#       axis.text.y  = element_blank(),
#       panel.border = element_blank()
#     )
#   
#   ggsave(file.path(PLOT_DIR, "pileup_global_read_distribution.pdf"), p.dist, width = 12, height = 14)  
# }
