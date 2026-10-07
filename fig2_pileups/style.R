library(showtext)
library(dplyr)

font_add(
  "Arial",
  regular    = "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
  bold       = "/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf",
  italic     = "/usr/share/fonts/truetype/msttcorefonts/Arial_Italic.ttf",
  bolditalic = "/usr/share/fonts/truetype/msttcorefonts/Arial_Bold_Italic.ttf"
)
showtext_auto()

colors_forpub <- tribble(
  ~group,    ~dark,     ~low,     ~saturated,~mid,       ~mid_high,  ~high,
  "olive",  "#36461A", "#66771E", "#9BA415", "#C9C700",  "#DFDC67",  "#F4F1B3",
  "green",  "#163B1C", "#227130", "#459434", "#62B347",  "#A0CA79",  "#D9E8C6",
  "teal",   "#00394E", "#016879", "#019AA3", "#4BBCBD",  "#97D2D4",  "#CCE7EE",
  "blue",   "#012C5C", "#014E91", "#0272B2", "#5799D1",  "#9DCBEC",  "#C8E7FB",
  "purple", "#481951", "#792C74", "#A84E94", "#BB7CB4",  "#D3AAD1",  "#EBD6E9",
  "red",    "#741915", "#9C2826", "#C93E3F", "#DE6866",  "#EDA4A7",  "#FAD0CE",
  "orange", "#843200", "#B74F06", "#EC6F00", "#F59A45",  "#FBC07F",  "#FDE0BD",
  "yellow", "#6B5513", "#9C7717", "#CCA02C", "#EBC850",  "#F6DC8A",  "#FFF1C3",
  "skin",   "#442D1F", "#755040", "#936A57", "#BF997D",  "#DDBDA3",  "#F8E6D7",
  "clay",   "#3C3B2E", "#625F4A", "#888366", "#A8A386",  "#C6C2A6",  "#E4E0CE",
  "grey",   "#010E21", "#253247", "#49566D", "#6F7B91",  "#99A3B4",  "#C8CEDA"
)

nature_pal <- function(palette = "saturated") {
  if (!palette %in% names(colors_forpub)) stop(
    "Unknown palette '", palette, "'. Choose from: ",
    paste(setdiff(names(colors_forpub), "group"), collapse = ", ")
  )
  colors <- setNames(colors_forpub[[palette]], colors_forpub$group)
  function(n) {
    if (n > length(colors)) stop(
      "Nature palette has ", length(colors), " colours, ", n, " requested."
    )
    unname(colors[seq_len(n)])
  }
}

scale_color_nature <- function(palette = "saturated", ...) {
  discrete_scale("colour", nature_pal(palette), ...)
}

scale_fill_nature <- function(palette = "saturated", ...) {
  discrete_scale("fill", nature_pal(palette), ...)
}

theme_forpub <- function(base_size = 24, font_scale = 1) {
  title_size <- 24 * font_scale
  text_size <- 20 * font_scale
  axis_text_size <- 14 * font_scale
  axis_title_size <- 18 * font_scale
  legend_title_size <- 14 * font_scale
  legend_text_size <- 14 * font_scale
  
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank(),
      panel.background    = element_rect(fill = "#FFFFFF", colour = NA),
      plot.background     = element_rect(fill = "#FFFFFF", colour = NA),
      strip.background    = element_blank(),
      legend.background   = element_rect(fill = "#FFFFFF", colour = NA),
      legend.key          = element_rect(fill = "#FFFFFF", colour = NA),
      legend.title        = element_text(
        size = legend_title_size, colour = "#000000",
        margin = margin(t = 2.5, b = 2.5), hjust = 0),
      legend.text         = element_text(size = legend_text_size, colour = "#000000"),
      legend.spacing.y    = unit(0, "cm"),
      legend.box.margin   = margin(),
      legend.margin       = margin(),
      text                = element_text(colour = "#000000", size = text_size),
      axis.text           = element_text(colour = "#000000", size = axis_text_size),
      axis.title          = element_text(
        colour = "#000000", size = axis_title_size,
        margin = margin(t = 1, r = 1, b = 1, l = 1)),
      plot.title          = element_text(size = title_size, margin = margin(b = 5)),
      plot.title.position = "plot",
      axis.line  = element_line(colour = "#000000", linewidth = base_size / 48),
      axis.ticks = element_line(colour = "#000000", linewidth = base_size / 48),
      panel.border = element_blank()
    )
}