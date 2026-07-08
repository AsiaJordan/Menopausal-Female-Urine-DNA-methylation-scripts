#!/usr/bin/env Rscript
#
# _style.R — shared plot theme + palettes for the menopause urine project.
#
# Usage (from any script, same directory or one level below):
#   source(file.path(dirname(sys.frame(1)$ofile), "_style.R"))
#   ... %>% + project_theme() + scale_colour_manual(values = STATUS_COLORS)
#
# Conventions:
#   - Plots get NO title and NO subtitle. Every figure is captioned in the
#     manuscript; an in-plot title is noise. If you must, use a short one-line
#     title, never a multiline subtitle.
#   - Axis labels carry the full information ("PC1 (95.4% variance)" etc.).
#   - Colours: STATUS_COLORS (Pre/Post), CELL_COLORS (cell types).

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2", repos = "https://cloud.r-project.org")
  library(ggplot2)
  if (!requireNamespace("egg", quietly = TRUE))
    install.packages("egg", repos = "https://cloud.r-project.org")
})

#' Publication theme wrapping egg::theme_article()
#'
#' Falls back to a theme_bw-derived replica if the `egg` package isn't
#' available, so scripts still run in minimal environments.
project_theme <- function(base_size = 10, base_family = "") {
  thm <- tryCatch(
    egg::theme_article(base_size = base_size, base_family = base_family),
    error = function(e) {
      theme_bw(base_size = base_size, base_family = base_family) +
        theme(panel.grid = element_blank(),
              panel.border = element_rect(colour = "black", fill = NA,
                                            linewidth = 0.5),
              strip.background = element_rect(fill = NA, colour = NA))
    }
  )
  thm +
    theme(
      plot.title       = element_text(size = base_size,
                                        face  = "plain",
                                        hjust = 0,
                                        margin = margin(b = 4)),
      plot.subtitle    = element_blank(),
      plot.caption     = element_text(size = base_size - 2,
                                        colour = "grey40",
                                        hjust = 0),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1,
                                        colour = "grey15"),
      legend.title     = element_text(size = base_size - 1),
      legend.text      = element_text(size = base_size - 1),
      legend.key.height = grid::unit(0.4, "lines"),
      legend.key.width  = grid::unit(0.8, "lines"),
      legend.margin     = margin(0, 0, 0, 0),
      strip.text        = element_text(size = base_size - 1,
                                         face = "plain"),
      plot.margin       = margin(4, 6, 4, 4),
      # Explicit OPAQUE backgrounds. egg::theme_article() and many base
      # themes leave plot.background = element_rect(fill = NA), which
      # ggsave writes to PNG as fully-transparent pixels — those render
      # as BLACK in Preview / VS Code / dark-mode viewers. Force white.
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key        = element_rect(fill = "white", colour = NA),
      strip.background  = element_rect(fill = "white", colour = NA)
    )
}

#' Save a ggplot to PNG with an opaque white background by default.
#'
#' Drop-in wrapper around ggplot2::ggsave() that forces bg = "white".
#' Older ggplot2 versions (< 3.4.0) default bg to theme background, which
#' egg::theme_article() sets to NA → transparent PNG → black in viewers.
ggsave_white <- function(filename, plot = ggplot2::last_plot(), ..., bg = "white") {
  ggplot2::ggsave(filename = filename, plot = plot, ..., bg = bg)
}

# ── Palettes ──────────────────────────────────────────────────────────────
STATUS_COLORS <- c(Post = "#c1272d", Pre = "#1f8a70")

CELL_COLORS <- c(
  Epi    = "#2e7d32",
  Fib    = "#f2a900",
  IC     = "#c1272d",
  EC     = "#8e24aa",
  # blood subtypes
  Neutro = "#c1272d",
  Mono   = "#6a4b9b",
  CD4T   = "#3a6ea5",
  CD8T   = "#4a90d9",
  B      = "#66a61e",
  NK     = "#e6ab02",
  Eosino = "#a6761d"
)

CHIP_SHAPES <- c("207758450138" = 16, "208356980025" = 17)

# Convenient shortcut for labs() with just a (short) title.
short_labs <- function(...) {
  args <- list(...)
  args$subtitle <- NULL   # never emit a subtitle
  do.call(labs, args)
}
