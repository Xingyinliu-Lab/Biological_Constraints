#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2); library(ggforce)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: draw_fig3_panelA_rna05_rna06.R <fig_dir>")
fig_dir <- args[[1]]; dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

ink <- "#172433"; muted <- "#617180"; backbone <- "#7F8A94"; pair_col <- "#8EB7D2"
blue_dark <- "#1D5891"; green <- "#348864"; red <- "#C7514C"; purple <- "#7451A6"

rounded_box <- function(xmin, xmax, ymin, ymax, fill, colour, radius = 2.5, linewidth = 0.42) {
  d <- data.frame(x = c(xmin, xmax, xmax, xmin), y = c(ymin, ymin, ymax, ymax), id = 1)
  ggforce::geom_shape(data = d, aes(x, y, group = id), fill = fill, colour = colour,
                      linewidth = linewidth, radius = grid::unit(radius, "mm"))
}

hairpin <- function(cx, base_y, stem_h = 1.20, stem_w = 0.95, loop_h = 0.68,
                    n_stem = 6, n_loop = 9, id = "hp") {
  lx <- cx - stem_w / 2; rx <- cx + stem_w / 2
  sy <- seq(base_y, base_y + stem_h, length.out = n_stem)
  th <- seq(pi, 0, length.out = n_loop)
  lp <- data.frame(x = cx + (stem_w / 2) * cos(th), y = base_y + stem_h + loop_h * sin(th))
  path <- rbind(data.frame(x = lx - 0.28, y = base_y - 0.20),
                data.frame(x = rep(lx, n_stem), y = sy), lp[-c(1, nrow(lp)), ],
                data.frame(x = rep(rx, n_stem), y = rev(sy)),
                data.frame(x = rx + 0.28, y = base_y - 0.20))
  path$id <- id
  pairs <- data.frame(x = lx, xend = rx, y = sy, yend = sy, id = id)
  list(path = path, pairs = pairs, sy = sy, lx = lx, rx = rx)
}

draw_hairpin <- function(hp, node_size = 2.0, line_width = 0.50) list(
  geom_path(data = hp$path, aes(x, y, group = id), colour = backbone,
            linewidth = line_width, lineend = "round"),
  geom_segment(data = hp$pairs, aes(x, y, xend = xend, yend = yend),
               colour = pair_col, linewidth = 0.36, linetype = "22"),
  geom_point(data = hp$path, aes(x, y), shape = 21, size = node_size,
             stroke = 0.35, fill = "#F7F9FA", colour = backbone))

left_hp <- hairpin(2.90, 2.05, 1.18, 1.06, 0.67, id = "left")
left_p <- data.frame(x = left_hp$lx, y = left_hp$sy[4]); left_u <- data.frame(x = 2.90, y = 3.90)
partner_hp <- hairpin(11.35, 1.95, 1.32, 1.06, 0.72, n_stem = 7, n_loop = 11, id = "partner")
nonpartner_hp <- hairpin(15.45, 1.95, 1.32, 1.06, 0.72, n_stem = 7, n_loop = 11, id = "nonpartner")
rna06_paths <- rbind(partner_hp$path, nonpartner_hp$path)
rna06_pairs <- rbind(partner_hp$pairs, nonpartner_hp$pairs)
pair_y06 <- partner_hp$sy[4]
anchor <- data.frame(x = partner_hp$lx, y = pair_y06)
true_partner <- data.frame(x = partner_hp$rx, y = pair_y06)
nonpartner <- data.frame(x = nonpartner_hp$lx, y = nonpartner_hp$sy[4])

p <- ggplot() +
  rounded_box(0.20, 8.88, 0.28, 5.30, "#F3F8FC", "#86B8DF", 3.2, 0.48) +
  rounded_box(9.12, 17.80, 0.28, 5.30, "#F8F5FC", "#B79ED8", 3.2, 0.48) +
  annotate("text", x = 0.55, y = 4.96, label = "RNA05", hjust = 0, family = "Arial", fontface = "bold", size = 3.15, colour = blue_dark) +
  annotate("text", x = 1.88, y = 4.96, label = "|  Pairing State", hjust = 0, family = "Arial", fontface = "bold", size = 3.15, colour = ink) +
  annotate("text", x = 0.55, y = 4.62, label = "Paired position versus unpaired position", hjust = 0, family = "Arial", size = 2.20, colour = muted) +
  annotate("text", x = 9.47, y = 4.96, label = "RNA06", hjust = 0, family = "Arial", fontface = "bold", size = 3.15, colour = purple) +
  annotate("text", x = 10.82, y = 4.96, label = "|  Partner Specificity", hjust = 0, family = "Arial", fontface = "bold", size = 3.15, colour = ink) +
  annotate("text", x = 9.47, y = 4.62, label = "True partner versus a matched paired nonpartner", hjust = 0, family = "Arial", size = 2.20, colour = muted) +
  draw_hairpin(left_hp, 2.25, 0.58) +
  geom_point(data = left_p, aes(x, y), shape = 21, size = 4.0, stroke = 0.55, fill = green, colour = "white") +
  geom_point(data = left_u, aes(x, y), shape = 21, size = 4.0, stroke = 0.55, fill = "#F7F9FA", colour = ink) +
  annotate("text", x = 0.70, y = 3.05, label = "Paired\nposition", hjust = 0, family = "Arial", fontface = "bold", size = 2.10, colour = green) +
  annotate("segment", x = 1.62, xend = left_p$x - 0.10, y = 3.03, yend = left_p$y + 0.02, colour = green, linewidth = 0.48, arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "closed")) +
  annotate("text", x = 3.95, y = 3.92, label = "Unpaired\nposition", hjust = 0, family = "Arial", fontface = "bold", size = 2.10, colour = ink) +
  annotate("segment", x = 3.85, xend = left_u$x + 0.08, y = 3.82, yend = left_u$y - 0.03, colour = muted, linewidth = 0.45, arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "closed")) +
  rounded_box(5.35, 8.42, 3.08, 4.28, "#F7FBF8", "#B7D6C5", 2.0, 0.34) +
  rounded_box(5.35, 8.42, 1.64, 2.84, "#FAFBFC", "#CCD5DC", 2.0, 0.34) +
  annotate("text", x = 6.88, y = 4.03, label = "Paired position", family = "Arial", fontface = "bold", size = 2.05, colour = green) +
  annotate("segment", x = 5.85, xend = 7.90, y = 3.45, yend = 3.45, colour = backbone, linewidth = 0.56) +
  annotate("point", x = c(6.10, 6.48, 6.86, 7.24, 7.62), y = 3.45, shape = 21, size = 2.45, stroke = 0.36, fill = "#F7F9FA", colour = backbone) +
  annotate("point", x = 6.86, y = 3.45, shape = 21, size = 3.50, stroke = 0.45, fill = green, colour = "white") +
  annotate("segment", x = 6.86, xend = 6.86, y = 3.28, yend = 3.11, colour = pair_col, linewidth = 0.42, linetype = "22") +
  annotate("text", x = 6.88, y = 2.58, label = "Unpaired position", family = "Arial", fontface = "bold", size = 2.05, colour = ink) +
  annotate("path", x = c(5.95, 6.25, 6.55, 6.88, 7.21, 7.51, 7.81), y = c(1.94, 2.10, 2.20, 2.24, 2.20, 2.10, 1.94), colour = backbone, linewidth = 0.56) +
  annotate("point", x = c(5.95, 6.25, 6.55, 6.88, 7.21, 7.51, 7.81), y = c(1.94, 2.10, 2.20, 2.24, 2.20, 2.10, 1.94), shape = 21, size = 2.30, stroke = 0.34, fill = "#F7F9FA", colour = backbone) +
  annotate("point", x = 6.88, y = 2.24, shape = 21, size = 3.35, stroke = 0.42, fill = "#F7F9FA", colour = ink) +
  geom_path(data = rna06_paths, aes(x, y, group = id), colour = backbone,
            linewidth = 0.52, lineend = "round") +
  geom_segment(data = rna06_pairs, aes(x, y, xend = xend, yend = yend),
               colour = pair_col, linewidth = 0.36, linetype = "22") +
  geom_point(data = rna06_paths, aes(x, y), shape = 21, size = 1.90,
             stroke = 0.34, fill = "#F7F9FA", colour = backbone) +
  geom_point(data = anchor, aes(x, y), shape = 21, size = 3.75, stroke = 0.55, fill = green, colour = "white") +
  geom_point(data = true_partner, aes(x, y), shape = 21, size = 3.75, stroke = 0.55, fill = red, colour = "white") +
  geom_point(data = nonpartner, aes(x, y), shape = 21, size = 3.75, stroke = 0.55, fill = purple, colour = "white") +
  annotate("text", x = anchor$x, y = 1.54, label = "A\nanchor", hjust = 0.5,
           family = "Arial", fontface = "bold", size = 1.35, colour = green) +
  annotate("text", x = true_partner$x, y = 1.54, label = "P*\ntrue partner", hjust = 0.5,
           family = "Arial", fontface = "bold", size = 1.35, colour = red) +
  annotate("text", x = nonpartner$x, y = 1.54, label = "N\npaired elsewhere", hjust = 0.5,
           family = "Arial", fontface = "bold", size = 1.35, colour = purple) +
  annotate("segment", x = partner_hp$lx + 0.10, xend = partner_hp$rx - 0.10,
           y = pair_y06, yend = pair_y06, colour = red, linewidth = 0.62,
           linetype = "22") +
  annotate("text", x = 11.35, y = pair_y06 + 0.38, label = "true pair",
           family = "Arial", fontface = "bold", size = 1.55, colour = red) +
  annotate("segment", x = nonpartner_hp$lx + 0.10, xend = nonpartner_hp$rx - 0.10,
           y = nonpartner$y, yend = nonpartner$y, colour = purple, linewidth = 0.62,
           linetype = "22") +
  annotate("text", x = 13.50, y = 4.20, label = "All candidate positions are paired",
           family = "Arial", fontface = "bold", size = 1.95, colour = purple) +
  rounded_box(0.48, 8.60, 0.48, 1.24, "#EAF4FB", "#A7CDE8", 2.2, 0.34) +
  annotate("text", x = 4.54, y = 0.86, label = "Is pairing state accessible from the representation?", family = "Arial", fontface = "bold", size = 2.22, colour = ink) +
  rounded_box(9.40, 17.52, 0.48, 1.24, "#F1EBF8", "#C5B0DE", 2.2, 0.34) +
  annotate("text", x = 13.46, y = 0.86, label = "Is exact partner identity accessible beyond pairing state?", family = "Arial", fontface = "bold", size = 2.15, colour = ink) +
  coord_fixed(xlim = c(0, 18), ylim = c(0.08, 5.55), clip = "off") +
  theme_void(base_family = "Arial") + theme(plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(4, 5, 4, 5))

width_mm <- 183; height_mm <- 61
stem <- file.path(fig_dir, "Fig3_PanelA_RNA05_RNA06_relations")
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4, family = "Arial", bg = "white"); print(p); dev.off()
svglite::svglite(paste0(stem, ".svg"), width = width_mm / 25.4, height = height_mm / 25.4, bg = "white"); print(p); dev.off()
ragg::agg_png(paste0(stem, ".png"), width = width_mm / 25.4, height = height_mm / 25.4, units = "in", res = 450, background = "white"); print(p); dev.off()
