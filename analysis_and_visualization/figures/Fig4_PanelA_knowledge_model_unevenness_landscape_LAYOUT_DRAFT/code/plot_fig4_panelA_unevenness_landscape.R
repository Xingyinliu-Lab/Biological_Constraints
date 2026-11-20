#!/usr/bin/env Rscript

# Figure 4A: formal-score dispersion versus mean formal score.
# Uses released F results only; it does not invoke a scorer or forward model.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: plot_fig4_panelA_unevenness_landscape.R <performance_tsv> <data_dir> <fig_dir>")
}

input_path <- normalizePath(args[[1]], mustWork = TRUE)
data_dir <- normalizePath(args[[2]], mustWork = TRUE)
fig_dir <- normalizePath(args[[3]], mustWork = TRUE)

scores <- read.delim(input_path, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("Constraint", "Short_name", "Model", "F")
if (length(setdiff(required, names(scores)))) stop("Performance table is missing required formal-F fields")
scores$F <- as.numeric(scores$F)
if (any(!is.finite(scores$F))) stop("Non-finite F values in performance table")

constraint_order <- c(
  "DNA01", "DNA02-02", "DNA03", "DNA06",
  "RNA01-02", "RNA04", "RNA05", "RNA06", "RNA09", "RNA10",
  "PS01", "PS02A", "PS02B", "PS04", "PS05", "PS08", "PS09", "PS10",
  "ST01", "ST02", "ST04", "ST05"
)
domain_lookup <- c(
  DNA01 = "DNA", `DNA02-02` = "DNA", DNA03 = "DNA", DNA06 = "DNA",
  `RNA01-02` = "RNA", RNA04 = "RNA", RNA05 = "RNA", RNA06 = "RNA", RNA09 = "RNA", RNA10 = "RNA",
  PS01 = "Protein sequence", PS02A = "Protein sequence", PS02B = "Protein sequence", PS04 = "Protein sequence",
  PS05 = "Protein sequence", PS08 = "Protein sequence", PS09 = "Protein sequence", PS10 = "Protein sequence",
  ST01 = "Protein structure", ST02 = "Protein structure", ST04 = "Protein structure", ST05 = "Protein structure"
)
domain_order <- c("DNA", "RNA", "Protein sequence", "Protein structure")

scores <- scores[scores$Constraint %in% constraint_order, ]
counts <- table(scores$Constraint)
if (!all(constraint_order %in% names(counts)) || any(counts[constraint_order] != 5L)) {
  stop("Exactly five formal model scores are required for each constraint")
}

rows <- lapply(constraint_order, function(id) {
  dat <- scores[scores$Constraint == id, c("Model", "F", "Short_name")]
  f <- setNames(dat$F, dat$Model)
  domain <- unname(domain_lookup[[id]])
  data.frame(
    constraint = id,
    short_name = dat$Short_name[[1]],
    domain = domain,
    model_dispersion_F = max(f) - min(f),
    mean_F = mean(f),
    minimum_F = min(f),
    maximum_F = max(f),
    minimum_model = names(f)[which.min(f)],
    maximum_model = names(f)[which.max(f)],
    stringsAsFactors = FALSE
  )
})
landscape <- do.call(rbind, rows)
landscape$domain <- factor(landscape$domain, levels = domain_order)
landscape$label <- ifelse(
  landscape$model_dispersion_F >= 0.70 | landscape$mean_F <= 0.20 | landscape$mean_F >= 0.90,
  landscape$short_name,
  ""
)

write.table(landscape, file = file.path(data_dir, "Fig4_PanelA_knowledge_model_unevenness_landscape.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

domain_colours <- c(DNA = "#4E7FA3", RNA = "#4C997F", `Protein sequence` = "#BA8547", `Protein structure` = "#8066A0")

plot <- ggplot(landscape, aes(x = model_dispersion_F, y = mean_F, colour = domain)) +
  geom_point(size = 3.3, stroke = 0.3, alpha = 0.96) +
  ggrepel::geom_text_repel(
    aes(label = label), show.legend = FALSE, colour = "#1A252C", size = 2.65,
    min.segment.length = 0, box.padding = 0.3, point.padding = 0.2,
    seed = 20261022, max.overlaps = Inf, segment.colour = "#9EA7AC", segment.size = 0.28
  ) +
  scale_colour_manual(values = domain_colours, name = NULL) +
  scale_x_continuous(limits = c(0, 1.02), breaks = c(0, 0.25, 0.5, 0.75, 1.0), expand = expansion(mult = c(0.02, 0.06))) +
  scale_y_continuous(limits = c(-0.05, 1.02), breaks = c(0, 0.25, 0.5, 0.75, 1.0), expand = expansion(mult = c(0.02, 0.06))) +
  labs(
    title = "A | Knowledge accessibility and model dispersion",
    x = expression("Model dispersion, " * Delta * F),
    y = expression("Mean FM score, " * bar(F))
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 8.5, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "#30383E", fill = NA, linewidth = 0.4),
    axis.line = element_line(colour = "#1D252B", linewidth = 0.35),
    axis.ticks = element_line(colour = "#1D252B", linewidth = 0.35),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text = element_text(colour = "#182027", size = 7.0),
    axis.title = element_text(colour = "#182027", size = 8.4),
    plot.title = element_text(colour = "#111D25", size = 9.2, face = "bold", hjust = 0),
    legend.position = "top",
    legend.text = element_text(size = 7.0, colour = "#182027"),
    legend.key.width = grid::unit(5, "mm"),
    plot.margin = margin(5, 10, 5, 5, unit = "mm")
  )

stem <- file.path(fig_dir, "Fig4_PanelA_knowledge_model_unevenness_landscape_LAYOUT_DRAFT")
width_in <- 150 / 25.4
height_in <- 115 / 25.4
ggsave(paste0(stem, ".pdf"), plot, width = width_in, height = height_in, device = grDevices::cairo_pdf, dpi = 600)
svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
print(plot)
grDevices::dev.off()
ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in, units = "in", res = 500)
print(plot)
grDevices::dev.off()

