#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: draw_st05_boxplot_violin.R <paper_cases_ST05> <data_dir> <fig_dir>")
}

source_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
data_dir <- args[[2]]
fig_dir <- args[[3]]
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

models <- data.frame(
  directory = c("sslm", "proteinmpnn-ca", "esm-if1", "esm3-struct", "protoken"),
  model = c("SSLM", "ProteinMPNN-CA", "ESM-IF1", "ESM3-struct", "ProToken"),
  colour = c("#3F6F8F", "#4E8B72", "#7B67A8", "#C47B50", "#67717E"),
  stringsAsFactors = FALSE
)

read_model <- function(directory, model) {
  output_dir <- file.path(source_root, "models", directory, "output")
  ledger_path <- file.path(output_dir, "ST05_TRANSFORM_RESULT_LEDGER.tsv")
  result_path <- file.path(output_dir, "EVALUATION_RESULT.json")

  ledger <- read.delim(ledger_path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("base_object_id", "F_credit", "status")
  if (!all(required %in% names(ledger))) {
    stop("Missing required ledger fields for ", model)
  }
  if (nrow(ledger) != 960L || length(unique(ledger$base_object_id)) != 60L) {
    stop("Unexpected formal ST05 population for ", model)
  }
  if (anyNA(ledger$F_credit) || any(!ledger$F_credit %in% c(0, 1))) {
    stop("Invalid ST05 credit values for ", model)
  }

  structure_scores <- aggregate(
    F_credit ~ base_object_id,
    data = ledger,
    FUN = mean
  )
  names(structure_scores)[names(structure_scores) == "F_credit"] <- "structure_F"

  result <- jsonlite::fromJSON(result_path, simplifyVector = TRUE)
  ci <- as.numeric(result$F_95_CI)
  structure_scores$model <- model
  structure_scores$formal_F <- as.numeric(result$F)
  structure_scores$ci_low <- ci[[1]]
  structure_scores$ci_high <- ci[[2]]
  structure_scores
}

scores <- do.call(
  rbind,
  Map(read_model, models$directory, models$model)
)
rownames(scores) <- NULL

display_levels <- models$model
scores$model <- factor(scores$model, levels = display_levels)
if (nrow(scores) != 300L || anyNA(scores$structure_F)) {
  stop("Expected 300 complete structure-level scores")
}
if (any(scores$structure_F < 0 | scores$structure_F > 1)) {
  stop("Structure-level scores fall outside [0, 1]")
}

source_data <- file.path(data_dir, "ST05_structure_level_scores.tsv")
write.table(scores, source_data, sep = "\t", row.names = FALSE, quote = FALSE)

summary_data <- unique(scores[c("model", "formal_F", "ci_low", "ci_high")])
summary_data <- summary_data[order(summary_data$model), ]
if (nrow(summary_data) != 5L) {
  stop("Expected one formal result per model")
}

palette <- setNames(models$colour, models$model)

base_theme <- theme_bw(base_size = 7.2, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.34, colour = "#35434B"),
    axis.ticks = element_line(linewidth = 0.30, colour = "#35434B"),
    axis.text = element_text(colour = "#25313A"),
    axis.text.y = element_text(size = 6.4),
    axis.text.x = element_text(
      size = 6.6, face = "bold", angle = 25, hjust = 1, vjust = 1
    ),
    axis.title.x = element_text(size = 7.0, margin = margin(t = 5)),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    aspect.ratio = 0.82,
    plot.title = element_text(size = 8.2, face = "bold", colour = "#25313A"),
    plot.subtitle = element_text(size = 5.9, colour = "#68757D", margin = margin(b = 6)),
    plot.margin = margin(7, 10, 7, 7)
  )

common_layers <- list(
  scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.015, 0.025))
  ),
  scale_fill_manual(values = palette),
  coord_cartesian(clip = "off"),
  labs(
    title = "Overall rigid-body invariance",
    subtitle = "n = 60 base structures per model",
    x = NULL,
    y = "Structure-level formal score"
  ),
  base_theme
)

box_plot <- ggplot(scores, aes(x = model, y = structure_F, fill = model)) +
  geom_boxplot(
    aes(colour = model),
    width = 0.54,
    linewidth = 0.40,
    outlier.shape = NA,
    alpha = 0.70
  ) +
  geom_point(
    aes(colour = model),
    position = position_jitter(width = 0.13, height = 0, seed = 20261021),
    size = 0.55,
    alpha = 0.22
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 16,
    size = 0.70,
    colour = "#111111"
  ) +
  scale_colour_manual(values = palette) +
  common_layers

violin_plot <- ggplot(scores, aes(x = model, y = structure_F, fill = model)) +
  geom_violin(
    aes(colour = model),
    width = 0.82,
    scale = "width",
    trim = TRUE,
    adjust = 0.78,
    linewidth = 0.40,
    alpha = 0.58
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 16,
    size = 0.70,
    colour = "#111111"
  ) +
  scale_colour_manual(values = palette) +
  common_layers

save_plot <- function(plot, stem, width_mm = 112, height_mm = 95) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  grDevices::cairo_pdf(
    file.path(fig_dir, paste0(stem, ".pdf")),
    width = width_in,
    height = height_in,
    family = "Arial",
    bg = "white",
    onefile = FALSE
  )
  print(plot)
  dev.off()

  svglite::svglite(
    file.path(fig_dir, paste0(stem, ".svg")),
    width = width_in,
    height = height_in,
    bg = "white"
  )
  print(plot)
  dev.off()

  ragg::agg_png(
    file.path(fig_dir, paste0(stem, ".png")),
    width = width_in,
    height = height_in,
    units = "in",
    res = 450,
    background = "white"
  )
  print(plot)
  dev.off()
}

save_plot(box_plot, "Fig2_PanelB_ST05_boxplot")
save_plot(violin_plot, "Fig2_PanelB_ST05_violin")

cat(source_data, "\n", sep = "")
