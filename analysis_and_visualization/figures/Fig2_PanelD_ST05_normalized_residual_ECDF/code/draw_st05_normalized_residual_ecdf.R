#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: draw_st05_normalized_residual_ecdf.R <paper_cases_ST05> <data_dir> <fig_dir>")
}

source_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
data_dir <- args[[2]]
fig_dir <- args[[3]]
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

models <- data.frame(
  directory = c("sslm", "proteinmpnn-ca", "esm-if1", "esm3-struct", "protoken"),
  model = c("SSLM", "ProteinMPNN-CA", "ESM-IF1", "ESM3-struct", "ProToken"),
  stringsAsFactors = FALSE
)

categories <- c(
  "Rotation 45", "Rotation 90", "Rotation 180",
  "Translation", "General combined"
)

palette <- c(
  "Rotation 45" = "#4E79A7",
  "Rotation 90" = "#4E9678",
  "Rotation 180" = "#8064A2",
  "Translation" = "#D58245",
  "General combined" = "#C45D5D"
)

classify_transform <- function(transform_id) {
  ifelse(
    grepl("^ROT_[XYZ]_45$", transform_id), "Rotation 45",
    ifelse(
      grepl("^ROT_[XYZ]_90$", transform_id), "Rotation 90",
      ifelse(
        grepl("^ROT_[XYZ]_180$", transform_id), "Rotation 180",
        ifelse(
          grepl("^TRANS_", transform_id), "Translation",
          ifelse(grepl("^COMBO_", transform_id), "General combined", NA_character_)
        )
      )
    )
  )
}

read_model <- function(directory, model) {
  path <- file.path(
    source_root, "models", directory, "output", "ST05_TRANSFORM_RESULT_LEDGER.tsv"
  )
  ledger <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("relation_id", "transform_id", "normalized_residual", "F_credit")
  if (!all(required %in% names(ledger))) {
    stop("Missing required ledger fields for ", model)
  }

  ledger <- ledger[ledger$transform_id != "IDENTITY", required]
  if (nrow(ledger) != 900L) {
    stop("Expected 900 non-identity relations for ", model)
  }
  if (anyNA(ledger$normalized_residual) || any(ledger$normalized_residual < 0)) {
    stop("Invalid normalized residual for ", model)
  }
  if (anyNA(ledger$F_credit) || any(!ledger$F_credit %in% c(0, 1))) {
    stop("Invalid credit for ", model)
  }

  ledger$category <- classify_transform(ledger$transform_id)
  if (anyNA(ledger$category)) {
    stop("Unclassified transform for ", model)
  }
  counts <- table(factor(ledger$category, levels = categories))
  if (!all(as.integer(counts) == 180L)) {
    stop("Expected 180 relations in every transformation class for ", model)
  }

  ledger$model <- model
  ledger[, c(
    "model", "relation_id", "transform_id", "category",
    "normalized_residual", "F_credit"
  )]
}

raw_data <- do.call(rbind, Map(read_model, models$directory, models$model))
rownames(raw_data) <- NULL
raw_data$model <- factor(raw_data$model, levels = models$model)
raw_data$category <- factor(raw_data$category, levels = categories)

zero_display_x <- -5
raw_data$exact_zero <- raw_data$normalized_residual == 0
raw_data$log10_normalized_residual <- ifelse(
  raw_data$exact_zero,
  zero_display_x,
  log10(raw_data$normalized_residual)
)

source_path <- file.path(data_dir, "ST05_normalized_residual_ecdf_source.tsv")
write.table(raw_data, source_path, sep = "\t", row.names = FALSE, quote = FALSE)

make_ecdf <- function(x, model, category) {
  x <- sort(x)
  n <- length(x)
  data.frame(
    model = model,
    category = category,
    x = c(zero_display_x - 0.08, x),
    cumulative_fraction = c(0, seq_len(n) / n),
    stringsAsFactors = FALSE
  )
}

split_key <- interaction(raw_data$model, raw_data$category, drop = TRUE, lex.order = TRUE)
ecdf_data <- do.call(rbind, lapply(split(raw_data, split_key), function(d) {
  make_ecdf(
    d$log10_normalized_residual,
    as.character(d$model[[1]]),
    as.character(d$category[[1]])
  )
}))
rownames(ecdf_data) <- NULL
ecdf_data$model <- factor(ecdf_data$model, levels = models$model)
ecdf_data$category <- factor(ecdf_data$category, levels = categories)

threshold_points <- aggregate(
  F_credit ~ model + category,
  data = raw_data,
  FUN = mean
)
names(threshold_points)[names(threshold_points) == "F_credit"] <- "F_satisfaction_fraction"
threshold_points$model <- factor(threshold_points$model, levels = models$model)
threshold_points$category <- factor(threshold_points$category, levels = categories)

write.table(
  threshold_points,
  file.path(data_dir, "ST05_transform_category_threshold_points.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

plot <- ggplot(
  ecdf_data,
  aes(x = x, y = cumulative_fraction, colour = category, group = category)
) +
  geom_step(linewidth = 0.58, direction = "hv") +
  geom_vline(
    xintercept = 0,
    linetype = "22",
    linewidth = 0.42,
    colour = "#252D33"
  ) +
  geom_point(
    data = threshold_points,
    aes(x = 0, y = F_satisfaction_fraction, colour = category),
    inherit.aes = FALSE,
    shape = 16,
    size = 1.15
  ) +
  facet_wrap(~model, nrow = 1) +
  scale_colour_manual(values = palette, breaks = categories) +
  scale_x_continuous(
    limits = c(-5.10, 2.35),
    breaks = c(-5, -4, -2, 0, 2),
    labels = c("q = 0", "-4", "-2", "0", "2"),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0, 0.015))
  ) +
  labs(
    title = "Normalized residual distributions",
    x = expression(log[10](residual / tolerance)),
    y = "Cumulative fraction",
    colour = NULL
  ) +
  theme_bw(base_size = 6.5, base_family = "Arial") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "#E4E7E9", linewidth = 0.25),
    panel.border = element_rect(colour = "#59646B", linewidth = 0.35),
    axis.ticks = element_line(colour = "#354047", linewidth = 0.30),
    axis.text = element_text(colour = "#25313A", size = 5.7),
    axis.text.x = element_text(size = 5.2),
    axis.title = element_text(colour = "#25313A", size = 6.4),
    strip.background = element_blank(),
    strip.text = element_text(
      colour = "#25313A", face = "bold", size = 6.1,
      margin = margin(t = 0, b = 1.5)
    ),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 5.7, colour = "#25313A"),
    legend.key.width = grid::unit(8, "mm"),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.margin = margin(0, 0, 2, 0),
    plot.title = element_text(size = 8.2, face = "bold", colour = "#25313A"),
    plot.margin = margin(5, 6, 5, 6),
    panel.spacing.x = grid::unit(2.2, "mm")
  ) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE))

width_mm <- 183
height_mm <- 52
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
stem <- file.path(fig_dir, "Fig2_PanelD_ST05_normalized_residual_ECDF")

grDevices::cairo_pdf(
  paste0(stem, ".pdf"),
  width = width_in, height = height_in,
  family = "Arial", bg = "white", onefile = FALSE
)
print(plot)
dev.off()

svglite::svglite(
  paste0(stem, ".svg"),
  width = width_in, height = height_in,
  bg = "white"
)
print(plot)
dev.off()

ragg::agg_png(
  paste0(stem, ".png"),
  width = width_in, height = height_in,
  units = "in", res = 450, background = "white"
)
print(plot)
dev.off()

cat(source_path, "\n", sep = "")
