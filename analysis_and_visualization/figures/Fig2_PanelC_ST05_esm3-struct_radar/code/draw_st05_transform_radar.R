#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: draw_st05_transform_radar.R <paper_cases_ST05> <data_dir> <fig_dir>")
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

all_categories <- c(
  "Identity", "Rotation 45", "Rotation 90", "Rotation 180",
  "Translation", "General combined"
)
expected_counts <- c(60L, 180L, 180L, 180L, 180L, 180L)
categories <- all_categories[all_categories != "Identity"]

classify_transform <- function(transform_id) {
  ifelse(
    transform_id == "IDENTITY", "Identity",
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
  )
}

read_model <- function(directory, model) {
  path <- file.path(
    source_root, "models", directory, "output", "ST05_TRANSFORM_RESULT_LEDGER.tsv"
  )
  ledger <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("base_object_id", "transform_id", "F_credit")
  if (!all(required %in% names(ledger))) {
    stop("Missing required ledger fields for ", model)
  }
  if (nrow(ledger) != 960L || length(unique(ledger$base_object_id)) != 60L) {
    stop("Unexpected formal ST05 population for ", model)
  }
  if (anyNA(ledger$F_credit) || any(!ledger$F_credit %in% c(0, 1))) {
    stop("Invalid ST05 credit values for ", model)
  }

  ledger$category <- classify_transform(ledger$transform_id)
  if (anyNA(ledger$category)) {
    stop("Unclassified transform ID for ", model)
  }

  counts <- table(factor(ledger$category, levels = all_categories))
  if (!identical(as.integer(counts), expected_counts)) {
    stop("Unexpected transform-category population for ", model)
  }

  out <- aggregate(F_credit ~ category, data = ledger, FUN = mean)
  names(out)[names(out) == "F_credit"] <- "F_satisfaction_fraction"
  out$model <- model
  out$n_relations <- as.integer(counts[match(out$category, names(counts))])
  out
}

radar_data <- do.call(rbind, Map(read_model, models$directory, models$model))
radar_data <- radar_data[radar_data$category != "Identity", ]
radar_data$category <- factor(radar_data$category, levels = categories)
radar_data <- radar_data[order(radar_data$model, radar_data$category), ]
rownames(radar_data) <- NULL

write.table(
  radar_data,
  file.path(data_dir, "ST05_transform_category_satisfaction.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

angles <- pi / 2 - 2 * pi * (seq_along(categories) - 1) / length(categories)
axis_data <- data.frame(
  category = factor(categories, levels = categories),
  angle = angles,
  x = cos(angles),
  y = sin(angles),
  label_x = 1.22 * cos(angles),
  label_y = 1.22 * sin(angles),
  stringsAsFactors = FALSE
)

ring_levels <- c(0.25, 0.50, 0.75, 1.00)
ring_data <- do.call(rbind, lapply(ring_levels, function(level) {
  data.frame(
    level = factor(level),
    order = seq_along(categories),
    x = level * cos(angles),
    y = level * sin(angles)
  )
}))

make_radar <- function(model_name, colour) {
  values <- radar_data[radar_data$model == model_name, ]
  values <- values[match(categories, as.character(values$category)), ]
  values$angle <- angles
  values$x <- values$F_satisfaction_fraction * cos(angles)
  values$y <- values$F_satisfaction_fraction * sin(angles)
  axis_plot <- axis_data
  category_labels <- c(
    "Rotation 45", "Rotation 90", "Rotation 180", "Translation",
    "General combined"
  )
  axis_plot$display_label <- paste0(
    category_labels, "\n", sprintf("%.2f", values$F_satisfaction_fraction)
  )

  ggplot() +
    geom_polygon(
      data = ring_data,
      aes(x = x, y = y, group = level),
      fill = NA, colour = "#D8DDE0", linewidth = 0.30
    ) +
    geom_segment(
      data = axis_data,
      aes(x = 0, y = 0, xend = x, yend = y),
      colour = "#D8DDE0", linewidth = 0.30
    ) +
    geom_polygon(
      data = values,
      aes(x = x, y = y, group = 1),
      fill = colour, colour = colour, alpha = 0.20, linewidth = 0.75
    ) +
    geom_point(
      data = values,
      aes(x = x, y = y),
      shape = 21, size = 2.2, stroke = 0.50,
      fill = "white", colour = colour
    ) +
    geom_text(
      data = axis_plot,
      aes(x = label_x, y = label_y, label = display_label),
      family = "Arial", fontface = "bold", size = 2.35,
      lineheight = 0.92, colour = "#25313A"
    ) +
    coord_equal(xlim = c(-1.38, 1.38), ylim = c(-1.30, 1.35), clip = "off") +
    labs(
      title = model_name,
      subtitle = "Satisfaction fraction by transformation class"
    ) +
    theme_void(base_family = "Arial", base_size = 7.2) +
    theme(
      plot.title = element_text(
        size = 9.0, face = "bold", colour = "#25313A", hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 6.1, colour = "#68757D", hjust = 0.5, margin = margin(b = 5)
      ),
      plot.margin = margin(8, 12, 8, 12),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

save_radar <- function(plot, stem, width_mm = 92, height_mm = 88) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  grDevices::cairo_pdf(
    file.path(fig_dir, paste0(stem, ".pdf")),
    width = width_in, height = height_in,
    family = "Arial", bg = "white", onefile = FALSE
  )
  print(plot)
  dev.off()

  ragg::agg_png(
    file.path(fig_dir, paste0(stem, ".png")),
    width = width_in, height = height_in,
    units = "in", res = 450, background = "white"
  )
  print(plot)
  dev.off()
}

for (i in seq_len(nrow(models))) {
  plot <- make_radar(models$model[[i]], models$colour[[i]])
  stem <- paste0("Fig2_PanelC_ST05_", models$directory[[i]], "_radar")
  save_radar(plot, stem)
}

cat(file.path(data_dir, "ST05_transform_category_satisfaction.tsv"), "\n", sep = "")
