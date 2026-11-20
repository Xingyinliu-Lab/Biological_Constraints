#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: draw_fig3_panelB_rna05_rna06_violin.R <paper_cases_rna> <data_dir> <fig_dir>")
}

source_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
data_dir <- args[[2]]
fig_dir <- args[[3]]
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

models <- data.frame(
  directory = c("ernie-rna", "gb-rna", "orthrus-4-track", "rinalmo-giga", "rna-fm"),
  model = c("ERNIE-RNA", "GB-RNA", "Orthrus-4-track", "RiNALMo-Giga", "RNA-FM"),
  stringsAsFactors = FALSE
)

constraints <- data.frame(
  directory = c("RNA05", "RNA06"),
  constraint = c("RNA05", "RNA06"),
  label = c("Pairing State", "Partner Specificity"),
  offset = c(-0.18, 0.18),
  stringsAsFactors = FALSE
)

read_one <- function(constraint_dir, constraint, constraint_label,
                     model_dir, model_label, model_index, offset) {
  output_dir <- file.path(source_root, constraint_dir, "models", model_dir, "output")
  ledger_path <- file.path(output_dir, "CLUSTER_SCORE_LEDGER.tsv")
  result_path <- file.path(output_dir, "FORMAL_RESULT.json")
  ledger <- read.delim(ledger_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"cluster_F" %in% names(ledger)) stop("Missing cluster_F: ", ledger_path)
  cluster_id_field <- setdiff(names(ledger), "cluster_F")
  if (length(cluster_id_field) != 1L) stop("Unexpected cluster ledger schema: ", ledger_path)
  if (nrow(ledger) != 20L || anyNA(ledger$cluster_F) ||
      any(ledger$cluster_F < -1 | ledger$cluster_F > 1)) {
    stop("Unexpected formal cluster population: ", ledger_path)
  }
  result <- jsonlite::fromJSON(result_path, simplifyVector = TRUE)
  f_ci <- as.numeric(result$F_95_CI)
  if (length(f_ci) != 2L) stop("Missing formal F confidence interval: ", result_path)
  data.frame(
    constraint = constraint,
    constraint_label = constraint_label,
    model = model_label,
    cluster_id = as.character(ledger[[cluster_id_field]]),
    cluster_F = as.numeric(ledger$cluster_F),
    formal_F = as.numeric(result$F),
    F_ci_low = f_ci[[1]],
    F_ci_high = f_ci[[2]],
    model_index = model_index,
    x = model_index + offset,
    violin_colour = ifelse(constraint == "RNA05", "#5F91B2", "#876FA2"),
    violin_fill = ifelse(
      constraint == "RNA05", "#C9DDEA", "#C9BED8"
    ),
    stringsAsFactors = FALSE
  )
}

pieces <- list()
k <- 1L
for (i in seq_len(nrow(models))) {
  for (j in seq_len(nrow(constraints))) {
    pieces[[k]] <- read_one(
      constraints$directory[[j]], constraints$constraint[[j]], constraints$label[[j]],
      models$directory[[i]], models$model[[i]], i,
      constraints$offset[[j]]
    )
    k <- k + 1L
  }
}
scores <- do.call(rbind, pieces)
rownames(scores) <- NULL
if (nrow(scores) != 200L || anyDuplicated(scores[c("constraint", "model", "cluster_id")])) {
  stop("Expected 200 unique model-by-constraint cluster scores")
}

rna05_medians <- aggregate(
  cluster_F ~ model,
  data = scores[scores$constraint == "RNA05", , drop = FALSE],
  FUN = median
)
model_order <- rna05_medians$model[order(-rna05_medians$cluster_F, rna05_medians$model)]
scores$model_index <- match(scores$model, model_order)
scores$x <- scores$model_index + ifelse(scores$constraint == "RNA05", -0.18, 0.18)

source_data <- file.path(data_dir, "RNA05_RNA06_cluster_scores.tsv")
write.table(scores, source_data, sep = "\t", row.names = FALSE, quote = FALSE)

medians <- aggregate(cluster_F ~ model + model_index + constraint + x,
                     data = scores, FUN = median)
medians$model <- factor(medians$model, levels = model_order)
medians$constraint <- factor(medians$constraint, levels = constraints$constraint)
medians <- medians[order(medians$model, medians$constraint), ]

scores$model <- factor(scores$model, levels = model_order)
scores$constraint <- factor(scores$constraint, levels = constraints$constraint)

violin_width <- 0.32
p <- ggplot() +
  annotate("rect", xmin = 0.56, xmax = 0.73, ymin = 1.045, ymax = 1.095,
           fill = "#C9DDEA", colour = "#5F91B2", linewidth = 0.28) +
  annotate("text", x = 0.78, y = 1.07, label = "RNA05",
           hjust = 0, vjust = 0.5, family = "Arial", fontface = "bold",
           size = 2.15, colour = "#456D86") +
  annotate("rect", xmin = 1.65, xmax = 1.82, ymin = 1.045, ymax = 1.095,
           fill = "#C9BED8", colour = "#876FA2", linewidth = 0.28) +
  annotate("text", x = 1.87, y = 1.07, label = "RNA06",
           hjust = 0, vjust = 0.5, family = "Arial", fontface = "bold",
           size = 2.15, colour = "#5C4C70") +
  geom_hline(yintercept = 0, colour = "#68737B", linewidth = 0.34,
             linetype = "32") +
  geom_violin(
    data = scores,
    aes(x = x, y = cluster_F, group = interaction(model, constraint),
        fill = violin_fill, colour = violin_colour),
    width = violin_width, scale = "width", trim = TRUE, adjust = 0.82,
    linewidth = 0.38
  ) +
  geom_line(
    data = medians,
    aes(x = x, y = cluster_F, group = model),
    colour = "#161616", linewidth = 0.42, lineend = "round"
  ) +
  geom_point(
    data = medians,
    aes(x = x, y = cluster_F),
    shape = 16, size = 1.05, colour = "#111111"
  ) +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_x_continuous(
    breaks = seq_len(nrow(models)), labels = model_order,
    limits = c(0.55, 5.45), expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(-0.62, 1.10), breaks = c(-0.5, 0, 0.5, 1),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(x = NULL, y = "Dependency-cluster F score") +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 7.2, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.34, colour = "#35434B"),
    axis.ticks = element_line(linewidth = 0.30, colour = "#35434B"),
    axis.text = element_text(colour = "#25313A"),
    axis.text.x = element_text(size = 6.6, face = "bold", angle = 20,
                               hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 6.4),
    axis.title.y = element_text(size = 7.0, margin = margin(r = 5)),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(8, 9, 7, 7)
  )

save_plot <- function(plot, stem, width_mm = 80, height_mm = 84) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  grDevices::cairo_pdf(file.path(fig_dir, paste0(stem, ".pdf")), width = w,
                       height = h, family = "Arial", bg = "white", onefile = FALSE)
  print(plot); dev.off()
  svglite::svglite(file.path(fig_dir, paste0(stem, ".svg")), width = w,
                   height = h, bg = "white")
  print(plot); dev.off()
  ragg::agg_png(file.path(fig_dir, paste0(stem, ".png")), width = w,
                height = h, units = "in", res = 450, background = "white")
  print(plot); dev.off()
}

save_plot(p, "Fig3_PanelB_RNA05_RNA06_violin")
cat(source_data, "\n", sep = "")
