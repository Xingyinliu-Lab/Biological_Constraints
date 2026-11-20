#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: draw_fig3_panelC_shared_outer_clusters.R <paper_cases_rna> <data_dir> <fig_dir>")
}

source_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
data_dir <- args[[2]]
fig_dir <- args[[3]]
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

models <- data.frame(
  directory = c("gb-rna", "ernie-rna", "rinalmo-giga", "rna-fm", "orthrus-4-track"),
  model = c("GB-RNA", "ERNIE-RNA", "RiNALMo-Giga", "RNA-FM", "Orthrus-4-track"),
  stringsAsFactors = FALSE
)

read_cluster_scores <- function(constraint, model_dir) {
  path <- file.path(source_root, constraint, "models", model_dir, "output", "CLUSTER_SCORE_LEDGER.tsv")
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"cluster_F" %in% names(x) || ncol(x) != 2L || nrow(x) != 20L) {
    stop("Expected 20-row F-only cluster ledger: ", path)
  }
  id_col <- setdiff(names(x), "cluster_F")
  data.frame(
    shared_outer_cluster = sub("^.*:", "", as.character(x[[id_col]])),
    cluster_F = as.numeric(x$cluster_F),
    stringsAsFactors = FALSE
  )
}

all_scores <- lapply(seq_len(nrow(models)), function(i) {
  rna05 <- read_cluster_scores("RNA05", models$directory[[i]])
  rna06 <- read_cluster_scores("RNA06", models$directory[[i]])
  names(rna05)[2] <- "RNA05_cluster_F"
  names(rna06)[2] <- "RNA06_cluster_F"
  joined <- merge(rna05, rna06, by = "shared_outer_cluster", all = FALSE, sort = TRUE)
  if (nrow(joined) != 20L || anyNA(joined)) {
    stop("Expected exactly 20 shared outer clusters for ", models$model[[i]])
  }
  joined$model <- models$model[[i]]
  joined
})
all_scores <- do.call(rbind, all_scores)
rownames(all_scores) <- NULL
write.table(all_scores, file.path(data_dir, "RNA05_RNA06_shared_outer_cluster_F.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

correlations <- do.call(rbind, lapply(models$model, function(model_name) {
  d <- all_scores[all_scores$model == model_name, , drop = FALSE]
  test <- suppressWarnings(cor.test(d$RNA05_cluster_F, d$RNA06_cluster_F,
                                    method = "spearman", exact = FALSE))
  data.frame(model = model_name, n_shared_outer_clusters = nrow(d),
             spearman_rho = unname(test$estimate), p_value = test$p.value,
             stringsAsFactors = FALSE)
}))
write.table(correlations, file.path(data_dir, "RNA05_RNA06_shared_outer_cluster_spearman.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

save_plot <- function(plot, stem, width_mm = 72, height_mm = 72) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  grDevices::cairo_pdf(file.path(fig_dir, paste0(stem, ".pdf")), width = w, height = h,
                       family = "Arial", bg = "white", onefile = FALSE)
  print(plot); dev.off()
  svglite::svglite(file.path(fig_dir, paste0(stem, ".svg")), width = w, height = h, bg = "white")
  print(plot); dev.off()
  ragg::agg_png(file.path(fig_dir, paste0(stem, ".png")), width = w, height = h,
                units = "in", res = 450, background = "white")
  print(plot); dev.off()
}

for (i in seq_len(nrow(models))) {
  d <- all_scores[all_scores$model == models$model[[i]], , drop = FALSE]
  stat <- correlations[correlations$model == models$model[[i]], , drop = FALSE]
  p_label <- if (stat$p_value < 0.001) "P < 0.001" else sprintf("P = %.3f", stat$p_value)
  rho_label <- sprintf("Spearman %s = %.2f\n%s", "\u03c1", stat$spearman_rho, p_label)
  p <- ggplot(d, aes(x = RNA05_cluster_F, y = RNA06_cluster_F)) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.38, colour = "#9AA6AE") +
    geom_hline(yintercept = 0, linewidth = 0.28, linetype = "32", colour = "#C1C8CD") +
    geom_vline(xintercept = 0, linewidth = 0.28, linetype = "32", colour = "#C1C8CD") +
    geom_point(shape = 21, size = 2.45, stroke = 0.45,
               fill = "#6E8FA5", colour = "#28343B", alpha = 0.84) +
    annotate("text", x = -0.54, y = 0.97, label = rho_label, hjust = 0, vjust = 1,
             family = "Arial", size = 2.55, colour = "#39464E") +
    labs(title = models$model[[i]], x = "RNA05 cluster F", y = "RNA06 cluster F") +
    scale_x_continuous(limits = c(-0.6, 1.02), breaks = c(-0.5, 0, 0.5, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.6, 1.02), breaks = c(-0.5, 0, 0.5, 1), expand = c(0, 0)) +
    coord_equal(clip = "off") +
    theme_bw(base_size = 7.0, base_family = "Arial") +
    theme(
      plot.title = element_text(size = 8.2, face = "bold", hjust = 0.5, colour = "#24313A", margin = margin(b = 3)),
      axis.title = element_text(size = 7.0, colour = "#25313A"),
      axis.text = element_text(size = 6.2, colour = "#25313A"),
      axis.line = element_line(linewidth = 0.36, colour = "#35434B"),
      axis.ticks = element_line(linewidth = 0.32, colour = "#35434B"),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(5, 5, 5, 6)
    )
  save_plot(p, paste0("Fig3_PanelC_RNA05_RNA06_shared_outer_clusters_", models$directory[[i]]))
}
