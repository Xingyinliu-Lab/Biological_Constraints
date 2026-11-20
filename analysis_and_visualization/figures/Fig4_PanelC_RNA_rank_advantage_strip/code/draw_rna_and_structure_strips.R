suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
  library(ragg)
})

# Figure contract: each strip shows a selected, CI-supported model-pair
# reversal. The left coordinate encodes rank within the complete five-model
# domain panel; the right coordinate encodes the signed formal F difference.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script with Rscript")
panel_dir <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), ".."), mustWork = TRUE)
fig4c_dir <- normalizePath(file.path(panel_dir, "..", "..", "analysis", "fig4c_statistical_tests"), mustWork = TRUE)
dir.create(file.path(panel_dir, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(panel_dir, "fig"), recursive = TRUE, showWarnings = FALSE)

scores <- read.delim(file.path(fig4c_dir, "data", "Supplementary_Table_S4_Complete_Model_Performance.tsv"),
                     check.names = FALSE, stringsAsFactors = FALSE)
scores$rank <- ave(-scores$F, scores$Constraint, FUN = function(x) rank(x, ties.method = "min"))

draw_case <- function(domain_name, constraints, relation_labels, model_a, model_b, tag) {
  case <- scores[scores$Constraint %in% constraints & scores$Model %in% c(model_a, model_b),
                 c("Constraint", "Short_name", "Model", "F", "F_95_CI", "rank"), drop = FALSE]
  if (nrow(case) != 2L * length(constraints)) stop("Missing selected model-pair records")
  case$row <- match(case$Constraint, constraints)
  case$pair_label <- ifelse(case$Model == model_a, "A", "B")

  delta <- do.call(rbind, lapply(constraints, function(id) {
    d <- case[case$Constraint == id, ]
    fa <- d$F[d$Model == model_a]
    fb <- d$F[d$Model == model_b]
    z <- fa - fb
    data.frame(
      Constraint = id, relation = relation_labels[[id]], row = match(id, constraints),
      F_model_A = fa, F_model_B = fb, delta_F = z, delta_label = sprintf("%+.2f", z),
      direction = ifelse(z >= 0, "A higher", "B higher"), stringsAsFactors = FALSE
    )
  }))
  slots <- expand.grid(row = seq_along(constraints), rank = 1:5)
  slots$occupied_by_pair <- mapply(function(row, rank) any(case$row == row & case$rank == rank), slots$row, slots$rank)

  write.table(case, file.path(panel_dir, "data", paste0(tag, "_source_data.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(delta, file.path(panel_dir, "data", paste0(tag, "_pair_differences.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  delta_to_x <- function(x) 9 + x / 0.8 * 2
  colors <- c("A" = "#356C7D", "B" = "#B56A3A", "A higher" = "#356C7D", "B higher" = "#B56A3A")
  relation_text <- data.frame(row = seq_along(constraints), relation = unname(relation_labels[constraints]))
  plot <- ggplot() +
    annotate("text", x = 3, y = -0.05, label = "Rank among five models",
             family = "Arial", size = 2.9, fontface = "bold", colour = "#223A4B") +
    annotate("text", x = 1, y = 0.25, label = "highest", family = "Arial", size = 2.25, hjust = 0, colour = "#66727A") +
    annotate("text", x = 5, y = 0.25, label = "lowest", family = "Arial", size = 2.25, hjust = 1, colour = "#66727A") +
    annotate("text", x = 9, y = -0.05, label = "Pairwise advantage (Delta F)",
             family = "Arial", size = 2.9, fontface = "bold", colour = "#223A4B") +
    annotate("text", x = 7.25, y = 0.25, label = model_b, family = "Arial", size = 1.95, hjust = 0, colour = colors[["B"]]) +
    annotate("text", x = 10.75, y = 0.25, label = model_a, family = "Arial", size = 1.95, hjust = 1, colour = colors[["A"]]) +
    geom_text(data = relation_text, aes(x = -4.3, y = row, label = relation), hjust = 0,
              size = 2.75, family = "Arial", colour = "#252525") +
    geom_point(data = slots[!slots$occupied_by_pair, ], aes(x = rank, y = row), shape = 21,
               size = 3.15, stroke = 0.45, fill = "#E2E6E8", colour = "#C7CED1") +
    geom_point(data = case, aes(x = rank, y = row, fill = pair_label), shape = 21,
               size = 4.35, stroke = 0.45, colour = "white") +
    geom_text(data = case, aes(x = rank, y = row, label = pair_label), family = "Arial",
              size = 2.45, fontface = "bold", colour = "white", vjust = 0.35) +
    geom_vline(xintercept = 9, linewidth = 0.42, linetype = "dashed", colour = "#7B858B") +
    geom_segment(data = delta, aes(x = 9, xend = delta_to_x(delta_F), y = row, yend = row, colour = direction),
                 linewidth = 1.35, lineend = "round") +
    geom_point(data = delta, aes(x = delta_to_x(delta_F), y = row, fill = direction), shape = 21,
               size = 4.0, stroke = 0.45, colour = "white") +
    geom_text(data = delta,
              aes(x = ifelse(delta_F >= 0, delta_to_x(delta_F) + 0.22, delta_to_x(delta_F) - 0.22),
                  y = row, label = delta_label, colour = direction),
              hjust = ifelse(delta$delta_F >= 0, 0, 1), family = "Arial", size = 2.75, fontface = "bold") +
    scale_fill_manual(values = colors, guide = "none") +
    scale_colour_manual(values = colors, guide = "none") +
    scale_x_continuous(limits = c(-4.5, 12.3), breaks = c(1:5, 7, 8, 9, 10, 11),
                       labels = c("1", "2", "3", "4", "5", "-0.8", "-0.4", "0", "+0.4", "+0.8"), expand = c(0, 0)) +
    scale_y_reverse(limits = c(3.45, -0.25), breaks = NULL) +
    theme_classic(base_size = 7.2, base_family = "Arial") +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
          axis.text.y = element_blank(), axis.text.x = element_text(size = 6.3, colour = "#333333", margin = margin(t = 1)),
          panel.grid = element_blank(), plot.margin = margin(6, 5, 4, 4))

  base <- file.path(panel_dir, "fig", tag)
  svglite(paste0(base, ".svg"), width = 128 / 25.4, height = 58 / 25.4); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = 128 / 25.4, height = 58 / 25.4, family = "Arial"); print(plot); dev.off()
  ragg::agg_png(paste0(base, ".png"), width = 128 / 25.4, height = 58 / 25.4, units = "in", res = 600); print(plot); dev.off()
}

draw_case(
  domain_name = "RNA",
  constraints = c("RNA01-02", "RNA06", "RNA09"),
  relation_labels = c(`RNA01-02` = "RNA01-02  Edit adjacency", RNA06 = "RNA06  Partner specificity", RNA09 = "RNA09  Remote structure"),
  model_a = "RNA-FM", model_b = "RiNALMo-Giga", tag = "Fig4_PanelC_RNA_rank_advantage_strip"
)

draw_case(
  domain_name = "protein structure",
  constraints = c("ST02", "ST04", "ST05"),
  relation_labels = c(ST02 = "ST02  Structural neighborhood", ST04 = "ST04  Long-range contacts", ST05 = "ST05  Rigid-body invariance"),
  model_a = "ESM-IF1", model_b = "SSLM", tag = "Fig4_PanelC_structure_rank_advantage_strip"
)
