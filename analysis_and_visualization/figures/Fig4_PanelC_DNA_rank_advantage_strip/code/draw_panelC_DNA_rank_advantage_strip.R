suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
  library(ragg)
})

# Figure contract
# Claim: the rank and F advantage of the same model pair reverse across DNA
# relations. Rank is shown within the full five-model panel; Delta F preserves
# the magnitude and direction of the pairwise reversal.
# Archetype: quantitative case strip; editable SVG/PDF and PNG preview.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script with Rscript")
panel_dir <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), ".."), mustWork = TRUE)
fig4c_dir <- normalizePath(file.path(panel_dir, "..", "..", "analysis", "fig4c_statistical_tests"), mustWork = TRUE)
dir.create(file.path(panel_dir, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(panel_dir, "fig"), recursive = TRUE, showWarnings = FALSE)

scores <- read.delim(file.path(fig4c_dir, "data", "Supplementary_Table_S4_Complete_Model_Performance.tsv"),
                     check.names = FALSE, stringsAsFactors = FALSE)
scores$rank <- ave(-scores$F, scores$Constraint, FUN = function(x) rank(x, ties.method = "min"))

constraints <- c("DNA01", "DNA02-02", "DNA06")
labels <- c(
  DNA01 = "DNA01  Reverse complement",
  `DNA02-02` = "DNA02-02  Edit adjacency",
  DNA06 = "DNA06  Codon consequence"
)
models <- c("Caduceus-PS", "DNABERT-2")

case <- scores[scores$Constraint %in% constraints & scores$Model %in% models,
               c("Constraint", "Short_name", "Model", "F", "F_95_CI", "rank"), drop = FALSE]
if (nrow(case) != 6L) stop("Expected the selected model pair for three constraints")
case$row <- match(case$Constraint, constraints)
case$rank_label <- ifelse(case$Model == "Caduceus-PS", "C", "D")

delta <- do.call(rbind, lapply(constraints, function(id) {
  d <- case[case$Constraint == id, ]
  c_score <- d$F[d$Model == "Caduceus-PS"]
  d_score <- d$F[d$Model == "DNABERT-2"]
  data.frame(
    Constraint = id,
    relation = labels[[id]],
    row = match(id, constraints),
    delta_F = c_score - d_score,
    delta_label = sprintf("%+.2f", c_score - d_score),
    direction = ifelse(c_score > d_score, "Caduceus-PS higher", "DNABERT-2 higher"),
    stringsAsFactors = FALSE
  )
}))

# Rank strip includes all positions; highlighted points identify the model pair.
rank_slots <- expand.grid(row = 1:3, rank = 1:5)
rank_slots$occupied_by_pair <- mapply(function(row, rank) any(case$row == row & case$rank == rank),
                                      rank_slots$row, rank_slots$rank)

write.table(case, file.path(panel_dir, "data", "panelC_DNA_rank_advantage_source_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(delta, file.path(panel_dir, "data", "panelC_DNA_rank_advantage_deltas.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Common coordinates: ranks occupy x = 1..5. Delta F is mapped to x = 7..11,
# with the zero reference at x = 9 and a fixed [-0.8, +0.8] display range.
delta_to_x <- function(x) 9 + x / 0.8 * 2
colors <- c("Caduceus-PS" = "#356C7D", "DNABERT-2" = "#B56A3A")
direction_colors <- c("Caduceus-PS higher" = "#356C7D", "DNABERT-2 higher" = "#B56A3A")

theme_strip <- function() {
  theme_classic(base_size = 7.2, base_family = "Arial") +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.title = element_blank(), axis.text.y = element_blank(),
      axis.text.x = element_text(size = 6.3, colour = "#333333", margin = margin(t = 1)),
      panel.grid = element_blank(),
      plot.margin = margin(6, 5, 4, 4)
    )
}

p <- ggplot() +
  # Left heading and directional rank scale.
  annotate("text", x = 3, y = -0.05, label = "Rank among five models", family = "Arial",
           size = 2.9, fontface = "bold", colour = "#223A4B") +
  annotate("text", x = 1, y = 0.25, label = "highest", family = "Arial", size = 2.25,
           hjust = 0, colour = "#66727A") +
  annotate("text", x = 5, y = 0.25, label = "lowest", family = "Arial", size = 2.25,
           hjust = 1, colour = "#66727A") +
  # Right heading and model-pair difference scale.
  annotate("text", x = 9, y = -0.05, label = "Pairwise advantage (Delta F)", family = "Arial",
           size = 2.9, fontface = "bold", colour = "#223A4B") +
  annotate("text", x = 7.25, y = 0.25, label = "DNABERT-2", family = "Arial", size = 1.95,
           hjust = 0, colour = colors[["DNABERT-2"]]) +
  annotate("text", x = 10.75, y = 0.25, label = "Caduceus-PS", family = "Arial", size = 1.95,
           hjust = 1, colour = colors[["Caduceus-PS"]]) +
  # Relation labels.
  geom_text(data = delta, aes(x = -4.3, y = row, label = relation), hjust = 0,
            size = 2.75, family = "Arial", colour = "#252525") +
  # Five available rank slots, with the three other models visually de-emphasized.
  geom_point(data = rank_slots[!rank_slots$occupied_by_pair, ], aes(x = rank, y = row),
             shape = 21, size = 3.15, stroke = 0.45, fill = "#E2E6E8", colour = "#C7CED1") +
  geom_point(data = case, aes(x = rank, y = row, fill = Model), shape = 21,
             size = 4.35, stroke = 0.45, colour = "white") +
  geom_text(data = case, aes(x = rank, y = row, label = rank_label), family = "Arial",
            size = 2.45, fontface = "bold", colour = "white", vjust = 0.35) +
  # A single baseline and signed advantage vector per relation.
  geom_vline(xintercept = 9, linewidth = 0.42, linetype = "dashed", colour = "#7B858B") +
  geom_segment(data = delta,
               aes(x = 9, xend = delta_to_x(delta_F), y = row, yend = row, colour = direction),
               linewidth = 1.35, lineend = "round") +
  geom_point(data = delta, aes(x = delta_to_x(delta_F), y = row, fill = direction),
             shape = 21, size = 4.0, stroke = 0.45, colour = "white") +
  geom_text(data = delta,
            aes(x = ifelse(delta_F >= 0, delta_to_x(delta_F) + 0.22, delta_to_x(delta_F) - 0.22),
                y = row, label = delta_label, colour = direction),
            hjust = ifelse(delta$delta_F >= 0, 0, 1), family = "Arial", size = 2.75, fontface = "bold") +
  scale_fill_manual(values = c(colors, direction_colors), guide = "none") +
  scale_colour_manual(values = c(colors, direction_colors), guide = "none") +
  scale_x_continuous(
    limits = c(-4.5, 12.3),
    breaks = c(1:5, 7, 8, 9, 10, 11),
    labels = c("1", "2", "3", "4", "5", "-0.8", "-0.4", "0", "+0.4", "+0.8"),
    expand = c(0, 0)
  ) +
  scale_y_reverse(limits = c(3.45, -0.25), breaks = NULL) +
  theme_strip()

base <- file.path(panel_dir, "fig", "Fig4_PanelC_DNA_rank_advantage_strip")
svglite(paste0(base, ".svg"), width = 128 / 25.4, height = 58 / 25.4)
print(p)
dev.off()
grDevices::cairo_pdf(paste0(base, ".pdf"), width = 128 / 25.4, height = 58 / 25.4, family = "Arial")
print(p)
dev.off()
ragg::agg_png(paste0(base, ".png"), width = 128 / 25.4, height = 58 / 25.4, units = "in", res = 600)
print(p)
dev.off()
