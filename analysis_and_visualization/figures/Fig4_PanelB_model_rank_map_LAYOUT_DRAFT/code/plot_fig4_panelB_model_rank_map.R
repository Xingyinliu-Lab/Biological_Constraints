suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(svglite)
  library(ragg)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script with Rscript")
panel_dir <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), ".."), mustWork = TRUE)
fig4c_dir <- normalizePath(file.path(panel_dir, "..", "..", "analysis", "fig4c_statistical_tests"), mustWork = TRUE)
input_file <- file.path(fig4c_dir, "data", "Supplementary_Table_S4_Complete_Model_Performance.tsv")
out_dir <- panel_dir
dir.create(file.path(out_dir, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "fig"), recursive = TRUE, showWarnings = FALSE)

domain_map <- c(
  DNA01 = "DNA", `DNA02-02` = "DNA", DNA03 = "DNA", DNA06 = "DNA",
  `RNA01-02` = "RNA", RNA04 = "RNA", RNA05 = "RNA", RNA06 = "RNA", RNA09 = "RNA", RNA10 = "RNA",
  PS01 = "Protein sequence", PS02A = "Protein sequence", PS02B = "Protein sequence", PS04 = "Protein sequence",
  PS05 = "Protein sequence", PS08 = "Protein sequence", PS09 = "Protein sequence", PS10 = "Protein sequence",
  ST01 = "Protein structure", ST02 = "Protein structure", ST04 = "Protein structure", ST05 = "Protein structure"
)
domain_levels <- c("DNA", "RNA", "Protein sequence", "Protein structure")
domain_cols <- c(
  "DNA" = "#557FA3", "RNA" = "#4B947D",
  "Protein sequence" = "#B0804C", "Protein structure" = "#81699C"
)
display_name <- c(
  "Caduceus-PS" = "Cad-PS", "CrossDNA-8M" = "CrossDNA", "DNABERT-2" = "DNABERT-2", "JanusDNA" = "JanusDNA", "NT-v2" = "NT-v2",
  "ERNIE-RNA" = "ERNIE-RNA", "GB-RNA" = "GB-RNA", "Orthrus-4-track" = "Orthrus", "RiNALMo-Giga" = "RiNALMo", "RNA-FM" = "RNA-FM",
  "CARP-640M" = "CARP", "ESM-2-650M" = "ESM-2", "ESM3-seq" = "ESM3-seq", "ESM-C-600M" = "ESM-C", "ProtT5" = "ProtT5",
  "ESM-IF1" = "ESM-IF1", "ESM3-struct" = "ESM3-struct", "ProteinMPNN-CA" = "MPNN-CA", "ProToken" = "ProToken", "SSLM" = "SSLM"
)

scores <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
scores$domain <- unname(domain_map[scores$Constraint])
stopifnot(!anyNA(scores$domain), all(vapply(split(scores, scores$Constraint), nrow, integer(1)) == 5L))
scores$domain <- factor(scores$domain, levels = domain_levels)

scores$rank <- ave(scores$F, scores$Constraint, FUN = function(x) rank(-x, ties.method = "average"))

rank_map <- do.call(rbind, lapply(split(scores, interaction(scores$domain, scores$Model, drop = TRUE)), function(d) {
  n_constraints <- length(unique(d$Constraint))
  data.frame(
    domain = as.character(d$domain[1]),
    Model = d$Model[1],
    display_model = unname(display_name[d$Model[1]]),
    n_constraints = n_constraints,
    rank_1 = sum(d$rank == 1),
    rank_2 = sum(d$rank == 2),
    rank_3 = sum(d$rank == 3),
    rank_4 = sum(d$rank == 4),
    rank_5 = sum(d$rank == 5),
    mean_rank = mean(d$rank),
    rank_dispersion = sd(d$rank),
    stringsAsFactors = FALSE
  )
}))
rank_map$domain <- factor(rank_map$domain, levels = domain_levels)
rank_map <- rank_map[order(rank_map$domain, rank_map$mean_rank, rank_map$rank_dispersion), ]

write.table(rank_map,
  file.path(out_dir, "data", "Fig4_PanelB_model_rank_map.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

theme_set(theme_bw(base_size = 8, base_family = "Arial") + theme(
  panel.grid = element_blank(),
  panel.border = element_blank(),
  axis.line = element_line(linewidth = 0.35, colour = "#28323A"),
  axis.ticks = element_line(linewidth = 0.35, colour = "#28323A"),
  axis.text = element_text(colour = "#28323A", size = 7),
  axis.title = element_text(colour = "#1B252C", size = 8),
  plot.title = element_text(face = "bold", size = 10, colour = "#15212B", hjust = 0),
  plot.subtitle = element_text(size = 7.5, colour = "#52616B", hjust = 0),
  plot.caption = element_text(size = 6.2, colour = "#52616B", hjust = 0),
  legend.position = "top",
  legend.title = element_blank(),
  legend.text = element_text(size = 7),
  plot.margin = margin(5.5, 12, 5.5, 5.5)
))

p <- ggplot(rank_map, aes(x = mean_rank, y = rank_dispersion, colour = domain, label = display_model)) +
  geom_vline(xintercept = 3, colour = "#B7C0C6", linetype = "dashed", linewidth = 0.35) +
  geom_hline(yintercept = 1, colour = "#B7C0C6", linetype = "dashed", linewidth = 0.35) +
  annotate("text", x = 1.12, y = 2.43, label = "\u2190  better overall performance", hjust = 0,
    family = "Arial", size = 3.0, colour = "#52616B") +
  annotate("text", x = 1.12, y = 0.14, label = "leading\nstable", hjust = 0,
    family = "Arial", size = 2.55, colour = "#87939A") +
  annotate("text", x = 1.12, y = 2.05, label = "leading\nrelation-dependent", hjust = 0,
    family = "Arial", size = 2.55, colour = "#87939A") +
  annotate("text", x = 4.88, y = 0.14, label = "trailing\nstable", hjust = 1,
    family = "Arial", size = 2.55, colour = "#87939A") +
  annotate("text", x = 4.88, y = 2.05, label = "trailing\nrelation-dependent", hjust = 1,
    family = "Arial", size = 2.55, colour = "#87939A") +
  geom_point(size = 3.5, stroke = 0, alpha = 0.97) +
  geom_text_repel(
    family = "Arial", size = 2.85, colour = "#1E2932",
    min.segment.length = 0, segment.colour = "#B8C0C5",
    box.padding = 0.25, point.padding = 0.25, max.overlaps = Inf,
    seed = 20261022, show.legend = FALSE
  ) +
  scale_colour_manual(values = domain_cols, breaks = domain_levels) +
  scale_x_continuous(limits = c(0.9, 5.1), breaks = 1:5, expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "b  Model rank map",
    subtitle = "Left: higher average rank. Up: greater rank variation across biological constraints.",
    x = "Mean rank across constraints (1 = highest)",
    y = "Rank dispersion (s.d. across constraints)",
    caption = "Each point is one model within its biological domain. Rankings use formal F scores within each constraint."
  ) +
  coord_cartesian(clip = "off")

base <- file.path(out_dir, "fig", "Fig4_PanelB_model_rank_map_LAYOUT_DRAFT")
ggsave(paste0(base, ".png"), p, width = 150, height = 115, units = "mm", dpi = 450, bg = "white")
svglite(paste0(base, ".svg"), width = 150 / 25.4, height = 115 / 25.4, bg = "white")
print(p)
dev.off()
grDevices::cairo_pdf(paste0(base, ".pdf"), width = 150 / 25.4, height = 115 / 25.4, family = "Arial")
print(p)
dev.off()
