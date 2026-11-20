#!/usr/bin/env Rscript
# Figure 5C | Comparator rank among the five formal models for each constraint.
# Ranks are computed only within a constraint (five models plus its one formal comparator).

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(svglite)
  library(ragg)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else ""
panel_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(panel_dir)) stop("Run this script with Rscript from its packaged location")
fig4c_dir <- normalizePath(file.path(panel_dir, "..", "..", "analysis", "fig4c_statistical_tests"), mustWork = TRUE)
source_s4 <- file.path(fig4c_dir, "data", "Supplementary_Table_S4_Complete_Model_Performance.tsv")
source_s3 <- file.path(panel_dir, "data", "Supplementary_Table_S3_Comparator_Overview.tsv")
data_dir <- file.path(panel_dir, "data")
fig_dir <- file.path(panel_dir, "fig")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# One currently formal comparator per Constraint. F values are copied from the
# corresponding released comparator result; Level refers to comparator knowledge access.
comparators <- tribble(
  ~Constraint, ~Comparator, ~Comparator_F, ~Level,
  "DNA01", "4-mer RC permutation", 1.000000, "Level 0",
  "DNA02-02", "Hamming", 1.000000, "Level 0",
  "DNA03", "Exact core containment", 1.000000, "Level 0",
  "DNA06", "Genetic code", 1.000000, "Level 0",
  "RNA01-02", "Exact Hamming", 1.000000, "Level 0",
  "RNA04", "CDS genetic code", 1.000000, "Level 0",
  "RNA05", "ViennaRNA pairing", 0.702226, "Level 1",
  "RNA06", "ViennaRNA partner probability", 0.947340, "Level 1",
  "RNA09", "RNApaln", 0.787880, "Level 1",
  "RNA10", "Positional k-mer ridge", 0.331200, "Level 1",
  "PS01", "NW+BLOSUM62", 0.812950, "Level 1",
  "PS02A", "NW+BLOSUM62", 0.349632, "Level 1",
  "PS02B", "NW+BLOSUM62", 0.379730, "Level 1",
  "PS04", "NW+BLOSUM62", 0.237720, "Level 1",
  "PS05", "NW+BLOSUM62", 0.288054, "Level 1",
  "PS08", "Hamming", 1.000000, "Level 0",
  "PS09", "IUPred2A", 0.591430, "Level 1",
  "PS10", "TANGO", 0.630843, "Level 1",
  "ST01", "Geometry white-box", 1.000000, "Level 0",
  "ST02", "TM-align", 0.617970, "Level 1",
  "ST04", "Geometry white-box", 1.000000, "Level 0",
  "ST05", "Geometry white-box", 1.000000, "Level 0"
) %>% mutate(
  Domain = case_when(
    startsWith(Constraint, "DNA") ~ "DNA",
    startsWith(Constraint, "RNA") ~ "RNA",
    startsWith(Constraint, "PS") ~ "Protein sequence",
    TRUE ~ "Protein structure"
  )
)

models <- read_tsv(source_s4, show_col_types = FALSE) %>%
  select(Constraint, Model, F) %>%
  mutate(source = "Model")

comparator_rows <- comparators %>%
  transmute(Constraint, Model = "Comparator", F = Comparator_F, source = "Comparator")

rank_data <- bind_rows(models, comparator_rows) %>%
  group_by(Constraint) %>%
  mutate(within_constraint_rank = rank(-F, ties.method = "average"), panel_size = n()) %>%
  ungroup() %>%
  left_join(comparators %>% select(Constraint, Comparator, Comparator_F, Level, Domain), by = "Constraint")

comparator_plot_data <- rank_data %>%
  filter(source == "Comparator") %>%
  mutate(
    Constraint = factor(Constraint, levels = comparators$Constraint)
  )

# Purposeful blank gutters, rather than boxes, separate the four biological domains.
domain_offsets <- c("DNA" = 0, "RNA" = 1, "Protein sequence" = 2, "Protein structure" = 3)
comparator_plot_data <- comparator_plot_data %>%
  group_by(Domain) %>%
  mutate(x = row_number() + sum(vapply(names(domain_offsets)[seq_len(match(first(Domain), names(domain_offsets))) - 1],
                                        function(z) sum(comparators$Domain == z), numeric(1))) + domain_offsets[first(Domain)]) %>%
  ungroup()

short_names <- read_tsv(source_s3, show_col_types = FALSE) %>%
  distinct(Constraint, Short_name)
comparator_plot_data <- comparator_plot_data %>% left_join(short_names, by = "Constraint")

write_tsv(comparator_plot_data, file.path(data_dir, "comparator_rank_source_data.tsv"))

domain_cols <- c(
  "DNA" = "#4E789D",
  "RNA" = "#4C998A",
  "Protein sequence" = "#C58C4C",
  "Protein structure" = "#806AA4"
)

domain_spans <- comparator_plot_data %>%
  group_by(Domain) %>%
  summarise(xmin = min(x) - 0.35, xmax = max(x) + 0.35, .groups = "drop") %>%
  mutate(xmid = (xmin + xmax) / 2)

p <- ggplot(comparator_plot_data, aes(x = x, y = within_constraint_rank)) +
  geom_hline(yintercept = 1:6, linewidth = 0.28, colour = "#D9DEE3") +
  geom_segment(data = domain_spans,
               aes(x = xmin, xend = xmax, y = 0.47, yend = 0.47, colour = Domain),
               inherit.aes = FALSE, linewidth = 1.15, lineend = "butt") +
  geom_text(data = domain_spans, aes(x = xmid, y = 0.15, label = Domain, colour = Domain),
            inherit.aes = FALSE, family = "Arial", fontface = "bold", size = 3.15) +
  geom_segment(aes(xend = x, y = 6, yend = within_constraint_rank),
               colour = "#C6CED6", linewidth = 0.42, lineend = "butt") +
  geom_point(aes(fill = Level), shape = 21, size = 3.15, colour = "#FFFFFF", stroke = 0.1) +
  scale_fill_manual(
    values = c("Level 0" = "#96A8B7", "Level 1" = "#314B68"),
    labels = c("Level 0: intrinsic", "Level 1: explicit abstract"),
    name = NULL
  ) +
  scale_colour_manual(values = domain_cols, guide = "none") +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = NA, size = 3.2))) +
  scale_x_continuous(
    breaks = comparator_plot_data$x,
    labels = comparator_plot_data$Short_name,
    expand = expansion(mult = c(0.018, 0.018))
  ) +
  scale_y_reverse(
    breaks = 1:6,
    limits = c(6.15, -0.25),
    labels = c("1  highest", 2:5, "6  lowest"),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Comparator rank",
    title = "Comparator rank across prespecified biological constraints"
  ) +
  theme_classic(base_family = "Arial", base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10.5, colour = "#17263A", hjust = 0),
    axis.title.y = element_text(size = 8.7, colour = "#263547", margin = margin(r = 9)),
    axis.text.x = element_text(size = 5.8, colour = "#2C3B4A", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 7.4, colour = "#425466"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 7.8),
    legend.key.width = unit(11, "pt"),
    legend.key.height = unit(9, "pt"),
    plot.margin = margin(7, 10, 5, 9)
  )

base <- file.path(fig_dir, "Fig5_PanelC_comparator_rank_landscape")
ggsave(paste0(base, ".pdf"), p, width = 210, height = 70, units = "mm", device = cairo_pdf)
svglite(paste0(base, ".svg"), width = 210 / 25.4, height = 70 / 25.4)
print(p)
dev.off()
ragg::agg_png(paste0(base, ".png"), width = 210 / 25.4, height = 70 / 25.4, units = "in", res = 600)
print(p)
dev.off()
