#!/usr/bin/env Rscript
# Figure 5C draft | Formal-score distributions for the unique Comparator and five models.

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(svglite)
  library(ragg)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else ""
panel_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(panel_dir)) stop("Run this script with Rscript from its packaged location")
fig4c_dir <- normalizePath(file.path(panel_dir, "..", "..", "analysis", "fig4c_statistical_tests"), mustWork = TRUE)
data_dir <- file.path(panel_dir, "data")
fig_dir <- file.path(panel_dir, "fig")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

source_s4 <- file.path(fig4c_dir, "data", "Supplementary_Table_S4_Complete_Model_Performance.tsv")

comparators <- tribble(
  ~Constraint, ~Comparator_F,
  "DNA01", 1.000000, "DNA02-02", 1.000000, "DNA03", 1.000000, "DNA06", 1.000000,
  "RNA01-02", 1.000000, "RNA04", 1.000000, "RNA05", 0.702226, "RNA06", 0.947340,
  "RNA09", 0.787880, "RNA10", 0.331200,
  "PS01", 0.812950, "PS02A", 0.349632, "PS02B", 0.379730, "PS04", 0.237720,
  "PS05", 0.288054, "PS08", 1.000000, "PS09", 0.591430, "PS10", 0.630843,
  "ST01", 1.000000, "ST02", 0.617970, "ST04", 1.000000, "ST05", 1.000000
) %>% mutate(
  Domain = case_when(
    startsWith(Constraint, "DNA") ~ "DNA",
    startsWith(Constraint, "RNA") ~ "RNA",
    startsWith(Constraint, "PS") ~ "Protein sequence",
    TRUE ~ "Protein structure"
  )
)

domain_cols <- c(
  "DNA" = "#5B82A5", "RNA" = "#569C91",
  "Protein sequence" = "#C78E4B", "Protein structure" = "#8A73A8"
)

score_data <- read_tsv(source_s4, show_col_types = FALSE) %>%
  select(Constraint, Model, F) %>%
  mutate(Method = Model, method_type = "Model") %>%
  select(-Model) %>%
  bind_rows(comparators %>% transmute(Constraint, F = Comparator_F, Method = "Comparator", method_type = "Comparator")) %>%
  left_join(comparators %>% select(Constraint, Domain), by = "Constraint")

write_tsv(score_data, file.path(data_dir, "comparator_model_score_source_data.tsv"))

for (domain_name in names(domain_cols)) {
  d <- score_data %>% filter(Domain == domain_name)
  model_order <- d %>%
    filter(method_type == "Model") %>%
    group_by(Method) %>%
    summarise(median_F = median(F), .groups = "drop") %>%
    arrange(desc(median_F), Method) %>%
    pull(Method)
  ordered_methods <- c("Comparator", model_order)
  d <- d %>% mutate(Method = factor(Method, levels = ordered_methods))
  lim <- range(c(d$F, 0), na.rm = TRUE)
  pad <- max(0.06, diff(lim) * 0.07)
  lim <- c(max(-1, lim[1] - pad), min(1, lim[2] + pad))

  p <- ggplot(d, aes(x = Method, y = F)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.42, colour = "#7E8994") +
    geom_vline(xintercept = 1.5, linewidth = 0.45, colour = "#84909C") +
    geom_boxplot(
      aes(fill = method_type), width = 0.42, outlier.shape = NA,
      alpha = 0.13, linewidth = 0.36, colour = "#4E5D6D"
    ) +
    geom_point(
      data = d %>% filter(method_type == "Model"),
      position = position_jitter(width = 0.10, height = 0, seed = 20261023),
      shape = 21, size = 2.55, stroke = 0.16, colour = "white",
      fill = domain_cols[[domain_name]], alpha = 0.80
    ) +
    geom_point(
      data = d %>% filter(method_type == "Comparator"),
      position = position_jitter(width = 0.10, height = 0, seed = 20261023),
      shape = 21, size = 2.55, stroke = 0.16, colour = "white",
      fill = "#66798A", alpha = 0.92
    ) +
    scale_x_discrete(limits = ordered_methods, drop = FALSE) +
    scale_fill_manual(values = c("Comparator" = "#66798A", "Model" = domain_cols[[domain_name]]), guide = "none") +
    annotate("text", x = 6.35, y = 0, label = "F = 0", hjust = 1.05, vjust = -0.55,
             family = "Arial", size = 2.65, colour = "#5B6875") +
    coord_cartesian(ylim = lim, clip = "off") +
    labs(title = domain_name,
         subtitle = "Comparator at left; models ordered by median formal score",
         x = NULL, y = "Formal score, F") +
    theme_classic(base_family = "Arial", base_size = 8.5) +
    theme(
      plot.title = element_text(face = "bold", size = 10.5, colour = domain_cols[[domain_name]], hjust = 0),
      plot.subtitle = element_text(size = 7.2, colour = "#52606D", margin = margin(b = 5)),
      axis.title.y = element_text(size = 8.2, colour = "#334155", margin = margin(r = 7)),
      axis.text.y = element_text(size = 7.3, colour = "#52606D"),
      axis.text.x = element_text(size = 7.6, colour = "#1F2933", margin = margin(t = 5)),
      axis.ticks = element_blank(), axis.line = element_blank(),
      plot.margin = margin(7, 8, 6, 8)
    )

  stem <- gsub(" ", "_", domain_name)
  base <- file.path(fig_dir, paste0("Fig5_PanelC_comparator_score_", stem))
  ggsave(paste0(base, ".pdf"), p, width = 150, height = 70, units = "mm", device = cairo_pdf)
  svglite(paste0(base, ".svg"), width = 150 / 25.4, height = 70 / 25.4)
  print(p)
  dev.off()
  ragg::agg_png(paste0(base, ".png"), width = 150 / 25.4, height = 70 / 25.4, units = "in", res = 600)
  print(p)
  dev.off()
}
