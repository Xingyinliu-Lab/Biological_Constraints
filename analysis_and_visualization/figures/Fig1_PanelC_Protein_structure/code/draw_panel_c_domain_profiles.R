#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: draw_panel_c_domain_profiles.R <S4.tsv> <output_dir>")
}

input_tsv <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.delim(input_tsv, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("Constraint", "Short_name", "Model", "F")
if (!all(required %in% names(d))) {
  stop("S4 table is missing required columns")
}
if (anyNA(d[c("Constraint", "Short_name", "Model", "F")]) ||
    any(!is.finite(d$F)) || any(d$F < -1 | d$F > 1) ||
    anyDuplicated(d[c("Constraint", "Model")])) {
  stop("Invalid or duplicated S4 rows")
}

d$Score <- d$F
d$Domain <- ifelse(
  grepl("^DNA", d$Constraint), "DNA",
  ifelse(grepl("^RNA", d$Constraint), "RNA",
         ifelse(grepl("^PS", d$Constraint), "Protein sequence", "Protein structure"))
)

domain_spec <- list(
  DNA = list(
    constraints = c("DNA01", "DNA02-02", "DNA03", "DNA06"),
    file = "Fig1_PanelC_DNA_profile.pdf", height_mm = 58
  ),
  RNA = list(
    constraints = c("RNA01-02", "RNA04", "RNA05", "RNA06", "RNA09", "RNA10"),
    file = "Fig1_PanelC_RNA_profile.pdf", height_mm = 60
  ),
  `Protein sequence` = list(
    constraints = c("PS01", "PS02A", "PS02B", "PS04", "PS05", "PS08", "PS09", "PS10"),
    file = "Fig1_PanelC_Protein_sequence_profile.pdf", height_mm = 62
  ),
  `Protein structure` = list(
    constraints = c("ST01", "ST02", "ST04", "ST05"),
    file = "Fig1_PanelC_Protein_structure_profile.pdf", height_mm = 58
  )
)

# Five restrained qualitative hues: adjacent trajectories remain identifiable
# in both the plot and its direct legend.
model_colours <- c("#2E6F9E", "#C26050", "#4B9679", "#936EAA", "#BD8438")
model_shapes <- c(16, 17, 15, 18, 8)

make_profile <- function(domain_name, spec) {
  x <- d[d$Domain == domain_name, , drop = FALSE]
  if (!setequal(unique(x$Constraint), spec$constraints)) {
    stop("Unexpected constraint set for ", domain_name)
  }
  if (any(table(x$Constraint) != 5L) || length(unique(x$Model)) != 5L) {
    stop(domain_name, " must contain exactly five models per constraint")
  }

  model_order <- unique(x$Model)
  short_name_map <- setNames(x$Short_name, x$Constraint)
  axis_labels <- gsub(" ", "\n", short_name_map[spec$constraints], fixed = TRUE)
  axis_labels <- gsub("[-–]", "\n", axis_labels)
  names(axis_labels) <- spec$constraints

  x$Constraint <- factor(x$Constraint, levels = spec$constraints)
  x$Model <- factor(x$Model, levels = model_order)
  colour_map <- setNames(model_colours, model_order)
  shape_map <- setNames(model_shapes, model_order)

  p <- ggplot(x, aes(Constraint, Score, group = Model, colour = Model, shape = Model)) +
    geom_hline(yintercept = c(0.5, 1.0), colour = "#E4E7E9",
               linewidth = 0.28) +
    geom_hline(yintercept = 0, colour = "#5E6870", linewidth = 0.36,
               linetype = "32") +
    annotate("text", x = 1, y = 0, label = "chance level", hjust = 0,
             vjust = -0.55, size = 1.65, colour = "#5E6870", family = "Arial") +
    geom_line(linewidth = 0.52, alpha = 0.88, lineend = "round") +
    geom_point(size = 1.75, stroke = 0) +
    scale_colour_manual(values = colour_map, drop = FALSE) +
    scale_shape_manual(values = shape_map, drop = FALSE) +
    scale_x_discrete(labels = axis_labels, expand = expansion(add = 0.35)) +
    scale_y_continuous(
      limits = c(-0.1, 1.02),
      breaks = c(0, 0.5, 1.0),
      labels = c("0", "0.5", "1.0"),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = domain_name,
      x = NULL,
      y = "Formal score",
      colour = NULL
    ) +
    guides(colour = guide_legend(
      nrow = 2, byrow = TRUE,
      override.aes = list(linewidth = 0.68, size = 1.9)
    ), shape = "none") +
    theme_classic(base_family = "Arial", base_size = 6.1) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(size = 7.6, face = "bold", colour = "#20282E",
                                margin = margin(b = 2.5)),
      axis.title.y = element_text(size = 5.5, colour = "#424A50",
                                  margin = margin(r = 3.5)),
      axis.text.y = element_text(size = 5.2, colour = "#4C545A"),
      axis.text.x = element_text(size = 5.0, colour = "#30383E",
                                 lineheight = 0.86, margin = margin(t = 3.5)),
      axis.line = element_line(colour = "#737B80", linewidth = 0.3),
      axis.ticks.y = element_line(colour = "#737B80", linewidth = 0.28),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.justification = "left",
      legend.box.just = "left",
      legend.margin = margin(0, 0, 2.5, 0),
      legend.key.width = grid::unit(5.2, "mm"),
      legend.key.height = grid::unit(2.4, "mm"),
      legend.text = element_text(size = 4.65, colour = "#3F474C"),
      plot.margin = margin(3.5, 6.0, 3.5, 4)
    )

  output_pdf <- file.path(output_dir, spec$file)
  grDevices::cairo_pdf(
    output_pdf, width = 89 / 25.4, height = spec$height_mm / 25.4,
    family = "Arial", onefile = TRUE
  )
  print(p)
  grDevices::dev.off()

  output_svg <- sub("\\.pdf$", ".svg", output_pdf)
  svglite::svglite(output_svg, width = 89 / 25.4, height = spec$height_mm / 25.4)
  print(p)
  grDevices::dev.off()

  output_png <- sub("\\.pdf$", ".png", output_pdf)
  ragg::agg_png(output_png, width = 89 / 25.4, height = spec$height_mm / 25.4,
                 units = "in", res = 450, background = "white")
  print(p)
  grDevices::dev.off()
  message(output_pdf, "\t", output_svg, "\t", output_png)
}

invisible(Map(make_profile, names(domain_spec), domain_spec))
