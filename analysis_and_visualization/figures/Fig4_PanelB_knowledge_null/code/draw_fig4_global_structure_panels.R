suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

# Figure contract
# Core conclusion: overall model standing, knowledge heterogeneity and
# model-by-knowledge interaction coexist across the prespecified panel.
# Domains are retained as blocking strata, not used as the narrative axis.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script with Rscript")
panel_dir <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), ".."), mustWork = TRUE)
fig4c_dir <- normalizePath(file.path(panel_dir, "..", "..", "analysis", "fig4c_statistical_tests"), mustWork = TRUE)
out <- panel_dir
dir.create(file.path(out, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "fig"), recursive = TRUE, showWarnings = FALSE)

score_file <- file.path(out, "data", "formal_F_scores_with_rank_and_CI.tsv")
reversal_file <- file.path(out, "data", "global_pairwise_CI_supported_reversals.tsv")
knowledge_null_file <- file.path(fig4c_dir, "results", "global_knowledge_unevenness_null_draws.tsv")
knowledge_result_file <- file.path(fig4c_dir, "results", "global_knowledge_unevenness_permutation.tsv")
interaction_null_file <- file.path(fig4c_dir, "results", "global_model_knowledge_interaction_null_draws.tsv")
interaction_result_file <- file.path(fig4c_dir, "results", "global_model_knowledge_interaction_wild_bootstrap.tsv")

scores <- read.delim(score_file, check.names = FALSE, stringsAsFactors = FALSE)
reversals <- read.delim(reversal_file, check.names = FALSE, stringsAsFactors = FALSE)
knowledge_null <- read.delim(knowledge_null_file, check.names = FALSE, stringsAsFactors = FALSE)[[1]]
knowledge_result <- read.delim(knowledge_result_file, check.names = FALSE, stringsAsFactors = FALSE)
interaction_null <- read.delim(interaction_null_file, check.names = FALSE, stringsAsFactors = FALSE)[[1]]
interaction_result <- read.delim(interaction_result_file, check.names = FALSE, stringsAsFactors = FALSE)

stopifnot(nrow(scores) == 110L, nrow(reversals) == 40L)
stopifnot(sum(reversals$ci_supported_bidirectional_reversal) == 23L)

pal <- c(
  ink = "#1E2B35", muted = "#667681", grid = "#D9E0E4", null = "#B8C6CE",
  null_dark = "#7E929F", observed = "#B54A45", blue = "#507F9F",
  teal = "#4E927F", ochre = "#B78343", violet = "#806B99",
  pale_blue = "#EAF1F5", pale_teal = "#EAF3F0", pale_ochre = "#F5EFE6"
)
domain_cols <- c(
  "DNA" = pal[["blue"]], "RNA" = pal[["teal"]],
  "Protein sequence" = pal[["ochre"]], "Protein structure" = pal[["violet"]]
)
domain_order <- c("DNA", "RNA", "Protein sequence", "Protein structure")
constraint_order <- c(
  "DNA01", "DNA02-02", "DNA03", "DNA06",
  "RNA01-02", "RNA04", "RNA05", "RNA06", "RNA09", "RNA10",
  "PS01", "PS02A", "PS02B", "PS04", "PS05", "PS08", "PS09", "PS10",
  "ST01", "ST02", "ST04", "ST05"
)
constraint_label <- c(
  DNA01 = "Reverse\ncomplement", `DNA02-02` = "Edit\nadjacency", DNA03 = "Context\nextension", DNA06 = "Codon\nconsequence",
  `RNA01-02` = "Edit\nadjacency", RNA04 = "Coding\nconsequence", RNA05 = "Pairing\nstate", RNA06 = "Partner\nspecificity", RNA09 = "Remote\nstructure", RNA10 = "Exon\nrecognition",
  PS01 = "Family\nidentity", PS02A = "Clan\nhomology", PS02B = "Remote\nhomology", PS04 = "Topology", PS05 = "Structural\nresemblance", PS08 = "Edit\nadjacency", PS09 = "Order–\ndisorder", PS10 = "Amyloid\npropensity",
  ST01 = "Global fold\nsimilarity", ST02 = "Structural\nneighbourhood", ST04 = "Long-range\ncontacts", ST05 = "Rigid-body\ninvariance"
)
model_order <- list(
  DNA = c("Caduceus-PS", "CrossDNA-8M", "DNABERT-2", "JanusDNA", "NT-v2"),
  RNA = c("GB-RNA", "RiNALMo-Giga", "ERNIE-RNA", "RNA-FM", "Orthrus-4-track"),
  `Protein sequence` = c("ESM3-seq", "ProtT5", "CARP-640M", "ESM-2-650M", "ESM-C-600M"),
  `Protein structure` = c("SSLM", "ProteinMPNN-CA", "ESM-IF1", "ESM3-struct", "ProToken")
)
model_label <- c(
  "Caduceus-PS" = "Caduceus-PS", "CrossDNA-8M" = "CrossDNA", "DNABERT-2" = "DNABERT-2", "JanusDNA" = "JanusDNA", "NT-v2" = "NT-v2",
  "GB-RNA" = "GB-RNA", "RiNALMo-Giga" = "RiNALMo", "ERNIE-RNA" = "ERNIE-RNA", "RNA-FM" = "RNA-FM", "Orthrus-4-track" = "Orthrus",
  "ESM3-seq" = "ESM3-seq", "ProtT5" = "ProtT5", "CARP-640M" = "CARP", "ESM-2-650M" = "ESM-2", "ESM-C-600M" = "ESM-C",
  "SSLM" = "SSLM", "ProteinMPNN-CA" = "MPNN-CA", "ESM-IF1" = "ESM-IF1", "ESM3-struct" = "ESM3-struct", "ProToken" = "ProToken"
)

theme_pub <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = pal[["ink"]]),
      axis.ticks = element_line(linewidth = 0.3, colour = pal[["ink"]]),
      axis.text = element_text(colour = pal[["ink"]]),
      axis.title = element_text(colour = pal[["ink"]]),
      plot.title = element_text(face = "bold", colour = pal[["ink"]], size = base_size + 0.8, hjust = 0),
      plot.subtitle = element_text(colour = pal[["muted"]], size = base_size - 0.2, hjust = 0),
      plot.margin = margin(5, 6, 5, 5),
      panel.grid = element_blank()
    )
}

save_panel <- function(plot, name, width_mm, height_mm, dpi = 600) {
  base <- file.path(out, "fig", name)
  svglite::svglite(paste0(base, ".svg"), width = width_mm / 25.4, height = height_mm / 25.4, bg = "white")
  print(plot); dev.off()
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4, family = "Arial")
  print(plot); dev.off()
  ragg::agg_png(paste0(base, ".png"), width = width_mm / 25.4, height = height_mm / 25.4, units = "in", res = dpi, background = "white")
  print(plot); dev.off()
}

# Panel A: conceptual decomposition. This is explanatory, not a fitted model.
component_data <- data.frame(
  x = c(1, 2, 3),
  title = c("Knowledge\naccessibility", "Overall model\nstanding", "Knowledge-specific\nadvantage"),
  symbol = c("K(k)", "M(m)", "I(m,k)"),
  question = c("Do Constraints\ndiffer?", "Is model order\nshared?", "Does advantage depend\non knowledge?"),
  stringsAsFactors = FALSE
)
write.table(component_data, file.path(out, "data", "Fig4_PanelA_effect_components.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pA <- ggplot(component_data, aes(x = x, y = 1)) +
  annotate("text", x = 0.14, y = 1, label = "F(m,k)  =", hjust = 0, family = "Arial", fontface = "bold", size = 4.2, colour = pal[["ink"]]) +
  geom_tile(aes(fill = factor(x)), width = 0.78, height = 0.56, colour = "white", linewidth = 0.8) +
  geom_text(aes(label = symbol), y = 1.09, family = "Arial", fontface = "bold", size = 4.1, colour = pal[["ink"]]) +
  geom_text(aes(label = title), y = 0.98, family = "Arial", fontface = "bold", size = 2.9, lineheight = 0.95, colour = pal[["ink"]]) +
  geom_text(aes(label = question), y = 0.82, family = "Arial", size = 2.35, lineheight = 0.95, colour = pal[["muted"]]) +
  annotate("text", x = c(1.5, 2.5), y = 1, label = "+", family = "Arial", fontface = "bold", size = 5, colour = pal[["muted"]]) +
  scale_fill_manual(values = c(pal[["pale_blue"]], pal[["pale_teal"]], pal[["pale_ochre"]]), guide = "none") +
  coord_cartesian(xlim = c(0.1, 3.45), ylim = c(0.68, 1.33), clip = "off") +
  labs(title = "Three effects shape each formal score") +
  theme_void(base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 9, colour = pal[["ink"]], margin = margin(b = 4)), plot.margin = margin(5, 5, 5, 5))
save_panel(pA, "Fig4_PanelA_effect_components", 178, 47)

# Panel B-D: observed statistics against their calibrated null distributions.
set.seed(20261022)
rank_mats <- lapply(domain_order, function(dom) {
  d <- scores[scores$domain == dom, ]
  # Preserve the original analysis order exactly for deterministic replay.
  cons <- unique(d$Constraint)
  mods <- sort(unique(d$Model))
  mat <- do.call(rbind, lapply(cons, function(id) {
    z <- d[d$Constraint == id, c("Model", "rank")]
    z$rank[match(mods, z$Model)]
  }))
  colnames(mat) <- mods
  rownames(mat) <- cons
  mat
})
rank_variability <- function(mats) {
  sum(vapply(mats, function(mat) sum((mat - rep(colMeans(mat), each = nrow(mat)))^2), numeric(1)))
}
observed_order <- rank_variability(rank_mats)
order_null_path <- file.path(out, "data", "Fig4_PanelD_model_order_null_draws.tsv")
if (file.exists(order_null_path)) {
  order_null <- read.delim(order_null_path, check.names = FALSE)[[1]]
} else {
  order_null <- numeric(50000L)
  for (i in seq_along(order_null)) {
    null_mats <- lapply(rank_mats, function(mat) t(apply(mat, 1L, sample, size = ncol(mat), replace = FALSE)))
    order_null[i] <- rank_variability(null_mats)
  }
}
order_p <- (1 + sum(order_null <= observed_order)) / (length(order_null) + 1)
stopifnot(abs(observed_order - 136.458333333333) < 1e-9)
stopifnot(abs(median(order_null) - 180.708333333333) < 1e-9)
stopifnot(abs(order_p - 0.00217995640087198) < 1e-12)
write.table(data.frame(null_rank_trajectory_variability = order_null), order_null_path, sep = "\t", quote = FALSE, row.names = FALSE)

null_stats <- data.frame(
  panel = c("B", "C", "D"),
  effect = c("Between-Constraint heterogeneity", "Model × Constraint interaction", "Overall model-order component"),
  observed = c(knowledge_result$observed_constraint_mean_SS, interaction_result$observed_interaction_SS, observed_order),
  null_median = c(knowledge_result$null_median, interaction_result$null_median_interaction_SS, median(order_null)),
  null_q025 = c(knowledge_result$null_q025, interaction_result$null_q025_interaction_SS, unname(quantile(order_null, 0.025))),
  null_q975 = c(knowledge_result$null_q975, interaction_result$null_q975_interaction_SS, unname(quantile(order_null, 0.975))),
  p_value = c(knowledge_result$p_value, interaction_result$p_value, order_p),
  stringsAsFactors = FALSE
)
write.table(null_stats, file.path(out, "data", "Fig4_PanelsBCD_null_statistics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

draw_null <- function(null, observed, title, subtitle, xlab, p_label, direction = "right") {
  d <- data.frame(value = null)
  q <- quantile(null, c(0.025, 0.975))
  x_range <- range(c(null, observed))
  span <- diff(x_range)
  p <- ggplot(d, aes(x = value)) +
    annotate("rect", xmin = q[[1]], xmax = q[[2]], ymin = -Inf, ymax = Inf, fill = pal[["pale_blue"]], alpha = 0.72) +
    geom_histogram(aes(y = after_stat(density)), bins = 70, fill = pal[["null"]], colour = "white", linewidth = 0.12) +
    geom_density(colour = pal[["null_dark"]], linewidth = 0.55, adjust = 1.05) +
    geom_vline(xintercept = observed, colour = pal[["observed"]], linewidth = 0.9) +
    annotate("text", x = observed, y = Inf, label = sprintf("Observed\n%.3f", observed), vjust = 1.15,
             hjust = ifelse(direction == "right", 1.05, -0.05), family = "Arial", fontface = "bold", size = 2.7, colour = pal[["observed"]]) +
    annotate("text", x = ifelse(direction == "right", x_range[2], x_range[1]), y = -Inf,
             label = p_label, hjust = ifelse(direction == "right", 1, 0), vjust = -0.85,
             family = "Arial", fontface = "bold", size = 2.65, colour = pal[["ink"]], parse = TRUE) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.05))) +
    labs(title = title, subtitle = subtitle, x = xlab, y = "Null density") +
    theme_pub(7.0) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          plot.title = element_text(size = 7.8), plot.subtitle = element_text(size = 6.3),
          plot.margin = margin(5, 7, 7, 5))
  p
}

pB <- draw_null(knowledge_null, knowledge_result$observed_constraint_mean_SS,
  "Between-Constraint heterogeneity", "Do biological knowledge measurements differ?",
  "Between-Constraint statistic", "italic(P) < 10^{-6}", "right")
pC <- draw_null(interaction_null, interaction_result$observed_interaction_SS,
  "Model-by-Constraint interaction", "Does relative model advantage depend on knowledge?",
  "Interaction statistic", "italic(P) == 1.50 %*% 10^{-5}", "right")
pD <- draw_null(order_null, observed_order,
  "Overall model-order component", "Are model ranks more coherent than random?",
  "Rank-trajectory variability", "italic(P) == 0.00218", "left")
save_panel(pB, "Fig4_PanelB_knowledge_null", 63, 52)
save_panel(pC, "Fig4_PanelC_interaction_null", 63, 52)
save_panel(pD, "Fig4_PanelD_model_order_null", 63, 52)

# Panel E: complete rank structure, with prespecified Constraint order retained.
rank_source <- scores[, c("domain", "Constraint", "Short_name", "Model", "F", "rank")]
rank_source$domain <- factor(rank_source$domain, levels = domain_order)
rank_source$Constraint <- factor(rank_source$Constraint, levels = constraint_order)
rank_source$constraint_label <- gsub("\n", " ", unname(constraint_label[as.character(rank_source$Constraint)]))
rank_source$model_label <- unname(model_label[rank_source$Model])
write.table(rank_source, file.path(out, "data", "Fig4_PanelE_rank_matrix.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

draw_rank_domain <- function(dom) {
  d <- rank_source[rank_source$domain == dom, ]
  cons <- constraint_order[constraint_order %in% as.character(d$Constraint)]
  mods <- rev(model_order[[dom]])
  d$Constraint <- factor(as.character(d$Constraint), levels = cons)
  d$Model <- factor(d$Model, levels = mods)
  ggplot(d, aes(x = Constraint, y = Model, fill = rank)) +
    geom_tile(colour = "white", linewidth = 0.65) +
    geom_text(aes(label = rank), family = "Arial", size = 2.25,
              colour = ifelse(d$rank <= 2, "white", pal[["ink"]]), fontface = "bold") +
    scale_fill_gradientn(colours = c("#244F73", "#5E83A1", "#9FB4C4", "#D4DEE4", "#F0F3F5"),
                         values = scales::rescale(1:5), limits = c(1, 5), guide = "none") +
    scale_x_discrete(labels = function(x) unname(constraint_label[x]), expand = c(0, 0)) +
    scale_y_discrete(labels = function(x) unname(model_label[x]), expand = c(0, 0)) +
    labs(title = dom, x = NULL, y = NULL) +
    coord_fixed(ratio = 0.82, clip = "off") +
    theme_minimal(base_family = "Arial", base_size = 6.3) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(size = 5.35, angle = 38, hjust = 1, vjust = 1, lineheight = 0.9, colour = pal[["ink"]], margin = margin(t = 2)),
          axis.text.y = element_text(size = 5.8, colour = pal[["ink"]]),
          plot.title = element_text(face = "bold", size = 7.2, colour = domain_cols[[dom]], margin = margin(b = 2)),
          plot.margin = margin(2, 4, 5, 2))
}
pE <- (draw_rank_domain("DNA") | draw_rank_domain("RNA")) /
      (draw_rank_domain("Protein sequence") | draw_rank_domain("Protein structure")) +
      plot_annotation(title = "Model ranks change across biological knowledge") &
      theme(plot.title = element_text(family = "Arial", face = "bold", size = 9, colour = pal[["ink"]]))
save_panel(pE, "Fig4_PanelE_rank_matrices", 183, 122)

# Panel F: each pairwise cell is split by direction of interval-supported advantage.
reversals$domain <- factor(reversals$domain, levels = domain_order)
write.table(reversals, file.path(out, "data", "Fig4_PanelF_pairwise_reversals.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

triangle_data <- function(dom) {
  mods <- model_order[[dom]]
  dd <- reversals[reversals$domain == dom, ]
  polys <- list(); labels <- list(); idx <- 1L
  for (i in seq_len(length(mods) - 1L)) for (j in (i + 1L):length(mods)) {
    a <- mods[i]; b <- mods[j]
    row <- dd[(dd$model_a == a & dd$model_b == b) | (dd$model_a == b & dd$model_b == a), ]
    if (nrow(row) != 1L) stop("Missing pair in reversal table")
    if (row$model_a == a) {
      a_wins <- row$a_interval_supported_wins; b_wins <- row$b_interval_supported_wins
    } else {
      a_wins <- row$b_interval_supported_wins; b_wins <- row$a_interval_supported_wins
    }
    x0 <- j; y0 <- length(mods) - i + 1
    polys[[idx]] <- data.frame(domain = dom, pair = paste(a, b), direction = "row model higher",
      x = c(x0 - .46, x0 + .46, x0 - .46), y = c(y0 + .46, y0 + .46, y0 - .46), supported = a_wins > 0)
    idx <- idx + 1L
    polys[[idx]] <- data.frame(domain = dom, pair = paste(a, b), direction = "column model higher",
      x = c(x0 + .46, x0 + .46, x0 - .46), y = c(y0 - .46, y0 + .46, y0 - .46), supported = b_wins > 0)
    idx <- idx + 1L
    labels[[length(labels) + 1L]] <- data.frame(x = x0 - .18, y = y0 + .17, label = ifelse(a_wins > 0, a_wins, ""), direction = "row model higher")
    labels[[length(labels) + 1L]] <- data.frame(x = x0 + .18, y = y0 - .17, label = ifelse(b_wins > 0, b_wins, ""), direction = "column model higher")
  }
  list(poly = do.call(rbind, polys), label = do.call(rbind, labels))
}

draw_reversal_domain <- function(dom) {
  mods <- model_order[[dom]]; td <- triangle_data(dom)
  nrev <- sum(reversals$domain == dom & reversals$ci_supported_bidirectional_reversal)
  td$poly$state <- ifelse(!td$poly$supported, "No supported advantage", td$poly$direction)
  ggplot(td$poly, aes(x = x, y = y, group = interaction(pair, direction), fill = state)) +
    geom_polygon(colour = "white", linewidth = 0.42) +
    geom_text(data = td$label[td$label$label != "", ], aes(x = x, y = y, label = label, colour = direction),
              inherit.aes = FALSE, family = "Arial", fontface = "bold", size = 2.0, show.legend = FALSE) +
    scale_fill_manual(values = c("row model higher" = pal[["teal"]], "column model higher" = pal[["ochre"]], "No supported advantage" = "#E4E8EA"), guide = "none") +
    scale_colour_manual(values = c("row model higher" = "white", "column model higher" = "white"), guide = "none") +
    scale_x_continuous(breaks = 2:5, labels = unname(model_label[mods[2:5]]), limits = c(1.45, 5.55), expand = c(0, 0)) +
    scale_y_continuous(breaks = 2:5, labels = unname(model_label[rev(mods[1:4])]), limits = c(1.45, 5.55), expand = c(0, 0)) +
    coord_fixed(clip = "off") +
    labs(title = dom, subtitle = sprintf("%d/10 bidirectional reversals", nrev), x = NULL, y = NULL) +
    theme_void(base_family = "Arial") +
    theme(axis.text.x = element_text(size = 5.1, angle = 45, hjust = 1, colour = pal[["ink"]]),
          axis.text.y = element_text(size = 5.1, colour = pal[["ink"]]),
          axis.ticks = element_blank(),
          plot.title = element_text(face = "bold", size = 7.2, colour = domain_cols[[dom]], margin = margin(b = 0)),
          plot.subtitle = element_text(size = 5.6, colour = pal[["muted"]], margin = margin(b = 2)),
          plot.margin = margin(2, 5, 8, 6))
}
pF <- (draw_reversal_domain("DNA") | draw_reversal_domain("RNA") |
       draw_reversal_domain("Protein sequence") | draw_reversal_domain("Protein structure")) +
      plot_annotation(
        title = "Bidirectional model advantages occur in 23 of 40 pairs",
        subtitle = "Teal: row model higher; ochre: column model higher. Numbers count supporting Constraints."
      ) & theme(plot.title = element_text(family = "Arial", face = "bold", size = 9, colour = pal[["ink"]]),
                plot.subtitle = element_text(family = "Arial", size = 6.2, colour = pal[["muted"]]))
save_panel(pF, "Fig4_PanelF_reversal_matrices", 183, 65)

# Panels G-H: representative signed differences. Bounds are conservative
# differences formed from the two formal marginal 95% confidence intervals;
# they are not relabelled as paired-bootstrap intervals.
build_delta <- function(dom, constraints, model_a, model_b) {
  do.call(rbind, lapply(constraints, function(id) {
    a <- scores[scores$Constraint == id & scores$Model == model_a, ]
    b <- scores[scores$Constraint == id & scores$Model == model_b, ]
    stopifnot(nrow(a) == 1L, nrow(b) == 1L)
    data.frame(
      domain = dom, Constraint = id, label = gsub("\n", " ", unname(constraint_label[id])),
      model_a = model_a, model_b = model_b,
      F_a = a$F, F_b = b$F, delta_F = a$F - b$F,
      conservative_lower = a$F_ci_lower - b$F_ci_upper,
      conservative_upper = a$F_ci_upper - b$F_ci_lower,
      stringsAsFactors = FALSE
    )
  }))
}
delta_dna <- build_delta("DNA", c("DNA01", "DNA02-02", "DNA06"), "Caduceus-PS", "DNABERT-2")
delta_st <- build_delta("Protein structure", c("ST02", "ST04", "ST05"), "ESM-IF1", "SSLM")
write.table(delta_dna, file.path(out, "data", "Fig4_PanelG_DNA_reversal.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(delta_st, file.path(out, "data", "Fig4_PanelH_structure_reversal.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

draw_delta <- function(d, title, left_label, right_label, col) {
  d$plot_label <- unname(constraint_label[d$Constraint])
  d$plot_label <- factor(d$plot_label, levels = rev(d$plot_label))
  lim <- max(abs(c(d$conservative_lower, d$conservative_upper))) * 1.22
  ggplot(d, aes(x = delta_F, y = plot_label)) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#F5F1EB") +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#EDF3F5") +
    geom_vline(xintercept = 0, linewidth = 0.45, colour = pal[["muted"]]) +
    geom_segment(aes(x = conservative_lower, xend = conservative_upper, yend = plot_label), linewidth = 0.75, colour = col) +
    geom_point(shape = 21, size = 3.1, stroke = 0.5, fill = col, colour = "white") +
    geom_text(aes(label = sprintf("%+.2f", delta_F)), nudge_y = 0.25, family = "Arial", fontface = "bold", size = 2.45, colour = pal[["ink"]]) +
    annotate("text", x = -lim * .95, y = Inf, label = left_label, hjust = 0, vjust = 1.35, family = "Arial", size = 2.35, colour = pal[["muted"]]) +
    annotate("text", x = lim * .95, y = Inf, label = right_label, hjust = 1, vjust = 1.35, family = "Arial", size = 2.35, colour = pal[["muted"]]) +
    scale_x_continuous(limits = c(-lim, lim), expand = c(0, 0)) +
    labs(title = title, x = expression(Delta*F), y = NULL) +
    theme_pub(7.0) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(size = 6.2),
          panel.grid.major.y = element_line(colour = "white", linewidth = 0.55),
          plot.title = element_text(size = 7.8), plot.margin = margin(6, 7, 5, 5))
}
pG <- draw_delta(delta_dna, "Caduceus-PS versus DNABERT-2", "DNABERT-2 higher", "Caduceus-PS higher", domain_cols[["DNA"]])
pH <- draw_delta(delta_st, "ESM-IF1 versus SSLM", "SSLM higher", "ESM-IF1 higher", domain_cols[["Protein structure"]])
save_panel(pG, "Fig4_PanelG_DNA_reversal", 86, 57)
save_panel(pH, "Fig4_PanelH_structure_reversal", 86, 57)

cat("Created Figure 4 candidate panels A-H in", file.path(out, "fig"), "\n")
