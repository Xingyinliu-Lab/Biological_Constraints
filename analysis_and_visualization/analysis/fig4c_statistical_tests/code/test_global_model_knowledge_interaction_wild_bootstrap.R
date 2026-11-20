suppressPackageStartupMessages(library(utils))

# Fixed-panel, domain-stratified test of model-by-constraint interaction.
# H0: within each domain, Constraint changes a shared baseline only; each model's
# relative F advantage is invariant across the prespecified constraints.
# Calibration: cluster wild bootstrap at Constraint level, carrying all five
# model residuals together. This is conditional on the fixed 22-constraint panel.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript")
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), ".."), mustWork = TRUE)
input_file <- file.path(root, "data/Supplementary_Table_S4_Complete_Model_Performance.tsv")
out_dir <- file.path(root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

domain_map <- c(
  DNA01 = "DNA", `DNA02-02` = "DNA", DNA03 = "DNA", DNA06 = "DNA",
  `RNA01-02` = "RNA", RNA04 = "RNA", RNA05 = "RNA", RNA06 = "RNA", RNA09 = "RNA", RNA10 = "RNA",
  PS01 = "Protein sequence", PS02A = "Protein sequence", PS02B = "Protein sequence", PS04 = "Protein sequence",
  PS05 = "Protein sequence", PS08 = "Protein sequence", PS09 = "Protein sequence", PS10 = "Protein sequence",
  ST01 = "Protein structure", ST02 = "Protein structure", ST04 = "Protein structure", ST05 = "Protein structure"
)
domain_levels <- c("DNA", "RNA", "Protein sequence", "Protein structure")

scores <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
scores$domain <- factor(unname(domain_map[scores$Constraint]), levels = domain_levels)
if (anyNA(scores$domain) || any(table(scores$Constraint) != 5L)) stop("Incomplete model-by-constraint panel")

fit_domain <- function(dat, response = dat$F) {
  dat <- dat[order(dat$Constraint, dat$Model), ]
  dat$.response <- response
  fit <- lm(.response ~ Constraint + Model, data = dat)
  x <- model.matrix(fit)
  q <- qr.Q(qr(x), complete = FALSE)
  residual_projector <- diag(nrow(x)) - tcrossprod(q)
  residual <- as.vector(residual_projector %*% dat$.response)
  list(
    data = dat,
    residual = residual,
    residual_projector = residual_projector,
    constraint_index = match(dat$Constraint, unique(dat$Constraint)),
    rss = sum(residual^2)
  )
}

fits <- lapply(split(scores, scores$domain, drop = TRUE), fit_domain)
observed_T <- sum(vapply(fits, `[[`, numeric(1), "rss"))
observed_by_domain <- data.frame(
  domain = names(fits),
  interaction_SS = vapply(fits, `[[`, numeric(1), "rss"),
  constraints = vapply(fits, function(z) length(unique(z$data$Constraint)), integer(1)),
  stringsAsFactors = FALSE
)

# Rademacher wild weights are assigned per Constraint, not per model score.
# The precomputed residual projector is algebraically equivalent to refitting
# F ~ Constraint + Model in every draw, while preserving the same null design.
simulate_null_batch <- function(fits, batch_size) {
  total <- numeric(batch_size)
  for (fit in fits) {
    multipliers <- matrix(sample(c(-1, 1), length(unique(fit$data$Constraint)) * batch_size, replace = TRUE),
                          nrow = length(unique(fit$data$Constraint)), ncol = batch_size)
    perturbed_residual <- fit$residual * multipliers[fit$constraint_index, , drop = FALSE]
    projected <- fit$residual_projector %*% perturbed_residual
    total <- total + colSums(projected^2)
  }
  total
}

set.seed(20261022)
B <- 200000L
null_T <- numeric(B)
batch_size <- 5000L
for (start in seq.int(1L, B, by = batch_size)) {
  end <- min(B, start + batch_size - 1L)
  null_T[start:end] <- simulate_null_batch(fits, end - start + 1L)
}
p_value <- (1 + sum(null_T >= observed_T)) / (B + 1)

result <- data.frame(
  constraints = length(unique(scores$Constraint)),
  formal_F_cells = nrow(scores),
  domains = length(domain_levels),
  bootstrap_replicates = B,
  observed_interaction_SS = observed_T,
  null_median_interaction_SS = median(null_T),
  null_q025_interaction_SS = unname(quantile(null_T, 0.025)),
  null_q975_interaction_SS = unname(quantile(null_T, 0.975)),
  exceedances = sum(null_T >= observed_T),
  p_value = p_value,
  stringsAsFactors = FALSE
)
write.table(result, file.path(out_dir, "global_model_knowledge_interaction_wild_bootstrap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(observed_by_domain, file.path(out_dir, "global_model_knowledge_interaction_by_domain.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(null_interaction_SS = null_T), file.path(out_dir, "global_model_knowledge_interaction_null_draws.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

p_text <- if (sum(null_T >= observed_T) == 0L) sprintf("P < %.2g", 1 / (B + 1)) else sprintf("P = %.5g", p_value)
writeLines(c(
  "# Global model-knowledge interaction test",
  "",
  "## Null hypothesis",
  "Within each biological domain, a Constraint may shift the common F baseline, but the relative F advantage of each model is invariant across Constraints.",
  "",
  "## Statistic and calibration",
  "The statistic is the sum, across the four domains, of the residual sum of squares from F ~ Constraint + Model. It quantifies model-by-Constraint interaction after removing Constraint and domain-specific model main effects. A 200,000-draw Rademacher cluster wild bootstrap assigns one multiplier to all five residuals within each Constraint, preserving within-Constraint shared structure under the null.",
  "",
  sprintf("Observed interaction SS = %.8f; null median = %.8f; %s.", observed_T, median(null_T), p_text),
  "",
  "## Interpretation boundary",
  "This is a conditional inference for the fixed, prespecified 22-constraint panel. It tests constraint-specific variation in relative model advantage, not intrinsic biological difficulty and not a superpopulation of all possible biological relations.",
  "",
  "No Forward pass, primary Scorer run, or formal-result modification was performed."
), file.path(out_dir, "GLOBAL_MODEL_KNOWLEDGE_INTERACTION_SUMMARY.md"))
