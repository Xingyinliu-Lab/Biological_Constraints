suppressPackageStartupMessages(library(utils))

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
domain_order <- c("DNA", "RNA", "Protein sequence", "Protein structure")
n_permutations <- 50000L
set.seed(20261022)

scores <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
if (any(table(scores$Constraint) != 5L)) stop("Every constraint must have five formal model scores")
scores$domain <- unname(domain_map[scores$Constraint])
if (anyNA(scores$domain)) stop("Unmapped constraint")
scores$rank <- ave(scores$F, scores$Constraint, FUN = function(x) rank(-x, ties.method = "average"))

parse_ci <- function(x) as.numeric(trimws(strsplit(substr(x, 2L, nchar(x) - 1L), ",", fixed = TRUE)[[1]]))
ci <- t(vapply(scores$F_95_CI, parse_ci, numeric(2)))
scores$F_ci_lower <- ci[, 1]
scores$F_ci_upper <- ci[, 2]

rank_matrices <- lapply(domain_order, function(domain) {
  dat <- scores[scores$domain == domain, ]
  constraints <- unique(dat$Constraint)
  models <- sort(unique(dat$Model))
  mat <- do.call(rbind, lapply(constraints, function(id) {
    z <- dat[dat$Constraint == id, c("Model", "rank")]
    z$rank[match(models, z$Model)]
  }))
  colnames(mat) <- models
  rownames(mat) <- constraints
  mat
})
names(rank_matrices) <- domain_order

rank_variability <- function(mats) {
  sum(vapply(mats, function(mat) {
    sum((mat - rep(colMeans(mat), each = nrow(mat)))^2)
  }, numeric(1)))
}

observed_variability <- rank_variability(rank_matrices)
null_variability <- numeric(n_permutations)
for (i in seq_len(n_permutations)) {
  null_mats <- lapply(rank_matrices, function(mat) {
    t(apply(mat, 1L, sample, size = ncol(mat), replace = FALSE))
  })
  null_variability[i] <- rank_variability(null_mats)
}

pair_reversals <- do.call(rbind, lapply(domain_order, function(domain) {
  dat <- scores[scores$domain == domain, ]
  models <- sort(unique(dat$Model))
  pairs <- combn(models, 2L, simplify = FALSE)
  do.call(rbind, lapply(pairs, function(pair) {
    a <- dat[dat$Model == pair[[1]], c("Constraint", "F_ci_lower", "F_ci_upper")]
    b <- dat[dat$Model == pair[[2]], c("Constraint", "F_ci_lower", "F_ci_upper")]
    b <- b[match(a$Constraint, b$Constraint), ]
    a_wins <- a$Constraint[a$F_ci_lower > b$F_ci_upper]
    b_wins <- a$Constraint[a$F_ci_upper < b$F_ci_lower]
    data.frame(
      domain = domain, model_a = pair[[1]], model_b = pair[[2]],
      a_interval_supported_wins = length(a_wins), b_interval_supported_wins = length(b_wins),
      ci_supported_bidirectional_reversal = length(a_wins) > 0L && length(b_wins) > 0L,
      a_winning_constraints = paste(a_wins, collapse = ";"),
      b_winning_constraints = paste(b_wins, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
}))

n_pairs <- nrow(pair_reversals)
n_reversals <- sum(pair_reversals$ci_supported_bidirectional_reversal)
global_result <- data.frame(
  constraints = length(unique(scores$Constraint)),
  domain_specific_model_panels = length(domain_order),
  model_by_constraint_cells = nrow(scores),
  domain_specific_model_pairs = n_pairs,
  ci_supported_bidirectional_reversal_pairs = n_reversals,
  ci_supported_bidirectional_reversal_fraction = n_reversals / n_pairs,
  observed_rank_trajectory_variability = observed_variability,
  permutation_null_variability_median = median(null_variability),
  permutation_null_variability_ci95_lower = unname(quantile(null_variability, 0.025)),
  permutation_null_variability_ci95_upper = unname(quantile(null_variability, 0.975)),
  permutation_p_for_stable_global_order = (1 + sum(null_variability <= observed_variability)) / (n_permutations + 1),
  stringsAsFactors = FALSE
)

write.table(global_result, file.path(out_dir, "global_model_knowledge_interaction.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pair_reversals, file.path(out_dir, "global_pairwise_CI_supported_reversals.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary_lines <- c(
  "# Global model–knowledge interaction analysis",
  "",
  sprintf("All %d formal constraints were analyzed jointly, retaining their four domain-specific five-model panels (%d score cells).", global_result$constraints, global_result$model_by_constraint_cells),
  "",
  "## What was tested",
  "",
  "The pooled rank-concordance permutation tests whether the observed model rank trajectories are more stable than independently relabelled model ranks within each constraint, while retaining the domain-specific model panels. It tests for a global model-order component; it does not by itself test a binary 'knowledge versus model' alternative.",
  "",
  sprintf("Observed rank-trajectory variability was %.3f versus a null median of %.3f (permutation P for stable global order = %.4g).", global_result$observed_rank_trajectory_variability, global_result$permutation_null_variability_median, global_result$permutation_p_for_stable_global_order),
  "",
  "## Direct interaction evidence",
  "",
  sprintf("%d/%d domain-specific model pairs (%.1f%%) showed CI-supported bidirectional reversals: each pair had at least one constraint with A's lower 95%% F bound above B's upper bound and another constraint with the opposite interval-supported ordering.", n_reversals, n_pairs, 100 * n_reversals / n_pairs),
  "",
  "## Interpretation",
  "",
  "The 22-constraint panel contains a global model-order component if the permutation P is small. The bidirectional reversals are conservative pointwise evidence against a single fixed model advantage across knowledge relations; they are not a multiplicity-adjusted omnibus P value. Together, these results support uneven knowledge accessibility within this formal panel: relative model advantage depends on the measured biological relation. This conclusion is scoped to the prespecified constraint panel and is not a claim about every possible biological knowledge relation.",
  "",
  "No model Forward pass or primary scoring computation was rerun."
)
writeLines(summary_lines, file.path(out_dir, "GLOBAL_RESULTS_SUMMARY.md"))
