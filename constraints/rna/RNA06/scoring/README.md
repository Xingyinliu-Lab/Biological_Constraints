# RNA06 public scoring

`score_rna06.py` is the only formal RNA06 scoring implementation. It accepts nucleotide-position representations, validates coverage and coordinates, computes all frozen constraint-specific features and cross-fit predictions, and writes the complete F-only score chain and F confidence interval.

Formal fixed settings:

- feature: concatenated vector sum and Hadamard product;
- outer split: leave one dependency cluster out;
- inner split: four deterministic folds;
- lambda grid: `1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100`;
- bootstrap: outer dependency clusters, PCG64 seed `20261003`, 20,000 draws, linear percentiles.

Before fitting, the scorer validates FASTA/object identity, coordinates, reciprocal edges, frozen positive-pair membership, candidate-to-endpoint consistency, hierarchy identities, population counts and submitted numerical types. `CANDIDATE_SCORES.tsv` preserves every nonpartner comparison before endpoint and higher-level aggregation.

The full primary RefSet is loaded by default. Explicit RefSet path arguments remain available only for compatible research use and are recorded in the result identity. They do not create another formal population.

The scorer validates the frozen biological sequence/position targets and the distinct anchor, true-partner and nonpartner associations. Model-specific token-to-nucleotide mapping is performed and documented by the user during Forward; it is not implemented by the scorer.
