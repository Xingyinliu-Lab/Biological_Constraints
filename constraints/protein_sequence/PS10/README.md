# PS10: amyloidogenic peptide readability

PS10 tests whether peptide representations support discrimination of experimentally amyloidogenic and non-amyloidogenic records under a frozen study-aware 22-fold cross-fit. Each held-out study AUC is transformed to `fold_F = 2AUC - 1`; the primary result is the unweighted macro mean of the 22 fold-level F values, and the bootstrap resamples those 22 fold-level values.

`PRIMARY_ASSIGNMENTS.tsv` defines the formal evaluation. `SENSITIVITY_ASSIGNMENTS.tsv` is retained as a separate sensitivity design and must not replace the primary result.

## 1. RefSet to model Forward input

[`refset/PEPTIDE_REFSET_RECORDS.tsv`](refset/PEPTIDE_REFSET_RECORDS.tsv) is the complete ready-to-Forward population. For every unique `record_id`, encode the canonical peptide sequence in the `normalized_peptide_identity` field using the model being evaluated. Preserve each record ID exactly. This normalized 20-amino-acid sequence, whose length and SHA256 are recorded in the same row, is the sole formal Forward sequence; there is no separate `sequence` field in the released RefSet.

Do not filter records according to model length, tokenizer behavior, representation quality or model output. The released experimental classification, source-study membership and fold roles are fixed scientific inputs. The user does not construct folds or train the PS10 readout.

## 2. Model Forward output contract

Write one NPZ containing exactly:

```text
record_id   one-dimensional unique string array
embedding   finite real numeric matrix, shape (n_records, d), d > 0
```

The ID set must exactly match `PEPTIDE_REFSET_RECORDS.tsv`; row order is arbitrary when IDs and rows remain aligned. Sequence-level pooling is part of the model Forward and must be documented by the user. All PS10-specific scaling, weighting, cross-fit and statistics are performed by the public scorer.

## 3. Run the official primary scorer

Install [`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md), then run:

```bash
python constraints/protein_sequence/PS10/scoring/score_ps10.py \
  --embeddings work/PS10/my_model/embeddings.npz \
  --model-id my_model \
  --output-dir work/PS10/my_model/score
```

The scorer automatically loads the released record table, primary 22-fold assignment and frozen numeric method. It verifies the two formal class labels, canonical peptide length and SHA256, the complete record-by-fold role matrix, held-out study separation, and train/test class closure. The formal CLI does not expose alternative folds, regularization, weighting or bootstrap parameters.

## 4. Official outputs

- `INDIVIDUAL_SCORES.tsv`: each held-out `record_id`, fold, label, source-study identity and raw logistic decision score.
- `FOLD_SCORES.tsv`: train/test/purged counts, raw held-out AUC evidence and formal `fold_F = 2AUC - 1` for each of 22 held-out study folds.
- `SCORE_SUMMARY.json`: F, F 95% CI, coverage, fixed settings and input identities.
- `BOOTSTRAP_DRAWS.float64le.bin`: 20,000 fold-level bootstrap draws.
- `OUTPUT_MANIFEST.json`: output sizes and SHA256 identities.

The formal reporting score F is the unweighted macro mean of the 22 values `fold_F = 2AUC - 1`. The frozen readout is logistic regression with `C=1`, `lbfgs`, maximum 1,000 iterations, training-only unweighted mean/SD scaling and released raw source-study weights; those weights are not fold-wise normalized. The F 95% interval resamples the 22 fold-level F values using PCG64 seed `20261008`, 20,000 draws and linear percentiles.

`formal_scoring_protocol=true` means the frozen primary method completed. It does not automatically register a third-party result as a paper result.

## 5. Comparator

The TANGO 2.3.1 comparator is fully contained in [`comparators/tango-2.3.1/`](comparators/tango-2.3.1/), including executable code, method identity, raw record scores, ID-level held-out outputs and the canonical summary. It uses the same primary fold assignment and shared fold-level aggregation.

Paper model Forward implementations and registered FM outputs are kept outside this independent Constraint package; they call this same public scorer.
