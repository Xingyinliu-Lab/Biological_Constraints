# DNA06 paper cases

These five cases use the current validated paper Forward implementations. Each
implementation produces the single formal codon-local representation as
`row_ids + embeddings` and submits it to the same public DNA06 scorer.

The public-scorer parity run completed for all five models. Differences from
the registered values are at most floating-point serialization noise
(`3.2e-15`) and are below the frozen numeric comparison tolerance.

Each `models/<model_id>/output/` directory contains:

- `INDIVIDUAL_SCORES.tsv` and `SITE_SCORE_LEDGER.tsv`: ID-level measurements
  and credits;
- `PARENT_SCORE_LEDGER.tsv` and `CLUSTER_SCORE_LEDGER.tsv`: F-only hierarchy ledgers;
- `FORMAL_RESULT.json`: F, its percentile confidence interval, coverage, and QC;
- `OUTPUT_MANIFEST.json` and `RUN_COMPLETE`: output identity and completion.

Forward implementations remain in the case; all
constraint-specific measurement and scoring is provided only by
`constraints/dna/DNA06/scoring/score_dna06.py`.
