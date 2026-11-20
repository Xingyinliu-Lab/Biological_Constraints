# DNA06 Standard Genetic Code Table 1 reference

This comparator is the frozen DNA06 white-box reference. It derives the
synonymous versus ordinary-missense relation directly from the released WT/ALT
CDS sequences using Standard Genetic Code Table 1. It uses the same matched-site
population and scoring hierarchy as the public model scorer.

## Formal reference

```bash
python constraints/dna/DNA06/comparators/standard-genetic-code-table1/score_dna06_table1_level0_v010.py \
  --refset-dir constraints/dna/DNA06/refset \
  --output-dir work/DNA06/table1
```

The expected deterministic result is `F=1`. This is a constructive
knowledge ceiling, not evidence of learned representation knowledge.

The command writes `INDIVIDUAL_SCORES.tsv` for all 11,506 variants,
`PAIR_CREDIT_LEDGER.tsv` with every missense–synonymous margin, tie and F credit,
`SITE_SCORE_LEDGER.tsv` for all 3,920 matched sites,
`PARENT_SCORE_LEDGER.tsv` for all 560 parent-CDS contributions,
`CLUSTER_SCORE_LEDGER.tsv` for all 554 dependency clusters, the aggregate
`FORMAL_RESULT.json`, input/output manifests and `RUN_COMPLETE`. The released
canonical copies are in `results/table1/`.

## Raw-codon ridge companion

`results/raw-codon-ridge/FORMAL_RESULT.json` is the separately identified raw-codon
ridge companion result; it is not a second implementation of the DNA06 public
scorer.

```bash
python constraints/dna/DNA06/comparators/standard-genetic-code-table1/score_dna06_raw_codon_ridge_v010.py \
  --refset-dir constraints/dna/DNA06/refset \
  --output-dir work/DNA06/raw_codon_ridge
```

This command writes the same three ID-level ledger layers plus its aggregate
`FORMAL_RESULT.json`, output manifest and completion marker. Canonical copies
are in `results/raw-codon-ridge/`.
