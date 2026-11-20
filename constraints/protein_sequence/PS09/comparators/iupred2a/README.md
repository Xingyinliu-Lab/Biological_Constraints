# PS09 comparator: IUPred2A long mode

IUPred2A is the built-in PS09 sequence-based comparator. This directory contains its executable implementation, environment, ID-level outputs and canonical aggregate F result. It does not maintain a second PS09 aggregation algorithm.

## Produce endpoint scores

Obtain the legacy IUPred2A distribution from
[`iupred2a/iupred2a`](https://github.com/iupred2a/iupred2a) at commit
`ab0115546af598c4a85807616620ea462ea95983`. Set `--iupred-dir` to its
legacy directory, with the required data directory present beside the legacy
program as described in [`ENVIRONMENT.md`](ENVIRONMENT.md).

```bash
python constraints/protein_sequence/PS09/comparators/iupred2a/produce_iupred_leaf.py \
  --refset-dir constraints/protein_sequence/PS09/refset \
  --iupred-dir /path/to/iupred2a/legacy \
  --output-dir work/PS09/iupred2a/measurement
```

The producer runs the included legacy IUPred2A implementation in `long` mode and writes the complete endpoint score file expected by PS09. See [`ENVIRONMENT.md`](ENVIRONMENT.md) for the executed software requirements.

## Apply the PS09 formal aggregation

```bash
python constraints/protein_sequence/PS09/scoring/score_ps09_formal.py \
  --binding constraints/protein_sequence/PS09/scoring/PS09_NUMERIC_INPUT_BINDING.json \
  --membership constraints/protein_sequence/PS09/refset/FORMAL_BLOCK_MEMBERSHIP.tsv.gz \
  --blocks constraints/protein_sequence/PS09/refset/MATCHED_BLOCKS.tsv \
  --clusters constraints/protein_sequence/PS09/refset/MMSEQS2_CLUSTERS.tsv \
  --fasta constraints/protein_sequence/PS09/refset/SELECTED_PROTEINS.faa \
  --folds constraints/protein_sequence/PS09/refset/PS09_FROZEN_CLUSTER_FOLD_ASSIGNMENT.tsv \
  --leaf-scores work/PS09/iupred2a/measurement/FORMAL_LEAF_SCORES.tsv \
  --out work/PS09/iupred2a/score/SCORE_SUMMARY.json \
  --individual-scores work/PS09/iupred2a/score/INDIVIDUAL_SCORES.tsv \
  --pair-credits work/PS09/iupred2a/score/PAIR_CREDIT_LEDGER.tsv.gz \
  --block-scores work/PS09/iupred2a/score/BLOCK_SCORES.tsv \
  --cluster-scores work/PS09/iupred2a/score/CLUSTER_SCORES.tsv
```

## Canonical outputs

[`results/`](results/) contains the executed IUPred2A endpoint profile, formal endpoint scores, pair credits, endpoint/block/cluster score tables and `SCORE_SUMMARY.json`.
