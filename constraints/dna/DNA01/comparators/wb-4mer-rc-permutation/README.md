# DNA01 exact 4-mer RC permutation

This is the registered analytic white-box reference. It converts each frozen sequence to overlapping 4-mer frequencies in lexicographic `A,C,G,T` channel order and applies the exact channel permutation `k → RC(k)`. No checkpoint or fitted operator is used.

The comparator reuses DNA01's shared RefSet validation, smooth score semantics, hierarchy, bootstrap and publication functions. Because the permutation realizes reverse complement exactly, its canonical result is `F=1` with 95% CI `[1,1]`. It is a constructive ceiling, not evidence that a learned model discovered the relation.

## Reproduce the reference

From the repository root:

```bash
python constraints/dna/DNA01/comparators/wb-4mer-rc-permutation/code/run_comparator.py \
  --relations constraints/dna/DNA01/refset/DNA01_RELATION_TABLE.tsv \
  --anchors constraints/dna/DNA01/refset/DNA01_ANCHOR_TABLE.tsv \
  --objects constraints/dna/DNA01/refset/DNA01_OBJECT_TABLE.tsv \
  --fold-assignment constraints/dna/DNA01/refset/DNA01_FOLD_ASSIGNMENT_TABLE.tsv \
  --dependency-clusters constraints/dna/DNA01/refset/DNA01_DEPENDENCY_CLUSTER_TABLE.tsv \
  --repeat-control-panel constraints/dna/DNA01/refset/DNA01_REPEAT_CONTROL_PANEL.tsv \
  --fasta constraints/dna/DNA01/refset/records.fasta \
  --output work/DNA01/wb-4mer-rc-permutation
```

Compare the generated `FORMAL_RESULT.json` and `INDIVIDUAL_SCORES.tsv` with the canonical v0.5.4 files in [`results/`](results/). Method and fixed settings are in [`method/METHOD.json`](method/METHOD.json).

The canonical individual-score table contains all 1,460 formal relations and uses the same evaluation IDs and columns as model scoring, enabling relation-level comparison with the built-in reference.
