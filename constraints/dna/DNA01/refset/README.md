# DNA01 frozen RefSet

This directory is the ready-to-use DNA01 population. It contains the complete formal object universe, relationships, folds, dependency closure, repeat-control plan and source sequence payload required by the public scorer.

| File | Rows / records | Role |
|---|---:|---|
| `DNA01_OBJECT_TABLE.tsv` | 2,920 | Formal model inputs and sequence identity. |
| `DNA01_RELATION_TABLE.tsv` | 1,460 | Paired FWD/RC views and frozen fold/cluster membership. |
| `DNA01_ANCHOR_TABLE.tsv` | 11,680 | Eight same-view, dependency-isolated anchor targets per relation. |
| `DNA01_FOLD_ASSIGNMENT_TABLE.tsv` | 278 | Dependency-cluster to outer-fold assignment. |
| `DNA01_DEPENDENCY_CLUSTER_TABLE.tsv` | 280 | Formal parent-locus membership of the 278 dependency clusters. |
| `DNA01_REPEAT_CONTROL_PANEL.tsv` | 2,920 | Complete deterministic repeat-Forward population. |
| `records.fasta` | 2,920 | Ready-to-use model Forward input; exactly the formal object universe. |

## Construct model inputs

Use `records.fasta` directly. Submit every record to the model and preserve its FASTA identifier in the Forward output. `DNA01_OBJECT_TABLE.tsv` provides the orientation, length and sequence identity for inspection.

Do not regenerate, filter or relabel the population. Model context limitations do not authorize removing objects. After Forward, follow the NPZ contract and scoring command in [`../README.md`](../README.md).
