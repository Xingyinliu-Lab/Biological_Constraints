# PS09 formal RefSet

This directory contains the complete, ready-to-Forward PS09 evaluation population. It contains no historical candidate population or RefSet-construction workspace.

## Files used for model evaluation

| File | Role |
|---|---|
| `SELECTED_PROTEINS.faa` | The 81 protein sequences the user encodes. |
| `FORMAL_BLOCK_MEMBERSHIP.tsv.gz` | The 15,708 formal residue endpoints; provides protein accession, one-based position, amino acid, block, state and cluster. |
| `MATCHED_BLOCKS.tsv` | Frozen matched-block inventory and endpoint counts. |
| `MMSEQS2_CLUSTERS.tsv` | Frozen homology-cluster membership used for dependency-aware aggregation. |
| `PS09_FROZEN_CLUSTER_FOLD_ASSIGNMENT.tsv` | Frozen five-fold assignment used by the public readout. |

`DISPROT_EXPLICIT_TRANSITION_REGIONS.tsv.gz` and `PS09_USED_COORDINATE_CACHE_INDEX.tsv` preserve the released truth-source and coordinate evidence. They are not additional sequences that the user must encode.

## Forward handoff

Encode all records in `SELECTED_PROTEINS.faa`, then extract the representation at each membership row's `uniprot_position_1based`. Preserve the formal endpoint ID defined in the [constraint guide](../README.md). Do not change membership, labels, blocks, clusters or folds.

The membership table supplies the frozen biological protein position; model-specific token-to-residue mapping, windowing and residue pooling are the submitter's responsibility. They are not a second scorer input or a tokenizer-validation requirement.

The user supplies model outputs only. The PS09 public scorer performs validation, cross-fit, aggregation and uncertainty calculation.
