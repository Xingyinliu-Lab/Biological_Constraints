# RNA05 released RefSet

This directory contains the complete, minimal and ready-to-Forward RNA05 population.

| File | Role |
|---|---|
| `RNA05_06_unique_inputs.fasta` | 105 unique RNA sequences; IDs are sequence SHA256 values |
| `FORMAL_8BIN_ENDPOINT_LEDGER.tsv` | 6,206 formal position endpoints with block, state, parent and dependency-cluster mappings |
| `RNA05_PRIMARY_MATCHED_BLOCKS_8BINS.tsv` | 1,514 complete matched blocks and their frozen identities |
| `RNA05_PRIMARY_SOURCE_TRUTH_SUMMARY.json` | source-truth and population materialization summary |

Only the complete 8-bin population is formal. The ViennaRNA-compatible subset belongs to the comparator scope and is not the model-scoring RefSet.

Encode every required sequence/position and preserve IDs and zero-based biological coordinates. Model token indices are not RefSet coordinates; tokenizer span mapping and position-vector readout must be documented as Forward provenance. The full mapping requirements, Forward-output contract and scorer command are in the [RNA05 README](../README.md).
