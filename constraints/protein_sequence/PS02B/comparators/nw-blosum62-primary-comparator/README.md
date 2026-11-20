# PS02B NW-BLOSUM62 primary comparator

This is the frozen Level-1 raw-sequence explanatory reference for PS02B. It
uses EMBOSS Needle 6.6.0.0 with EBLOSUM62 and fixed gap penalties. Pair scores
do not use HSF/topology labels, profiles, retrieval, taxonomy, structure or
pretrained representations.

The method contract is
[`method/PS02B_NW_BLOSUM62_COMPARATOR_CONTRACT_v0.1.0_FROZEN.json`](method/PS02B_NW_BLOSUM62_COMPARATOR_CONTRACT_v0.1.0_FROZEN.json).
Environment and EMBOSS requirements are in [`ENVIRONMENT.md`](ENVIRONMENT.md).

## Run the complete comparator

From the repository root:

```bash
python constraints/protein_sequence/PS02B/comparators/nw-blosum62-primary-comparator/code/run_ps02b_nw_blosum62_formal.py \
  --needle /path/to/emboss-6.6.0.0/bin/needle \
  --output-dir work/PS02B/nw-blosum62 \
  --workers 8
```

The wrapper validates the frozen RefSet and comparator implementation, records
the SHA256 of the supplied Needle executable, computes pair and relation
measurements, and passes them to the same public PS02B scoring core used for
model representations. EMBOSS Needle 6.6.0.0 and the frozen algorithm
parameters remain required; byte identity with the historical executable is
provenance, not a public admission gate. `--workers` changes parallel execution
only. Existing output directories are not overwritten.

A complete EMBOSS Needle 6.6.0.0 run over all 8,208 formal relations matched
the released comparator at every ID-level relation measurement.

## Fixed inputs and method

- FASTA: `refset/endpoint_sequences.faa`, all 3,937 released Forward objects;
- relations: `refset/FORMAL_RELATION_LEDGER.tsv.gz`, all 8,208 contrasts;
- positive: same HSF as the anchor;
- negative: same topology but a different HSF;
- similarity: `2*s(x,y)/(s(x,x)+s(y,y))`;
- matrix: EBLOSUM62;
- gap open/extend: 10.0/0.5;
- terminal gaps penalized with end open/extend 10.0/0.5;
- relation credit, hierarchy and 20,000-draw bootstrap: the public PS02B score
  core with its frozen 351 groups.

## Outputs

| Path | Contents |
|---|---|
| `measurement/PAIR_SCORE_CACHE.tsv.gz` | raw Needle/self scores and normalized pair similarities |
| `measurement/FORMAL_LEAF_SCORES.tsv.gz` | relation-level positive and negative measurements |
| `measurement/CONTRAST_SCORE_LEDGER.tsv.gz` | relation ID, difference and F_credit (win 1; tie 0; loss -1) |
| `INDIVIDUAL_SCORES.tsv.gz` | public scorer ID-level output and dependency identities |
| `SCORE_SUMMARY.json` | comparator F/F 95% CI, coverage and diagnostics |
| `EXECUTION_RELEASE.json` | method, code, input and output identities |

[`results/`](results/) contains canonical measurements and the result generated
through the consolidated public scoring core. There is no comparator-local copy
of the PS02B score algorithm.
