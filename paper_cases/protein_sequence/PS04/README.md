# PS04 paper implementation cases

These five cases preserve the paper Forward implementations and registered
outputs, then submit the same representation contract
to the single public PS04 scorer.

| Model | Forward | Formal output |
|---|---|---|
| ESM-2 650M | `models/esm2-650m/forward/` | `models/esm2-650m/output/` |
| ProtT5 | `models/prott5/forward/` | `models/prott5/output/` |
| ESM-C 600M | `models/esmc-600m/forward/` | `models/esmc-600m/output/` |
| CARP-640M | `models/carp-640m/forward/` | `models/carp-640m/output/` |
| ESM3 sequence track | `models/esm3-seq/forward/` | `models/esm3-seq/output/` |

Representations are intentionally not duplicated in Git. They can be
regenerated with the retained Forward implementations and registered model
resources. Case-local historical scorer copies were not carried forward;
`scoring_call/COMMANDS.md` always invokes the Constraint-owned public entry.
Comparator implementations and results remain inside the PS04 release.

## Public-scorer parity

The registered full Forward outputs were submitted to the public scorer with
the frozen large RefSet. All five aggregate results match their registered
values exactly:

| Model | F | Parity |
|---|---:|---|
| ESM-2 650M | 0.2223862172111226 | PASS |
| ProtT5 | 0.2824558818637282 | PASS |
| ESM-C 600M | 0.0989775258284014 | PASS |
| CARP-640M | 0.2906705030261021 | PASS |
| ESM3 sequence | 0.2699058782037078 | PASS |

Each `output/` directory contains the public `RAW_MEASUREMENTS.tsv.gz`,
`INDIVIDUAL_SCORES.tsv.gz`, `FORMAL_RESULT.json`, and
`OUTPUT_MANIFEST.json`. The case owns no scoring algorithm.
