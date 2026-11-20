# DNA02-02 released RefSet

This directory contains the complete, minimal and directly forwardable DNA02-02 evaluation population. Candidate construction tables and historical materialization intermediates are intentionally excluded from the public Constraint.

## Runtime assets

| File | Role |
|---|---|
| `DNA02_02_SEQUENCE_OBJECTS.fasta` | The 54,020 DNA objects that every evaluated model must encode. |
| `SEQUENCE_OBJECT.tsv` | Stable object IDs, lengths, alphabet and sequence SHA256 values. |
| `CANONICAL_LADDER.tsv` | The 46,720 frozen 1/1/2 edit ladders and all relation/group identities required for scoring. |
| `DNA02_02_OUTER_FOLDS_v0.1.0.tsv` | Frozen outer-cluster fold assignment. |
| `DNA02_02_TAU_CONTROLS_v0.1.1.tsv` | Frozen repeat-calibration membership: 32 ladders in each formal length stratum. |

Use the Forward and representation contracts in [../README.md](../README.md). Do not substitute earlier candidate, small or construction populations.
