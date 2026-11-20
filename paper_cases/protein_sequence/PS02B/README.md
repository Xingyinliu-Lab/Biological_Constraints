# PS02B paper implementation cases

These five cases record how the paper models produced PS02B whole-domain
representations and submitted them to the same public Constraint scorer used
by external models.

| Model | Forward implementation | Public-scorer result |
|---|---|---|
| ESM-2 650M | `models/esm2-650m/forward/` | `models/esm2-650m/output/` |
| ProtT5 | `models/prott5/forward/` | `models/prott5/output/` |
| ESM-C 600M | `models/esmc-600m/forward/` | `models/esmc-600m/output/` |
| CARP-640M | `models/carp-640m/forward/` | `models/carp-640m/output/` |
| ESM3 sequence track | `models/esm3-seq/forward/` | `models/esm3-seq/output/` |

The large Forward representation files are not duplicated in Git. Each case
keeps its Forward implementation and registered outputs. Every
`scoring_call/COMMANDS.md` invokes
`constraints/protein_sequence/PS02B/scoring/score_ps02b.py`; no case owns a
second scoring implementation. Comparator execution and results remain inside
the PS02B Constraint release.
