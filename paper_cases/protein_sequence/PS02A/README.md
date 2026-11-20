# PS02A paper implementation cases

These five cases record how the paper models produced PS02A representations and
how those representations were submitted to the same public Constraint scorer
available to external users.

| Model | Forward implementation | Public-scorer result |
|---|---|---|---|
| ESM-2 650M | `models/esm2-650m/forward/` | `models/esm2-650m/output/` |
| ProtT5 | `models/prott5/forward/` | `models/prott5/output/` |
| ESM-C 600M | `models/esmc-600m/forward/` | `models/esmc-600m/output/` |
| CARP-640M | `models/carp-640m/forward/` | `models/carp-640m/output/` |
| ESM3 sequence track | `models/esm3-seq/forward/` | `models/esm3-seq/output/` |

The saved Forward representations are not duplicated in Git. The supplied Forward
implementation and model resource record are the reproduction route. Every
`scoring_call/COMMANDS.md` invokes
`constraints/protein_sequence/PS02A/scoring/score_ps02a.py` and contains no
case-local scoring implementation.

Comparator execution and results are intentionally absent here. They belong to
the PS02A Constraint release.
