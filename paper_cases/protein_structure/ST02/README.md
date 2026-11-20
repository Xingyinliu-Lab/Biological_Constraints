# ST02 paper implementation cases

Five paper-model cases preserve their public Forward code and formal F-only outputs. They invoke the same public ST02 scorer available to external users. Comparator execution and results remain under the Constraint.

## Cases

| Model | Forward guide | Public scoring call | Registered output |
|---|---|---|---|
| ESM-IF1 | [`models/esm-if1/forward/README.md`](models/esm-if1/forward/README.md) | [`models/esm-if1/scoring_call/COMMANDS.md`](models/esm-if1/scoring_call/COMMANDS.md) | [`models/esm-if1/output/FORMAL_RESULT.json`](models/esm-if1/output/FORMAL_RESULT.json) |
| ESM3-STRUCT | [`models/esm3-struct/forward/README.md`](models/esm3-struct/forward/README.md) | [`models/esm3-struct/scoring_call/COMMANDS.md`](models/esm3-struct/scoring_call/COMMANDS.md) | [`models/esm3-struct/output/FORMAL_RESULT.json`](models/esm3-struct/output/FORMAL_RESULT.json) |
| ProteinMPNN-CA | [`models/proteinmpnn-ca/forward/README.md`](models/proteinmpnn-ca/forward/README.md) | [`models/proteinmpnn-ca/scoring_call/COMMANDS.md`](models/proteinmpnn-ca/scoring_call/COMMANDS.md) | [`models/proteinmpnn-ca/output/FORMAL_RESULT.json`](models/proteinmpnn-ca/output/FORMAL_RESULT.json) |
| ProToken | [`models/protoken/forward/README.md`](models/protoken/forward/README.md) | [`models/protoken/scoring_call/COMMANDS.md`](models/protoken/scoring_call/COMMANDS.md) | [`models/protoken/output/FORMAL_RESULT.json`](models/protoken/output/FORMAL_RESULT.json) |
| SSLM | [`models/sslm/forward/README.md`](models/sslm/forward/README.md) | [`models/sslm/scoring_call/COMMANDS.md`](models/sslm/scoring_call/COMMANDS.md) | [`models/sslm/output/FORMAL_RESULT.json`](models/sslm/output/FORMAL_RESULT.json) |

All five cases consume the 509 released structures in `constraints/protein_structure/ST02/refset/coordinates/` and hand an NPZ containing `object_ids` `(509,)` and `embeddings` `(509,d)` to the unique public ST02 scorer. Model installation and checkpoint identity are recorded once under `paper_cases/model_resources/protein_structure/`.
