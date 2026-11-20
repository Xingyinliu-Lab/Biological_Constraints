# esm3-struct

This directory owns the shared model identity and executed-environment record for `esm3-struct`. Constraint-specific input construction, Forward execution, coordinate handling and representation export belong to the corresponding implementation under `paper_cases/`. Formal scoring belongs to the selected Constraint.

## Environment and checkpoint

| File | Purpose |
|---|---|
| [`ENVIRONMENT.md`](ENVIRONMENT.md) | Human-readable record of the environment used for the formal runs |
| [`ENVIRONMENT.json`](ENVIRONMENT.json) | Machine-readable executed package inventory |
| [`checkpoint_identity.json`](checkpoint_identity.json) | Official source, checkpoint/revision identities, SHA256 evidence, download and loading method |

Recorded model identity:

- Official source: https://huggingface.co/EvolutionaryScale/esm3-sm-open-v1
- Evidence status: `EXECUTED_MULTI_FILE_CHECKPOINT_IDENTITY_RECORDED`
- Upstream revision: `NOT_RECORDED_IN_FORMAL_EXECUTION`
- Executed checkpoint SHA256: `5ead5a135c658068db6a4f1b933e72d6110992c4668822e1c0e2dcc53e38acd9`, `467acbaee703ba3ccde6e75241a912a316952e5ff071355f85c1d33c68704f40`, `752362cd8ed1e7826ac366197f3629c47b9a3996457defff30a922b5fde01b6f`
- Loading interface: `ESM3.from_pretrained('esm3_sm_open_v1', device=device).eval(); experimental structure track supplied, native sequence interior masked, structure-token/backbone-coordinate residue states mean pooled`

Use the exact checkpoint and configuration identity recorded above when reproducing a paper result. A newer checkpoint, tokenizer, representation layer or pooling rule is a different execution and must not be presented as reproduction of a registered result. Historical local paths in the identity record are provenance only; users do not need to reproduce the author's directory layout.

## Constraint-specific Forward records

Constraint-specific RefSet handling, Forward code, representation export and scorer calls are recorded in the relevant directory under `paper_cases/<domain>/<constraint_id>/models/<model_id>/`. This shared model resource does not duplicate those commands or scoring outputs.

Use the selected Constraint README as the authority for its current Forward-output contract and public scorer. Historical local paths in a paper case are provenance only; they do not define the public filesystem layout.

## Forward-output handoff

Each paper case records:

1. the released RefSet objects consumed by the Forward path;
2. the executed extraction scripts;
3. the generated representation contract and retained file identity;
4. the call to the unique public Constraint scorer.

Formal paper outputs are stored with the corresponding paper case. They are not model-loading outputs and are not duplicated in this shared resource directory.
