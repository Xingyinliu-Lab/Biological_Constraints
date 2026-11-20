# sslm

This directory owns the shared model identity and executed-environment record for `sslm`. Constraint-specific input construction, Forward execution, coordinate handling and representation export belong to the corresponding implementation under `paper_cases/`. Formal scoring belongs to the selected Constraint.

## Environment and checkpoint

| File | Purpose |
|---|---|
| [`ENVIRONMENT.md`](ENVIRONMENT.md) | Human-readable record of the environment used for the formal runs |
| [`ENVIRONMENT.json`](ENVIRONMENT.json) | Machine-readable executed package inventory |
| [`checkpoint_identity.json`](checkpoint_identity.json) | Official source, checkpoint/revision identities, SHA256 evidence, download and loading method |

Recorded model identity:

- Official source: https://github.com/wr1102/CPE-Pro
- Evidence status: `EXECUTED_CHECKPOINT_HASH_RECORDED`
- Upstream revision: `NOT_RECORDED_IN_FORMAL_EXECUTION`
- Executed checkpoint SHA256: `1daf10751fa0469e1ab24481a04b8a1e6809c32d05b219156247cb880d54c198`
- Loading interface: `EsmTokenizer.from_pretrained(model_dir); EsmConfig.from_pretrained(model_dir); EsmModel(config); torch.load(model_dir/'pytorch_model.bin'); load the sequential encoder state dict with checkpoint-closure checks; mean pool valid 3Di token states`

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
