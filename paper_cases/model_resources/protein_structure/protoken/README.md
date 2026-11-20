# protoken

This directory owns the shared model identity and executed-environment record for `protoken`. Constraint-specific input construction, Forward execution, coordinate handling and representation export belong to the corresponding implementation under `paper_cases/`. Formal scoring belongs to the selected Constraint.

## Environment and checkpoint

| File | Purpose |
|---|---|
| [`ENVIRONMENT.md`](ENVIRONMENT.md) | Human-readable record of the environment used for the formal runs |
| [`ENVIRONMENT.json`](ENVIRONMENT.json) | Machine-readable executed package inventory |
| [`checkpoint_identity.json`](checkpoint_identity.json) | Official source, checkpoint/revision identities, SHA256 evidence, download and loading method |

Recorded model identity:

- Official source: https://github.com/issacAzazel/ProToken
- Evidence status: `EXECUTED_CHECKPOINT_HASH_RECORDED`
- Upstream revision: `NOT_RECORDED_IN_FORMAL_EXECUTION`
- Executed checkpoint SHA256: `4fed4b27d615fbea77373c2f02fa22dd0cc9351d24802d260f0e2a9781a6ffef`
- Loading interface: `Run the frozen ProToken stage-1 inference with encoder.yaml, decoder.yaml, vq.yaml and --load_ckpt_path protoken_params_100000.pkl; map token IDs through the frozen VQ codebook and mean pool codebook vectors`

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
