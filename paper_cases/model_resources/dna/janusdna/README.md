# janusdna

This directory owns the model identity and executed-environment record for `janusdna`. Constraint-specific input construction, forward execution, coordinate handling and representation export belong to the corresponding implementation under `paper_cases/`. Formal scoring belongs to the selected constraint.

## Environment and checkpoint

| File | Purpose |
|---|---|
| [`ENVIRONMENT.md`](ENVIRONMENT.md) | Human-readable record of the environment used for the formal runs |
| [`ENVIRONMENT.json`](ENVIRONMENT.json) | Machine-readable executed package inventory |
| [`checkpoint_identity.json`](checkpoint_identity.json) | Official source, checkpoint/revision identities, SHA256 evidence, download and loading method |

Recorded model identity:

- Official source: https://github.com/Qihao-Duan/JanusDNA
- Evidence status: `EXECUTED_CHECKPOINT_HASH_RECORDED`
- Upstream revision: `NOT_RECORDED_IN_FORMAL_EXECUTION`
- Executed checkpoint SHA256: `c1dc802dc88c298672f9c1ceb2fda009dab6a2eef2c0551d046fcda3b1b628e2`
- Loading interface: `Import JanusDNA classes from the frozen source tree; construct the model from config.yaml; torch.load(last.ckpt, map_location='cpu')['state_dict']; apply the recorded output-projection compatibility remap; load_state_dict(..., strict=False) with fail-closed missing/unexpected-key checks; eval().cuda()`

Use the exact checkpoint and configuration identity recorded above when reproducing a paper result. A newer checkpoint, tokenizer, representation layer or pooling rule is a different execution and must not be presented as reproduction of a registered result. Historical local paths in the identity record are provenance only; users do not need to reproduce the author's directory layout.

## Constraint-specific forward manuals

Select the constraint you want to run, then follow its execution manual. The manual—not this model index—is authoritative for RefSet inputs, command-line arguments, representation extraction and the file handed to the scorer.

| Constraint | Forward manual |
|---|---|
| DNA01 | [`paper_cases/dna/DNA01/models/janusdna/forward/README.md`](../../../dna/DNA01/models/janusdna/forward/README.md) |
| DNA06 | [`paper_cases/dna/DNA06/models/janusdna/forward/README.md`](../../../dna/DNA06/models/janusdna/forward/README.md) |
| DNA03 | [`paper_cases/dna/DNA03/models/janusdna/forward/README.md`](../../../dna/DNA03/models/janusdna/forward/README.md) |

If a manual identifies an executed historical script without a portable CLI, treat that limitation as authoritative. Do not infer a new command or silently change its model interface.

## Forward-output handoff

For every constraint, the execution manual must identify:

1. the released RefSet files consumed by the forward path;
2. the actual forward/extraction entry point and required parameters;
3. the generated representation or other scorer-input file;
4. the constraint scorer and argument that consume that file.

The exact handoff differs by constraint and may be an embedding archive, position-level representation, distance/feature table or another frozen interface. Do not infer the handoff from a formal result JSON. After the forward step completes, use the scorer documented under the selected constraint's `scoring/` directory.

Formal paper outputs are indexed under `paper_cases/`; they are not forward outputs and are intentionally not duplicated here.
