# caduceus-ps

This directory owns the model identity and executed-environment record for `caduceus-ps`. Constraint-specific input construction, forward execution, coordinate handling and representation export belong to the corresponding implementation under `paper_cases/`. Formal scoring belongs to the selected constraint.

## Environment and checkpoint

| File | Purpose |
|---|---|
| [`ENVIRONMENT.md`](ENVIRONMENT.md) | Human-readable record of the environment used for the formal runs |
| [`ENVIRONMENT.json`](ENVIRONMENT.json) | Machine-readable executed package inventory |
| [`checkpoint_identity.json`](checkpoint_identity.json) | Official source, checkpoint/revision identities, SHA256 evidence, download and loading method |

Recorded model identity:

- Official source: https://huggingface.co/kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16
- Evidence status: `EXECUTED_CHECKPOINT_HASH_RECORDED`
- Upstream revision: `d89eeb853136ea64da7feb3d0c8e909771b17ae6`
- Executed checkpoint SHA256: `a3e6976fe90460ff5d90d371b457d0263dc65e156ea5651b4a452e98228daef2`
- Loading interface: `AutoTokenizer.from_pretrained(model_path, trust_remote_code=True, local_files_only=True); AutoModel.from_pretrained(model_path, trust_remote_code=True, local_files_only=True).eval()`

Use the exact checkpoint and configuration identity recorded above when reproducing a paper result. A newer checkpoint, tokenizer, representation layer or pooling rule is a different execution and must not be presented as reproduction of a registered result. Historical local paths in the identity record are provenance only; users do not need to reproduce the author's directory layout.

## Constraint-specific forward manuals

Select the constraint you want to run, then follow its execution manual. The manual—not this model index—is authoritative for RefSet inputs, command-line arguments, representation extraction and the file handed to the scorer.

| Constraint | Forward manual |
|---|---|
| DNA01 | [`paper_cases/dna/DNA01/models/caduceus-ps/forward/README.md`](../../../dna/DNA01/models/caduceus-ps/forward/README.md) |
| DNA06 | [`paper_cases/dna/DNA06/models/caduceus-ps/forward/README.md`](../../../dna/DNA06/models/caduceus-ps/forward/README.md) |
| DNA03 | [`paper_cases/dna/DNA03/models/caduceus-ps/forward/README.md`](../../../dna/DNA03/models/caduceus-ps/forward/README.md) |

For DNA03, the case-specific forward implementation loads `AutoModelForMaskedLM`; that documented implementation supersedes the generic loading-interface summary above for this constraint.

If a manual identifies an executed historical script without a portable CLI, treat that limitation as authoritative. Do not infer a new command or silently change its model interface.

## Forward-output handoff

For every constraint, the execution manual must identify:

1. the released RefSet files consumed by the forward path;
2. the actual forward/extraction entry point and required parameters;
3. the generated representation or other scorer-input file;
4. the constraint scorer and argument that consume that file.

The exact handoff differs by constraint and may be an embedding archive, position-level representation, distance/feature table or another frozen interface. Do not infer the handoff from a formal result JSON. After the forward step completes, use the scorer documented under the selected constraint's `scoring/` directory.

Formal paper outputs are indexed under `paper_cases/`; they are not forward outputs and are intentionally not duplicated here.
