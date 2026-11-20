# orthrus-4-track

This directory owns the model identity and executed-environment record for `orthrus-4-track`. Constraint-specific input construction, forward execution, coordinate handling and representation export belong to the corresponding implementation under `results/`. Formal scoring belongs to the selected constraint.

## Environment and checkpoint

| File | Purpose |
|---|---|
| [`ENVIRONMENT.md`](ENVIRONMENT.md) | Human-readable record of the environment used for the formal runs |
| [`ENVIRONMENT.json`](ENVIRONMENT.json) | Machine-readable executed package inventory |
| [`checkpoint_identity.json`](checkpoint_identity.json) | Official source, checkpoint/revision identities, SHA256 evidence, download and loading method |

Recorded model identity:

- Official source: https://huggingface.co/antichronology/orthrus-4-track
- Evidence status: `VERIFIED_EXECUTED_CHECKPOINT_IDENTITY`
- Upstream revision: `1326c68932a1f1856e7537df66b131ee3870a0df`
- Executed checkpoint SHA256: `6c2c0f3ec348fead02e9db5b42f230624b9b13ec4cd84678a7322a60cb068b64`
- Loading interface: `AutoModel.from_pretrained(<MODEL_DIR>, trust_remote_code=True, local_files_only=True).to(device).eval(); inputs use model.seq_to_oh; outputs use representation or representation_unpooled`

Use the exact checkpoint and configuration identity recorded above when reproducing a paper result. A newer checkpoint, tokenizer, representation layer or pooling rule is a different execution and must not be presented as reproduction of a registered result. Historical local paths in the identity record are provenance only; users do not need to reproduce the author's directory layout.

## Constraint-specific forward manuals

Select the constraint you want to run, then follow its execution manual. The manual—not this model index—is authoritative for RefSet inputs, command-line arguments, representation extraction and the file handed to the scorer.

| Constraint | Forward manual |
|---|---|
| RNA01-02 | [`paper_cases/rna/RNA01-02/models/orthrus-4-track/forward/IMPLEMENTATION.md`](../../../rna/RNA01-02/models/orthrus-4-track/forward/IMPLEMENTATION.md) |
| RNA04 | [`paper_cases/rna/RNA04/models/orthrus-4-track/forward/README.md`](../../../rna/RNA04/models/orthrus-4-track/forward/README.md) |
| RNA05 | [`paper_cases/rna/RNA05/models/orthrus-4-track/forward/README.md`](../../../rna/RNA05/models/orthrus-4-track/forward/README.md) |
| RNA06 | [`paper_cases/rna/RNA06/models/orthrus-4-track/README.md`](../../../rna/RNA06/models/orthrus-4-track/README.md) |
| RNA09 | [`paper_cases/rna/RNA09/models/orthrus-4-track/forward/README.md`](../../../rna/RNA09/models/orthrus-4-track/forward/README.md) |
| RNA10 | [`paper_cases/rna/RNA10/models/orthrus-4-track/README.md`](../../../rna/RNA10/models/orthrus-4-track/README.md) |

If a manual identifies an executed historical script without a portable CLI, treat that limitation as authoritative. Do not infer a new command or silently change its model interface.

## Forward-output handoff

For every constraint, the execution manual must identify:

1. the released RefSet files consumed by the forward path;
2. the actual forward/extraction entry point and required parameters;
3. the generated representation or other scorer-input file;
4. the constraint scorer and argument that consume that file.

The exact handoff differs by constraint and may be an embedding archive, position-level representation, distance/feature table or another frozen interface. Do not infer the handoff from a formal result JSON. After the forward step completes, use the scorer documented under the selected constraint's `scoring/` directory.

Formal paper outputs are indexed under `results/`; they are not forward outputs and are intentionally not duplicated here.
