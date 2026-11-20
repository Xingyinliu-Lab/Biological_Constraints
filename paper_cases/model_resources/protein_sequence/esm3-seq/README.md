# esm3-seq

This directory owns the model identity and executed-environment record for `esm3-seq`. Constraint-specific input construction, forward execution, coordinate handling and representation export belong to the corresponding implementation case. Formal scoring belongs to the selected constraint.

## Environment and checkpoint

| File | Purpose |
|---|---|
| [`ENVIRONMENT.md`](ENVIRONMENT.md) | Human-readable record of the environment used for the formal runs |
| [`ENVIRONMENT.json`](ENVIRONMENT.json) | Machine-readable executed package inventory |
| [`checkpoint_identity.json`](checkpoint_identity.json) | Official source, checkpoint/revision identities, SHA256 evidence, download and loading method |

Recorded model identity:

- Official source: https://huggingface.co/EvolutionaryScale/esm3-sm-open-v1
- Evidence status: `EXECUTED_CHECKPOINT_HASH_RECORDED`
- Upstream revision: `NOT_RECORDED_IN_FORMAL_EXECUTION`
- Executed checkpoint SHA256: `752362cd8ed1e7826ac366197f3629c47b9a3996457defff30a922b5fde01b6f`
- Loading interface: `ESM3.from_pretrained('esm3_sm_open_v1', device=device).eval(); sequence track only; other tracks omitted/masked; residue states mean pooled`

Use the exact checkpoint and configuration identity recorded above when reproducing a paper result. A newer checkpoint, tokenizer, representation layer or pooling rule is a different execution and must not be presented as reproduction of a registered result. Historical local paths in the identity record are provenance only; users do not need to reproduce the author's directory layout.

## Constraint-specific forward manuals

Select the constraint you want to run, then follow its execution manual. The manual—not this model index—is authoritative for RefSet inputs, command-line arguments, representation extraction and the file handed to the scorer.

| Constraint | Forward manual |
|---|---|
| PS01 | [`paper_cases/protein_sequence/PS01/models/esm3-seq/forward/IMPLEMENTATION.md`](../../../protein_sequence/PS01/models/esm3-seq/forward/IMPLEMENTATION.md) |
| PS02A | [`paper_cases/protein_sequence/PS02A/models/esm3-seq/README.md`](../../../protein_sequence/PS02A/models/esm3-seq/README.md) |
| PS02B | [`paper_cases/protein_sequence/PS02B/models/esm3-seq/README.md`](../../../protein_sequence/PS02B/models/esm3-seq/README.md) |
| PS04 | [`paper_cases/protein_sequence/PS04/models/esm3-seq/README.md`](../../../protein_sequence/PS04/models/esm3-seq/README.md) |
| PS05 | [`paper_cases/protein_sequence/PS05/models/esm3-seq/README.md`](../../../protein_sequence/PS05/models/esm3-seq/README.md) |
| PS08 | [`paper_cases/protein_sequence/PS08/models/esm3-seq/README.md`](../../../protein_sequence/PS08/models/esm3-seq/README.md) |
| PS09 | [`paper_cases/protein_sequence/PS09/models/esm3-seq/README.md`](../../../protein_sequence/PS09/models/esm3-seq/README.md) |
| PS10 | [`paper_cases/protein_sequence/PS10/models/esm3-seq/README.md`](../../../protein_sequence/PS10/models/esm3-seq/README.md) |

If a manual identifies an executed historical script without a portable CLI, treat that limitation as authoritative. Do not infer a new command or silently change its model interface.

## Forward-output handoff

For every constraint, the execution manual must identify:

1. the released RefSet files consumed by the forward path;
2. the actual forward/extraction entry point and required parameters;
3. the generated representation or other scorer-input file;
4. the constraint scorer and argument that consume that file.

The exact handoff differs by constraint and may be an embedding archive, position-level representation, distance/feature table or another frozen interface. Do not infer the handoff from a formal result JSON. After the forward step completes, use the scorer documented under the selected constraint's `scoring/` directory.

Formal paper outputs are indexed under `paper_cases/`; they are not forward outputs and are intentionally not duplicated here.
