# carp-640m

This directory owns the model identity and executed-environment record for `carp-640m`. Constraint-specific input construction, forward execution, coordinate handling and representation export belong to the corresponding implementation case. Formal scoring belongs to the selected constraint.

## Environment and checkpoint

| File | Purpose |
|---|---|
| [`ENVIRONMENT.md`](ENVIRONMENT.md) | Human-readable record of the environment used for the formal runs |
| [`ENVIRONMENT.json`](ENVIRONMENT.json) | Machine-readable executed package inventory |
| [`checkpoint_identity.json`](checkpoint_identity.json) | Official source, checkpoint/revision identities, SHA256 evidence, download and loading method |

Recorded model identity:

- Official source: https://github.com/microsoft/protein-sequence-models
- Evidence status: `EXECUTED_CHECKPOINT_HASH_RECORDED`
- Upstream revision: `NOT_RECORDED_IN_FORMAL_EXECUTION`
- Executed checkpoint SHA256: `ae1d7e0cbe713c731cfde7ace213f2a73372c31d8df543089e6b6d371228c4f7`
- Loading interface: `from sequence_models.pretrained import load_model_and_alphabet; model, collater = load_model_and_alphabet('carp_640M'); model.eval().cuda(); final ByteNet layer residue states are mean pooled`

Use the exact checkpoint and configuration identity recorded above when reproducing a paper result. A newer checkpoint, tokenizer, representation layer or pooling rule is a different execution and must not be presented as reproduction of a registered result. Historical local paths in the identity record are provenance only; users do not need to reproduce the author's directory layout.

## Constraint-specific forward manuals

Select the constraint you want to run, then follow its execution manual. The manual—not this model index—is authoritative for RefSet inputs, command-line arguments, representation extraction and the file handed to the scorer.

| Constraint | Forward manual |
|---|---|
| PS01 | [`paper_cases/protein_sequence/PS01/models/carp-640m/forward/IMPLEMENTATION.md`](../../../protein_sequence/PS01/models/carp-640m/forward/IMPLEMENTATION.md) |
| PS02A | [`paper_cases/protein_sequence/PS02A/models/carp-640m/README.md`](../../../protein_sequence/PS02A/models/carp-640m/README.md) |
| PS02B | [`paper_cases/protein_sequence/PS02B/models/carp-640m/README.md`](../../../protein_sequence/PS02B/models/carp-640m/README.md) |
| PS04 | [`paper_cases/protein_sequence/PS04/models/carp-640m/README.md`](../../../protein_sequence/PS04/models/carp-640m/README.md) |
| PS05 | [`paper_cases/protein_sequence/PS05/models/carp-640m/README.md`](../../../protein_sequence/PS05/models/carp-640m/README.md) |
| PS08 | [`paper_cases/protein_sequence/PS08/models/carp-640m/README.md`](../../../protein_sequence/PS08/models/carp-640m/README.md) |
| PS09 | [`paper_cases/protein_sequence/PS09/models/carp-640m/README.md`](../../../protein_sequence/PS09/models/carp-640m/README.md) |
| PS10 | [`paper_cases/protein_sequence/PS10/models/carp-640m/README.md`](../../../protein_sequence/PS10/models/carp-640m/README.md) |

If a manual identifies an executed historical script without a portable CLI, treat that limitation as authoritative. Do not infer a new command or silently change its model interface.

## Forward-output handoff

For every constraint, the execution manual must identify:

1. the released RefSet files consumed by the forward path;
2. the actual forward/extraction entry point and required parameters;
3. the generated representation or other scorer-input file;
4. the constraint scorer and argument that consume that file.

The exact handoff differs by constraint and may be an embedding archive, position-level representation, distance/feature table or another frozen interface. Do not infer the handoff from a formal result JSON. After the forward step completes, use the scorer documented under the selected constraint's `scoring/` directory.

Formal paper outputs are indexed under `paper_cases/`; they are not forward outputs and are intentionally not duplicated here.
