# ESM-2 650M

This directory owns the model identity and executed-environment record for `esm2-650m`. Constraint-specific input construction, forward execution, coordinate handling and representation export belong to the corresponding implementation case. Formal scoring belongs to the selected constraint.

## Environment and checkpoint

Use these files before running a forward job:

| File | Purpose |
|---|---|
| [`INSTALL.md`](INSTALL.md) | Model installation, checkpoint acquisition, loading and identity checks |
| [`ENVIRONMENT.md`](ENVIRONMENT.md) | Human-readable record of the environment used for the formal runs |
| [`ENVIRONMENT.json`](ENVIRONMENT.json) | Machine-readable executed package inventory |
| [`checkpoint_identity.json`](checkpoint_identity.json) | Official source, executed checkpoint filename, SHA256 and loading interface |

The registered paper runs used:

```text
checkpoint: esm2_t33_650M_UR50D.pt
SHA256: ea9d0522b335a8778dea6535a65301f10208dece28cd5865482b0b1fc446168c
loader: esm.pretrained.load_model_and_alphabet_local
```

The exact upstream source-code revision was not recorded and is explicitly marked as such in `checkpoint_identity.json`. Do not invent one. A newer checkpoint, tokenizer, representation layer or pooling rule is a different execution and must not be presented as reproduction of a registered paper result.

## Constraint-specific forward manuals

Select the constraint you want to run, then follow its execution manual. The manual—not this model index—is authoritative for RefSet inputs, command-line arguments, representation extraction and the file handed to the scorer.

| Constraint | Forward manual |
|---|---|
| PS01 | [`paper_cases/protein_sequence/PS01/models/esm2-650m/forward/IMPLEMENTATION.md`](../../../protein_sequence/PS01/models/esm2-650m/forward/IMPLEMENTATION.md) |
| PS02A | [`paper_cases/protein_sequence/PS02A/models/esm2-650m/README.md`](../../../protein_sequence/PS02A/models/esm2-650m/README.md) |
| PS02B | [`paper_cases/protein_sequence/PS02B/models/esm2-650m/README.md`](../../../protein_sequence/PS02B/models/esm2-650m/README.md) |
| PS04 | [`paper_cases/protein_sequence/PS04/models/esm2-650m/README.md`](../../../protein_sequence/PS04/models/esm2-650m/README.md) |
| PS05 | [`paper_cases/protein_sequence/PS05/models/esm2-650m/README.md`](../../../protein_sequence/PS05/models/esm2-650m/README.md) |
| PS08 | [`paper_cases/protein_sequence/PS08/models/esm2-650m/README.md`](../../../protein_sequence/PS08/models/esm2-650m/README.md) |
| PS09 | [`paper_cases/protein_sequence/PS09/models/esm2-650m/README.md`](../../../protein_sequence/PS09/models/esm2-650m/README.md) |
| PS10 | [`paper_cases/protein_sequence/PS10/models/esm2-650m/README.md`](../../../protein_sequence/PS10/models/esm2-650m/README.md) |

## Forward-output handoff

For every constraint, the execution manual must identify:

1. the released RefSet files consumed by the forward path;
2. the actual forward/extraction entry point and required parameters;
3. the generated representation or other scorer-input file;
4. the constraint scorer and argument that consume that file.

The exact handoff differs by constraint and may be an embedding archive, position-level representation, distance/feature table or another frozen interface. Do not infer the handoff from a formal result JSON. After the forward step completes, use the scorer documented under `constraints/protein_sequence/<constraint_id>/scoring/`.

Formal paper outputs are indexed under `paper_cases/`; they are not forward outputs and are intentionally not duplicated here.
