# rinalmo-giga

This directory owns the model identity and executed-environment record for `rinalmo-giga`. Constraint-specific input construction, forward execution, coordinate handling and representation export belong to the corresponding implementation under `results/`. Formal scoring belongs to the selected constraint.

## Environment and checkpoint

| File | Purpose |
|---|---|
| [`ENVIRONMENT.md`](ENVIRONMENT.md) | Human-readable record of the environment used for the formal runs |
| [`ENVIRONMENT.json`](ENVIRONMENT.json) | Machine-readable executed package inventory |
| [`checkpoint_identity.json`](checkpoint_identity.json) | Official source, checkpoint/revision identities, SHA256 evidence, download and loading method |

Recorded model identity:

- Official source: https://github.com/lbcb-sci/RiNALMo
- Evidence status: `EXECUTED_CHECKPOINT_HASH_RECORDED`
- Upstream revision: `2c2c5c14a5ae609d8c560a5d9ca32e51e0288955`
- Executed checkpoint SHA256: `cd93c3f21eb3e767373c9491192686b5846247bd1110693e453c1dd0f321c0db`
- Loading interface: `from rinalmo.pretrained import get_pretrained_model; model, alphabet = get_pretrained_model('giga-v1'); model.eval().to(device); exclude boundary tokens and reconstruct nucleotide coordinates before pooling`

Use the exact checkpoint and configuration identity recorded above when reproducing a paper result. A newer checkpoint, tokenizer, representation layer or pooling rule is a different execution and must not be presented as reproduction of a registered result. Historical local paths in the identity record are provenance only; users do not need to reproduce the author's directory layout.

## Constraint-specific forward manuals

Select the constraint you want to run, then follow its execution manual. The manual—not this model index—is authoritative for RefSet inputs, command-line arguments, representation extraction and the file handed to the scorer.

| Constraint | Forward manual |
|---|---|
| RNA01-02 | [`paper_cases/rna/RNA01-02/models/rinalmo-giga/forward/IMPLEMENTATION.md`](../../../rna/RNA01-02/models/rinalmo-giga/forward/IMPLEMENTATION.md) |
| RNA04 | [`paper_cases/rna/RNA04/models/rinalmo-giga/forward/README.md`](../../../rna/RNA04/models/rinalmo-giga/forward/README.md) |
| RNA05 | [`paper_cases/rna/RNA05/models/rinalmo-giga/forward/README.md`](../../../rna/RNA05/models/rinalmo-giga/forward/README.md) |
| RNA06 | [`paper_cases/rna/RNA06/models/rinalmo-giga/README.md`](../../../rna/RNA06/models/rinalmo-giga/README.md) |
| RNA09 | [`paper_cases/rna/RNA09/models/rinalmo-giga/forward/README.md`](../../../rna/RNA09/models/rinalmo-giga/forward/README.md) |
| RNA10 | [`paper_cases/rna/RNA10/models/rinalmo-giga/README.md`](../../../rna/RNA10/models/rinalmo-giga/README.md) |

If a manual identifies an executed historical script without a portable CLI, treat that limitation as authoritative. Do not infer a new command or silently change its model interface.

## Forward-output handoff

For every constraint, the execution manual must identify:

1. the released RefSet files consumed by the forward path;
2. the actual forward/extraction entry point and required parameters;
3. the generated representation or other scorer-input file;
4. the constraint scorer and argument that consume that file.

The exact handoff differs by constraint and may be an embedding archive, position-level representation, distance/feature table or another frozen interface. Do not infer the handoff from a formal result JSON. After the forward step completes, use the scorer documented under the selected constraint's `scoring/` directory.

Formal paper outputs are indexed under `results/`; they are not forward outputs and are intentionally not duplicated here.
