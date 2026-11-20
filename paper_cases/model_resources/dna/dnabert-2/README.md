# dnabert-2

This directory owns the model identity and executed-environment record for `dnabert-2`. Constraint-specific input construction, forward execution, coordinate handling and representation export belong to the corresponding implementation under `paper_cases/`. Formal scoring belongs to the selected constraint.

## Environment and checkpoint

| File | Purpose |
|---|---|
| [`ENVIRONMENT.md`](ENVIRONMENT.md) | Human-readable record of the environment used for the formal runs |
| [`ENVIRONMENT.json`](ENVIRONMENT.json) | Machine-readable executed package inventory |
| [`checkpoint_identity.json`](checkpoint_identity.json) | Official source, checkpoint/revision identities, SHA256 evidence, download and loading method |

Recorded model identity:

- Official source: https://huggingface.co/zhihan1996/DNABERT-2-117M
- Evidence status: `EXECUTED_CHECKPOINT_HASH_RECORDED`
- Upstream revision: `7bce263b15377fc15361f52cfab88f8b586abda0`
- Executed checkpoint SHA256: `7ff39ec77a484dd01070a41bfd6e95cdd7247bec80fe357ab43a4be33687aeba`
- Loading interface: `AutoTokenizer.from_pretrained(model_path, trust_remote_code=True, local_files_only=True); AutoModel.from_pretrained(model_path, trust_remote_code=True, local_files_only=True).eval()`

Use the exact checkpoint and configuration identity recorded above when reproducing a paper result. A newer checkpoint, tokenizer, representation layer or pooling rule is a different execution and must not be presented as reproduction of a registered result. Historical local paths in the identity record are provenance only; users do not need to reproduce the author's directory layout.

## Constraint-specific forward manuals

Select the constraint you want to run, then follow its execution manual. The manual—not this model index—is authoritative for RefSet inputs, command-line arguments, representation extraction and the file handed to the scorer.

| Constraint | Forward manual |
|---|---|
| DNA01 | [`paper_cases/dna/DNA01/models/dnabert-2/forward/README.md`](../../../dna/DNA01/models/dnabert-2/forward/README.md) |
| DNA06 | [`paper_cases/dna/DNA06/models/dnabert-2/forward/README.md`](../../../dna/DNA06/models/dnabert-2/forward/README.md) |
| DNA03 | [`paper_cases/dna/DNA03/models/dnabert-2/forward/README.md`](../../../dna/DNA03/models/dnabert-2/forward/README.md) |

If a manual identifies an executed historical script without a portable CLI, treat that limitation as authoritative. Do not infer a new command or silently change its model interface.

## Forward-output handoff

For every constraint, the execution manual must identify:

1. the released RefSet files consumed by the forward path;
2. the actual forward/extraction entry point and required parameters;
3. the generated representation or other scorer-input file;
4. the constraint scorer and argument that consume that file.

The exact handoff differs by constraint and may be an embedding archive, position-level representation, distance/feature table or another frozen interface. Do not infer the handoff from a formal result JSON. After the forward step completes, use the scorer documented under the selected constraint's `scoring/` directory.

Formal paper outputs are indexed under `paper_cases/`; they are not forward outputs and are intentionally not duplicated here.
