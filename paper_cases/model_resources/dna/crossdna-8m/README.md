# crossdna-8m

This directory owns the model identity and executed-environment record for `crossdna-8m`. Constraint-specific input construction, forward execution, coordinate handling and representation export belong to the corresponding implementation under `paper_cases/`. Formal scoring belongs to the selected constraint.

## Environment and checkpoint

| File | Purpose |
|---|---|
| [`ENVIRONMENT.md`](ENVIRONMENT.md) | Human-readable record of the environment used for the formal runs |
| [`ENVIRONMENT.json`](ENVIRONMENT.json) | Machine-readable executed package inventory |
| [`checkpoint_identity.json`](checkpoint_identity.json) | Official source, checkpoint/revision identities, SHA256 evidence, download and loading method |

Recorded model identity:

- Official source: https://github.com/LuoGroup2023/CrossDNA
- Evidence status: `EXECUTED_MULTI_FILE_CHECKPOINT_IDENTITY_RECORDED`
- Upstream revision: `NOT_RECORDED_IN_FORMAL_EXECUTION`
- Executed checkpoint SHA256: `886fc527a7103aa27d4f7f41626a9a283f400983e91d5e09247243d79b1dac63`, `9fb0669559f4aff27b88a26fcd52f048b4d3ee36a7c3ea35ef07e049c355d250`
- Loading interface: `AutoTokenizer.from_pretrained(model_dir, trust_remote_code=True, local_files_only=True); AutoModelForMaskedLM.from_pretrained(model_dir, trust_remote_code=True, local_files_only=True).eval().half().cuda()`

Use the exact checkpoint and configuration identity recorded above when reproducing a paper result. A newer checkpoint, tokenizer, representation layer or pooling rule is a different execution and must not be presented as reproduction of a registered result. Historical local paths in the identity record are provenance only; users do not need to reproduce the author's directory layout.

## Constraint-specific forward manuals

Select the constraint you want to run, then follow its execution manual. The manual—not this model index—is authoritative for RefSet inputs, command-line arguments, representation extraction and the file handed to the scorer.

| Constraint | Forward manual |
|---|---|
| DNA01 | [`paper_cases/dna/DNA01/models/crossdna-8m/forward/README.md`](../../../dna/DNA01/models/crossdna-8m/forward/README.md) |
| DNA06 | [`paper_cases/dna/DNA06/models/crossdna-8m/forward/README.md`](../../../dna/DNA06/models/crossdna-8m/forward/README.md) |
| DNA03 | [`paper_cases/dna/DNA03/models/crossdna-8m/forward/README.md`](../../../dna/DNA03/models/crossdna-8m/forward/README.md) |

For DNA03, the case-specific forward implementation uses its frozen hand-coded nucleotide encoding and backbone readout; it does not use an `AutoTokenizer`.

If a manual identifies an executed historical script without a portable CLI, treat that limitation as authoritative. Do not infer a new command or silently change its model interface.

## Forward-output handoff

For every constraint, the execution manual must identify:

1. the released RefSet files consumed by the forward path;
2. the actual forward/extraction entry point and required parameters;
3. the generated representation or other scorer-input file;
4. the constraint scorer and argument that consume that file.

The exact handoff differs by constraint and may be an embedding archive, position-level representation, distance/feature table or another frozen interface. Do not infer the handoff from a formal result JSON. After the forward step completes, use the scorer documented under the selected constraint's `scoring/` directory.

Formal paper outputs are indexed under `paper_cases/`; they are not forward outputs and are intentionally not duplicated here.
