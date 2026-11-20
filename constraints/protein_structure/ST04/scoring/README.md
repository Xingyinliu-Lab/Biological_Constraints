# ST04 scoring

## Environment

Create the standalone scoring environment described in [`ENVIRONMENT.md`](ENVIRONMENT.md). It is intentionally separate from model CUDA environments.

## Programs

### `score_st04.py`

Run the frozen ST04 method on a compatible residue representation.

```bash
python constraints/protein_structure/ST04/scoring/score_st04.py \
  --relations <RELATIONS> \
  --refset-release <REFSET_RELEASE> \
  --numeric-method <NUMERIC_METHOD> \
  --split-release <SPLIT_RELEASE> \
  --embeddings <EMBEDDINGS> \
  --model-id <MODEL_ID> \
  --out-dir <OUT_DIR>
```

All compatible models are scored against the complete released RefSet.

### `st04_numeric_core_v0_1_1.py`

Frozen mathematical implementation imported by `score_st04.py`. It contains
the nested cross-fit, lambda selection, aggregation and bootstrap calculations;
it is not a second scorer and exposes no CLI.

### `st04_io_qc.py`

Internal I/O, schema validation, identity and output-writing support imported by
`score_st04.py`. It does not calculate an alternative score and exposes no CLI.

## Inputs

Use the exact RefSet paths and representation schema documented in [`../README.md`](../README.md). Inputs must cover the complete formal ID universe, contain finite numeric values and preserve frozen folds/dependency units.

## Outputs

Write outputs to `work/ST04/<model>/score/`. `EVALUATION_RESULT.json` contains
the aggregate F/CI and population identity. `RELATION_SCORE_LEDGER.tsv.gz`
is the required relation-level audit output: it records the frozen relation
identity, positive/negative logits, delta, tie classification, credit, and
query contribution. `OOF_QUERY_SCORES.jsonl` is the
formal ID-level output and records held-out raw positive/negative scores,
relation credits and query credits. `HIERARCHY_SCORES.json`, coverage, release
QC, selected lambdas and fit-parameter identities make the aggregation and
cross-fit traceable. These evidence files are public scoring outputs, not
disposable runtime caches.

## One public scoring pipeline
