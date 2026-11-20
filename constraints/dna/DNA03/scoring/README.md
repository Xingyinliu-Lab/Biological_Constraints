# DNA03 scoring

`score_dna03.py` is the only public DNA03 scoring entry. It accepts the primary and repeat core-local representation NPZ files documented in [`../README.md`](../README.md) and performs the complete frozen measurement after model Forward.

```bash
python constraints/dna/DNA03/scoring/score_dna03.py \
  --embeddings <PRIMARY_NPZ> \
  --repeat-embeddings <REPEAT_NPZ> \
  --model-id <MODEL_ID> \
  --output <EMPTY_OUTPUT_DIRECTORY>
```

The RefSet paths default to the assets released beside this program. `--queries`, `--negatives` and `--controls` remain available only for an explicitly compatible RefSet; the actual input identities are recorded in `FORMAL_RESULT.json`.

Internal modules:

- `dna03_tau_core.py`: validates the representation contract and computes tau from the frozen repeat controls.
- `dna03_numeric_core.py`: frozen distance, tie-credit, hierarchical aggregation and bootstrap functions.
- `score_dna03.py`: orchestration, population validation and formal output generation.

All fixed formal parameters are prescribed by the release. See [`ENVIRONMENT.md`](ENVIRONMENT.md) for the scoring-only environment and the constraint [README](../README.md) for input and output fields.

The formal output retains relation-level raw distances, margin, tie and `F_credit`; query, cluster-transition and transition F ledgers; and the raw primary/repeat control measurements used to select tau. The final transition ledger macro-average is checked against reported F.
