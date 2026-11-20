# PS10 public scoring

[`score_ps10.py`](score_ps10.py) is the sole user-facing primary scorer. It validates the formal record labels, peptide identities and study-aware fold closure plus the model NPZ, performs the frozen 22-fold readout, saves every held-out prediction and fold AUC, then computes formal F and F 95% CI.

```bash
python constraints/protein_sequence/PS10/scoring/score_ps10.py \
  --embeddings work/PS10/my_model/embeddings.npz \
  --model-id my_model \
  --output-dir work/PS10/my_model/score
```

[`ps10_refset.py`](ps10_refset.py) implements the shared RefSet semantic checks. [`ps10_scoring_core.py`](ps10_scoring_core.py) is the frozen numeric core shared with the TANGO comparator; neither is a second public scorer. [`PS10_NUMERIC_METHOD.json`](PS10_NUMERIC_METHOD.json) records the frozen method.

See [`../README.md`](../README.md) for the complete Forward contract, output columns and interpretation.
