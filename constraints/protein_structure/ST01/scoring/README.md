# ST01 public scoring package

`score_st01_direct.py` is the single public ST01 scorer for paper models, third-party models and compatible reference methods. It owns all ST01-specific cosine measurement, tie-aware contrast credit, hierarchy aggregation, bootstrap and public output generation.

## Model representation input

```bash
python constraints/protein_structure/ST01/scoring/score_st01_direct.py \
  --embeddings work/ST01/my_model/embeddings.npz \
  --output-dir work/ST01/my_model/score
```

The NPZ contract and output fields are defined in the [Constraint README](../README.md). The complete public mode writes ID-level, endpoint, anchor and cluster tables before the aggregate summary, and records the actual Forward input mode and SHA256 in `SCORE_SUMMARY.json`.

## Pair-score input

```bash
python constraints/protein_structure/ST01/scoring/score_st01_direct.py \
  --pair-scores compatible_pair_scores.tsv \
  --output-dir work/ST01/reference/score
```

This mode uses the same frozen aggregation and uncertainty procedure and writes the complete public F ledger chain.

The formal values for tie tolerance, bootstrap generator, seed and draw count come from the bundled Numeric Method and are not user-selectable scoring decisions.
