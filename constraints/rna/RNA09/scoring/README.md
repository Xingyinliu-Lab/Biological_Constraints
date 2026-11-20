# RNA09 scoring

Create the standalone environment described in [`ENVIRONMENT.md`](ENVIRONMENT.md), then run from the repository root:

```bash
python constraints/rna/RNA09/scoring/score_official.py \
  --embeddings <REPRESENTATIONS.npz> \
  --output-dir work/RNA09/<model>/score
```

The optional RefSet path arguments are input-location overrides only; they do not define another scoring mode. The representation NPZ must contain exactly the released endpoint sequence universe and finite real vectors.

The scorer writes `LEAF_LEDGER.tsv` (`F_credit`), anchor, source-parent, dependency-unit and group F ledgers, `F_BOOTSTRAP_MEANS.npy`, `FORMAL_RESULT.json`, and `OUTPUT_MANIFEST.json`. `F_95_CI` is a 20,000-draw percentile dependency-unit bootstrap interval.

All paper models, third-party models and comparator measurement producers use this one F-only credit, hierarchy and bootstrap implementation.
