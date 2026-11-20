# RNA05 scoring call for gb-rna

Run from the repository root after regenerating the Forward output:

```powershell
python constraints/rna/RNA05/scoring/score_rna05.py `
  --vectors work/RNA05/gb-rna/forward/REPRESENTATIONS.npz `
  --model-id gb-rna `
  --output work/RNA05/gb-rna/score
```

The scorer owns validation, frozen cross-fit, ID-level scoring, hierarchical aggregation and bootstrap statistics.
