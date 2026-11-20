# RNA05 scoring call for ernie-rna

Run from the repository root after regenerating the Forward output:

```powershell
python constraints/rna/RNA05/scoring/score_rna05.py `
  --vectors work/RNA05/ernie-rna/forward/REPRESENTATIONS.npz `
  --model-id ernie-rna `
  --output work/RNA05/ernie-rna/score
```

The scorer owns validation, frozen cross-fit, ID-level scoring, hierarchical aggregation and bootstrap statistics.
