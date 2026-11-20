# RNA05 scoring call for rna-fm

Run from the repository root after regenerating the Forward output:

```powershell
python constraints/rna/RNA05/scoring/score_rna05.py `
  --vectors work/RNA05/rna-fm/forward/REPRESENTATIONS.npz `
  --model-id rna-fm `
  --output work/RNA05/rna-fm/score
```

The scorer owns validation, frozen cross-fit, ID-level scoring, hierarchical aggregation and bootstrap statistics.
