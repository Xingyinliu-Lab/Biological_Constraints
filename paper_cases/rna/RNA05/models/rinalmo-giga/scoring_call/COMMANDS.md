# RNA05 scoring call for rinalmo-giga

Run from the repository root after regenerating the Forward output:

```powershell
python constraints/rna/RNA05/scoring/score_rna05.py `
  --vectors work/RNA05/rinalmo-giga/forward/REPRESENTATIONS.npz `
  --model-id rinalmo-giga `
  --output work/RNA05/rinalmo-giga/score
```

The scorer owns validation, frozen cross-fit, ID-level scoring, hierarchical aggregation and bootstrap statistics.
