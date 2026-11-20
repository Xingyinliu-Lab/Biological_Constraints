# RNA06 scoring call for rinalmo-giga

Run from the repository root after regenerating the Forward output:

```powershell
python constraints/rna/RNA06/scoring/score_rna06.py `
  --vectors work/RNA06/rinalmo-giga/forward/REPRESENTATIONS.npz `
  --model rinalmo-giga `
  --output work/RNA06/rinalmo-giga/score
```

The scorer loads the full formal RefSet automatically and owns feature construction, frozen cross-fit, ID-level scoring, hierarchical aggregation and bootstrap statistics.
