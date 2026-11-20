# RNA06 scoring call for rna-fm

Run from the repository root after regenerating the Forward output:

```powershell
python constraints/rna/RNA06/scoring/score_rna06.py `
  --vectors work/RNA06/rna-fm/forward/REPRESENTATIONS.npz `
  --model rna-fm `
  --output work/RNA06/rna-fm/score
```

The scorer loads the full formal RefSet automatically and owns feature construction, frozen cross-fit, ID-level scoring, hierarchical aggregation and bootstrap statistics.
