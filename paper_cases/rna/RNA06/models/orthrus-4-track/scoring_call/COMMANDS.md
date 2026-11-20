# RNA06 scoring call for orthrus-4-track

Run from the repository root after regenerating the Forward output:

```powershell
python constraints/rna/RNA06/scoring/score_rna06.py `
  --vectors work/RNA06/orthrus-4-track/forward/REPRESENTATIONS.npz `
  --model orthrus-4-track `
  --output work/RNA06/orthrus-4-track/score
```

The scorer loads the full formal RefSet automatically and owns feature construction, frozen cross-fit, ID-level scoring, hierarchical aggregation and bootstrap statistics.
