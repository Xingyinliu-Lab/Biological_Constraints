# RNA05 scoring call for orthrus-4-track

Run from the repository root after regenerating the Forward output:

```powershell
python constraints/rna/RNA05/scoring/score_rna05.py `
  --vectors work/RNA05/orthrus-4-track/forward/REPRESENTATIONS.npz `
  --model-id orthrus-4-track `
  --output work/RNA05/orthrus-4-track/score
```

The scorer owns validation, frozen cross-fit, ID-level scoring, hierarchical aggregation and bootstrap statistics.
