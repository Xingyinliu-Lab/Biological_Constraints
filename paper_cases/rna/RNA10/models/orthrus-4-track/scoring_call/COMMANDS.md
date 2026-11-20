# RNA10 F-only scoring

```bash
python constraints/rna/RNA10/scoring/build_measurement_from_vectors.py --truth-relations constraints/rna/RNA10/refset/truth_relations.jsonl --objects-wt constraints/rna/RNA10/refset/objects_wt.jsonl --objects-mutant constraints/rna/RNA10/refset/objects_mutant.jsonl --vectors work/RNA10/orthrus-4-track/forward/REPRESENTATIONS.npz --out work/RNA10/orthrus-4-track/score --model orthrus-4-track
python constraints/rna/RNA10/scoring/score_rna10.py --measurement-ledger work/RNA10/orthrus-4-track/score/MEASUREMENT_LEDGER.json --model-id orthrus-4-track --out work/RNA10/orthrus-4-track/score/FORMAL_RESULT.json
```
