# RNA10 F-only scoring

```bash
python constraints/rna/RNA10/scoring/build_measurement_from_vectors.py --truth-relations constraints/rna/RNA10/refset/truth_relations.jsonl --objects-wt constraints/rna/RNA10/refset/objects_wt.jsonl --objects-mutant constraints/rna/RNA10/refset/objects_mutant.jsonl --vectors work/RNA10/gb-rna/forward/REPRESENTATIONS.npz --out work/RNA10/gb-rna/score --model gb-rna
python constraints/rna/RNA10/scoring/score_rna10.py --measurement-ledger work/RNA10/gb-rna/score/MEASUREMENT_LEDGER.json --model-id gb-rna --out work/RNA10/gb-rna/score/FORMAL_RESULT.json
```
