# RNA10 F-only scoring

```bash
python constraints/rna/RNA10/scoring/build_measurement_from_vectors.py --truth-relations constraints/rna/RNA10/refset/truth_relations.jsonl --objects-wt constraints/rna/RNA10/refset/objects_wt.jsonl --objects-mutant constraints/rna/RNA10/refset/objects_mutant.jsonl --vectors work/RNA10/rinalmo-giga/forward/REPRESENTATIONS.npz --out work/RNA10/rinalmo-giga/score --model rinalmo-giga
python constraints/rna/RNA10/scoring/score_rna10.py --measurement-ledger work/RNA10/rinalmo-giga/score/MEASUREMENT_LEDGER.json --model-id rinalmo-giga --out work/RNA10/rinalmo-giga/score/FORMAL_RESULT.json
```
