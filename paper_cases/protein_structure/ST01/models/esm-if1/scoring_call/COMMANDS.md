# Public ST01 scoring call

After generating a compatible representation file, run the single public scorer:

```bash
python constraints/protein_structure/ST01/scoring/score_st01_direct.py \
  --embeddings work/ST01/esm-if1/embeddings.npz \
  --output-dir work/ST01/esm-if1/score
```

The scorer writes `INDIVIDUAL_SCORES.tsv`, endpoint/anchor/cluster ledgers, `SCORE_SUMMARY.json`, the bootstrap sample and `OUTPUT_MANIFEST.json`. The reconstructed registered output is stored under [`../output/`](../output/).
