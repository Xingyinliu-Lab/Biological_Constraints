# RNA04 / rna-fm public scoring call

After the Forward implementation produces `EMBEDDINGS.npz`, run from the repository root:

```bash
python constraints/rna/RNA04/scoring/score_rna04.py \
  --vectors work/RNA04/rna-fm/EMBEDDINGS.npz \
  --model-id rna-fm \
  --output work/RNA04/rna-fm/score
```

The canonical registered F-only outputs are stored under `output/RNA04_v0.2.2_F_ONLY/`.
