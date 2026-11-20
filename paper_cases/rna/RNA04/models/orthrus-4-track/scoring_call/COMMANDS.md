# RNA04 / orthrus-4-track public scoring call

After the Forward implementation produces `EMBEDDINGS.npz`, run from the repository root:

```bash
python constraints/rna/RNA04/scoring/score_rna04.py \
  --vectors work/RNA04/orthrus-4-track/EMBEDDINGS.npz \
  --model-id orthrus-4-track \
  --output work/RNA04/orthrus-4-track/score
```

The canonical registered F-only outputs are stored under `output/RNA04_v0.2.2_F_ONLY/`.
