# RNA04 / rinalmo-giga public scoring call

After the Forward implementation produces `EMBEDDINGS.npz`, run from the repository root:

```bash
python constraints/rna/RNA04/scoring/score_rna04.py \
  --vectors work/RNA04/rinalmo-giga/EMBEDDINGS.npz \
  --model-id rinalmo-giga \
  --output work/RNA04/rinalmo-giga/score
```

The canonical registered F-only outputs are stored under `output/RNA04_v0.2.2_F_ONLY/`.
