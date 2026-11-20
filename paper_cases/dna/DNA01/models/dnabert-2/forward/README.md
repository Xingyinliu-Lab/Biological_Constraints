
# Forward implementation

`scripts/encode_dna01_transformer_models.py` is the model-by-DNA01 Forward adapter retained from the final FULL evidence. It reads the released 2,920-object DNA01 RefSet and writes `EMBEDDINGS.npz`, `REPEAT_EMBEDDINGS.npz`, Forward validation and identity metadata.

The model environment and checkpoint identity are stored once at `paper_cases/model_resources/dna/dnabert-2/`. This directory performs model inference only; it does not calculate DNA01 scores.

To reproduce the Forward execution, follow [`../runplan/COMMANDS.md`](../runplan/COMMANDS.md). The frozen executed parameters remain in [`../runplan/RUNPLAN.json`](../runplan/RUNPLAN.json).
