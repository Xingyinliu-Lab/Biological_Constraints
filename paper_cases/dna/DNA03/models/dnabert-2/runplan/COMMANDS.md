# DNA03 / dnabert-2 commands

Use the checkpoint and environment recorded in `paper_cases/model_resources/dna/dnabert-2/`. The adapter reads the released FASTA and membership table and writes primary and repeat NPZ files directly in the DNA03 public contract. The output directory must not already exist.

```bash
mkdir -p work/DNA03/dnabert-2
DNA03_MODEL_PATH=/path/to/DNABERT-2-117M-snapshot \
python paper_cases/dna/DNA03/models/dnabert-2/forward/code/encode_dna03_transformer_models_v020.py \
  --model dnabert2 \
  --membership constraints/dna/DNA03/refset/CORE_VIEW_MEMBERSHIP.tsv \
  --fasta constraints/dna/DNA03/refset/DNA03_FORWARD_VIEWS.fasta \
  --output work/DNA03/dnabert-2/forward
```


Then call the public scorer:

```bash
python constraints/dna/DNA03/scoring/score_dna03.py \
  --embeddings work/DNA03/dnabert-2/forward/EMBEDDINGS.npz \
  --repeat-embeddings work/DNA03/dnabert-2/forward/REPEAT_EMBEDDINGS.npz \
  --model-id dnabert-2 \
  --output work/DNA03/dnabert-2/score
```
