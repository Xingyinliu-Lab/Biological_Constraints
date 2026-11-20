# DNA03 / caduceus-ps commands

Use the checkpoint and environment recorded in `paper_cases/model_resources/dna/caduceus-ps/`. The adapter reads the released FASTA and membership table and writes primary and repeat NPZ files directly in the DNA03 public contract. The output directory must not already exist.

```bash
mkdir -p work/DNA03/caduceus-ps
DNA03_MODEL_PATH=/path/to/caduceus-ps-snapshot \
python paper_cases/dna/DNA03/models/caduceus-ps/forward/code/encode_dna03_transformer_models_v020.py \
  --model caduceus_ps \
  --membership constraints/dna/DNA03/refset/CORE_VIEW_MEMBERSHIP.tsv \
  --fasta constraints/dna/DNA03/refset/DNA03_FORWARD_VIEWS.fasta \
  --output work/DNA03/caduceus-ps/forward
```


Then call the public scorer:

```bash
python constraints/dna/DNA03/scoring/score_dna03.py \
  --embeddings work/DNA03/caduceus-ps/forward/EMBEDDINGS.npz \
  --repeat-embeddings work/DNA03/caduceus-ps/forward/REPEAT_EMBEDDINGS.npz \
  --model-id caduceus-ps \
  --output work/DNA03/caduceus-ps/score
```
