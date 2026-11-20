# DNA03 / nt-v2 commands

Use the checkpoint and environment recorded in `paper_cases/model_resources/dna/nt-v2/`. The adapter reads the released FASTA and membership table and writes primary and repeat NPZ files directly in the DNA03 public contract. The output directory must not already exist.

```bash
mkdir -p work/DNA03/nt-v2
DNA03_MODEL_PATH=/path/to/nucleotide-transformer-v2-100m-multi-species-snapshot \
python paper_cases/dna/DNA03/models/nt-v2/forward/code/encode_dna03_transformer_models_v020.py \
  --model nt \
  --membership constraints/dna/DNA03/refset/CORE_VIEW_MEMBERSHIP.tsv \
  --fasta constraints/dna/DNA03/refset/DNA03_FORWARD_VIEWS.fasta \
  --output work/DNA03/nt-v2/forward
```


Then call the public scorer:

```bash
python constraints/dna/DNA03/scoring/score_dna03.py \
  --embeddings work/DNA03/nt-v2/forward/EMBEDDINGS.npz \
  --repeat-embeddings work/DNA03/nt-v2/forward/REPEAT_EMBEDDINGS.npz \
  --model-id nt-v2 \
  --output work/DNA03/nt-v2/score
```
