# DNA03 / janusdna commands

Use the source, config, checkpoint and environment recorded in `paper_cases/model_resources/dna/janusdna/`. The adapter reads the released FASTA and membership table and writes primary and repeat NPZ files directly in the DNA03 public contract. The output directory must not already exist.

```bash
mkdir -p work/DNA03/janusdna
python paper_cases/dna/DNA03/models/janusdna/forward/code/encode_dna03_janusdna_v020.py \
  --membership constraints/dna/DNA03/refset/CORE_VIEW_MEMBERSHIP.tsv \
  --fasta constraints/dna/DNA03/refset/DNA03_FORWARD_VIEWS.fasta \
  --source /path/to/JanusDNA \
  --config /path/to/config.yaml \
  --checkpoint /path/to/last.ckpt \
  --output work/DNA03/janusdna/forward
```


Then call the public scorer:

```bash
python constraints/dna/DNA03/scoring/score_dna03.py \
  --embeddings work/DNA03/janusdna/forward/EMBEDDINGS.npz \
  --repeat-embeddings work/DNA03/janusdna/forward/REPEAT_EMBEDDINGS.npz \
  --model-id janusdna \
  --output work/DNA03/janusdna/score
```
