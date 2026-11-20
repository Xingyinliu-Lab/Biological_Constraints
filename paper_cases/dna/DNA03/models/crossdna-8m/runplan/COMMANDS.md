# DNA03 / crossdna-8m commands

Use the checkpoint and environment recorded in `paper_cases/model_resources/dna/crossdna-8m/`. The adapter reads the released FASTA and membership table and writes primary and repeat NPZ files directly in the DNA03 public contract. The output directory must not already exist.

```bash
mkdir -p work/DNA03/crossdna-8m
python paper_cases/dna/DNA03/models/crossdna-8m/forward/code/encode_dna03_crossdna_v020.py \
  --membership constraints/dna/DNA03/refset/CORE_VIEW_MEMBERSHIP.tsv \
  --fasta constraints/dna/DNA03/refset/DNA03_FORWARD_VIEWS.fasta \
  --model-dir /path/to/crossdna_inference \
  --output work/DNA03/crossdna-8m/forward
```


Then call the public scorer:

```bash
python constraints/dna/DNA03/scoring/score_dna03.py \
  --embeddings work/DNA03/crossdna-8m/forward/EMBEDDINGS.npz \
  --repeat-embeddings work/DNA03/crossdna-8m/forward/REPEAT_EMBEDDINGS.npz \
  --model-id crossdna-8m \
  --output work/DNA03/crossdna-8m/score
```
