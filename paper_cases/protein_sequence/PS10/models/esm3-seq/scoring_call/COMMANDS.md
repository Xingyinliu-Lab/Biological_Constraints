# Public scoring call

```bash
python paper_cases/protein_sequence/PS10/models/esm3-seq/forward/code/ps10_esm3_handoff.py export-npz \
  --records constraints/protein_sequence/PS10/refset/PEPTIDE_REFSET_RECORDS.tsv \
  --embeddings work/PS10/esm3-seq/forward/raw/embeddings.safetensors \
  --index work/PS10/esm3-seq/forward/raw/embedding_index.tsv \
  --output work/PS10/esm3-seq/forward/embeddings.npz

python constraints/protein_sequence/PS10/scoring/score_ps10.py \
  --embeddings work/PS10/esm3-seq/forward/embeddings.npz \
  --model-id esm3-seq \
  --output-dir work/PS10/esm3-seq/score
```
