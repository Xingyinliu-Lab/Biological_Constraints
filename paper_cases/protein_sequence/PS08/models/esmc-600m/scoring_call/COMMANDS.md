# PS08 esmc-600m public scoring call

```bash
python paper_cases/protein_sequence/PS08/shared_forward/export_safetensors_to_npz.py \
  --embeddings work/PS08/esmc-600m/embeddings.safetensors \
  --index work/PS08/esmc-600m/embedding_index.tsv \
  --population constraints/protein_sequence/PS08/refset/PLANNED_POPULATION_MANIFEST.json \
  --output work/PS08/esmc-600m/REPRESENTATIONS.npz

python constraints/protein_sequence/PS08/scoring/representations_to_leaf.py \
  --population constraints/protein_sequence/PS08/refset/PLANNED_POPULATION_MANIFEST.json \
  --objects-fasta constraints/protein_sequence/PS08/refset/ENDPOINT_SEQUENCES.faa.gz \
  --representations work/PS08/esmc-600m/REPRESENTATIONS.npz \
  --output work/PS08/esmc-600m/FORMAL_LEAF_SCORES.tsv

python constraints/protein_sequence/PS08/scoring/score_ps08_formal.py \
  --population constraints/protein_sequence/PS08/refset/PLANNED_POPULATION_MANIFEST.json \
  --objects-fasta constraints/protein_sequence/PS08/refset/ENDPOINT_SEQUENCES.faa.gz \
  --leaf-ledger work/PS08/esmc-600m/FORMAL_LEAF_SCORES.tsv \
  --outdir work/PS08/esmc-600m/score
```
