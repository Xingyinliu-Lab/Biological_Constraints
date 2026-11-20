# Forward

From the repository root, run the released checkpoint with an explicit path:

```bash
python paper_cases/protein_sequence/PS10/models/esmc-600m/forward/code/forward_peptide_embeddings.py \
  --model esmc --checkpoint /path/to/biohub_ESMC-600M \
  --records constraints/protein_sequence/PS10/refset/PEPTIDE_REFSET_RECORDS.tsv \
  --output-dir work/PS10/esmc-600m/forward
python constraints/protein_sequence/PS10/scoring/score_ps10.py \
  --embeddings work/PS10/esmc-600m/forward/embeddings.npz --model-id esmc-600m \
  --output-dir work/PS10/esmc-600m/score
```
