# Forward

From the repository root, run the released checkpoint with an explicit path:

```bash
python paper_cases/protein_sequence/PS10/models/esm2-650m/forward/code/forward_peptide_embeddings.py \
  --model esm2 --checkpoint /path/to/esm2_t33_650M_UR50D.pt \
  --records constraints/protein_sequence/PS10/refset/PEPTIDE_REFSET_RECORDS.tsv \
  --output-dir work/PS10/esm2-650m/forward
python constraints/protein_sequence/PS10/scoring/score_ps10.py \
  --embeddings work/PS10/esm2-650m/forward/embeddings.npz --model-id esm2-650m \
  --output-dir work/PS10/esm2-650m/score
```
