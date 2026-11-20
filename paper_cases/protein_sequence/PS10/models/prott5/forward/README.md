# Forward

From the repository root, run the released checkpoint with an explicit path:

```bash
python paper_cases/protein_sequence/PS10/models/prott5/forward/code/forward_peptide_embeddings.py \
  --model prott5 --checkpoint /path/to/prot_t5_xl_half_uniref50-enc \
  --records constraints/protein_sequence/PS10/refset/PEPTIDE_REFSET_RECORDS.tsv \
  --output-dir work/PS10/prott5/forward
python constraints/protein_sequence/PS10/scoring/score_ps10.py \
  --embeddings work/PS10/prott5/forward/embeddings.npz --model-id prott5 \
  --output-dir work/PS10/prott5/score
```
