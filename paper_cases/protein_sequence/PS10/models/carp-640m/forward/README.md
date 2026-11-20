# Forward

From the repository root, point `--checkpoint` at `TORCH_HOME/hub/checkpoints/carp_640M.pt` (or provide that cache root with `--torch-home`):

```bash
python paper_cases/protein_sequence/PS10/models/carp-640m/forward/code/forward_peptide_embeddings.py \
  --model carp640 --checkpoint /path/to/torch_cache/hub/checkpoints/carp_640M.pt \
  --torch-home /path/to/torch_cache \
  --records constraints/protein_sequence/PS10/refset/PEPTIDE_REFSET_RECORDS.tsv \
  --output-dir work/PS10/carp-640m/forward
python constraints/protein_sequence/PS10/scoring/score_ps10.py \
  --embeddings work/PS10/carp-640m/forward/embeddings.npz --model-id carp-640m \
  --output-dir work/PS10/carp-640m/score
```
