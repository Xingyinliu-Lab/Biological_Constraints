# RNA09 Forward: rna-fm

The Python files here encode the released 36 endpoint sequences and write the sequence-level NPZ required by the public scorer. Run from the repository root:

```powershell
python paper_cases/rna/RNA09/models/rna-fm/forward/encode_rnafm_global.py `
  --fasta constraints/rna/RNA09/refset/scoring_input_v01/RNA09_PANEL_ENDPOINTS.fasta `
  --weights <RNAFM_WEIGHTS> `
  --source-root <RNAFM_SOURCE_ROOT> `
  --output work/RNA09/rna-fm/forward/REPRESENTATIONS.npz
```

Use the checkpoint identity and whole-sequence pooling rule recorded in `paper_cases/model_resources/rna/rna-fm/`. The resulting NPZ is the input to `../scoring_call/COMMANDS.md`.
