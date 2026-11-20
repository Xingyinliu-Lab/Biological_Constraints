# RNA09 Forward: rinalmo-giga

The Python files here encode the released 36 endpoint sequences and write the sequence-level NPZ required by the public scorer. Run from the repository root:

```powershell
python paper_cases/rna/RNA09/models/rinalmo-giga/forward/encode_official.py `
  --fasta constraints/rna/RNA09/refset/scoring_input_v01/RNA09_PANEL_ENDPOINTS.fasta `
  --repo <RINALMO_GIGA_REPOSITORY> `
  --out work/RNA09/rinalmo-giga/forward/REPRESENTATIONS.npz
```

Use the recorded repository revision and whole-sequence pooling rule from `paper_cases/model_resources/rna/rinalmo-giga/`. The resulting NPZ is the input to `../scoring_call/COMMANDS.md`.
