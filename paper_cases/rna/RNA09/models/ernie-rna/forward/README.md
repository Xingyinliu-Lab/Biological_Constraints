# RNA09 Forward: ernie-rna

The Python files here encode the released 36 endpoint sequences and write the sequence-level NPZ required by the public scorer. Run from the repository root:

```powershell
python paper_cases/rna/RNA09/models/ernie-rna/forward/encode_ernierna_global.py `
  --fasta constraints/rna/RNA09/refset/scoring_input_v01/RNA09_PANEL_ENDPOINTS.fasta `
  --model-dir <ERNIE_RNA_CHECKPOINT> `
  --output work/RNA09/ernie-rna/forward/REPRESENTATIONS.npz
```

Use the checkpoint identity and nucleotide-state pooling recorded in `paper_cases/model_resources/rna/ernie-rna/`. The resulting NPZ is the input to `../scoring_call/COMMANDS.md`.
