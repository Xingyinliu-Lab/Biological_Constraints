# RNA09 Forward: orthrus-4-track

The Python files here encode the released 36 endpoint sequences and write the sequence-level NPZ required by the public scorer. Run from the repository root:

```powershell
python paper_cases/rna/RNA09/models/orthrus-4-track/forward/encode_orthrus_4track.py `
  --fasta constraints/rna/RNA09/refset/scoring_input_v01/RNA09_PANEL_ENDPOINTS.fasta `
  --mode global `
  --output work/RNA09/orthrus-4-track/forward/REPRESENTATIONS.npz
```

The default `--model-id antichronology/orthrus-4-track` is the registered checkpoint; record any override. The resulting NPZ is the input to `../scoring_call/COMMANDS.md`.
