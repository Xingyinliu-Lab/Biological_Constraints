# RNA06 scoring environment

## Supported Python

- Python 3.10 or 3.11
- Standard-library modules are not listed as dependencies.

## Required packages for the released NPZ interface

- `numpy==1.26.4`
- `scipy==1.11.4`

`torch==2.5.1` is needed only when reading the historical PyTorch tensor-file
compatibility input. It is not required for the released NPZ contract.

## Create the environment

```bash
python -m venv .venv-rna06
# Linux/macOS
source .venv-rna06/bin/activate
# Windows PowerShell
# .venv-rna06\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install numpy==1.26.4 scipy==1.11.4
```

## Verify the entry point

Run the official scoring command shown in the constraint README. Installing this scoring environment does not install or validate a foundation-model checkpoint.

The package versions above are the project-owned scoring baseline. Do not merge incompatible model CUDA stacks into this environment unless the selected scorer explicitly requires them.
