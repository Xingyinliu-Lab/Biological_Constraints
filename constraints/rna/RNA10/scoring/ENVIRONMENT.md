# RNA10 scoring environment

## Supported Python

- Python 3.10 (the release replay used Python 3.10.20)
- Standard-library modules are not listed as dependencies.

## Required third-party packages

- `numpy==1.22.4`
- `scipy==1.8.0`
- NumPy/SciPy wheels linked to OpenBLAS ILP64 (`openblas64_`)

## Create the environment

```bash
python -m venv .venv-rna10
# Linux/macOS
source .venv-rna10/bin/activate
# Windows PowerShell
# .venv-rna10\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install numpy==1.22.4 scipy==1.8.0
```

## Verify the entry point

Run the official scoring command shown in the constraint README. Installing this scoring environment does not install or validate a foundation-model checkpoint.

SciPy 1.8.0 is part of the frozen RNA10 numerical method because the
weighted L2 logistic-ridge fits use its L-BFGS-B implementation and default
stopping tolerances. Use this scoring environment independently of the model
Forward environment.
