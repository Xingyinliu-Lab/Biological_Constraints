# RNA01-02 scoring environment

## Supported Python

- Python 3.10 or 3.11
- Standard-library modules are not listed as dependencies.

## Required third-party packages

- `numpy==1.26.4`

## Create the environment

```bash
python -m venv .venv-rna01-02
# Linux/macOS
source .venv-rna01-02/bin/activate
# Windows PowerShell
# .venv-rna01-02\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install numpy==1.26.4
```

## Verify the entry point

Run the official scoring command shown in the constraint README. Installing this scoring environment does not install or validate a foundation-model checkpoint.

The package versions above are the project-owned scoring baseline. Do not merge incompatible model CUDA stacks into this environment unless the selected scorer explicitly requires them.
