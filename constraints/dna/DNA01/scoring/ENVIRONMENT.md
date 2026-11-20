# DNA01 scoring environment

## Supported Python

- Python 3.10, 3.11 or 3.12
- Standard-library modules are not listed as dependencies.

## Required third-party packages

- `numpy==1.26.4`

## Create the environment

```bash
python -m venv .venv-dna01
# Linux/macOS
source .venv-dna01/bin/activate
# Windows PowerShell
# .venv-dna01\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install numpy==1.26.4
```

## Verify the entry point

Run the official scoring command shown in the constraint README. The built-in 4-mer comparator uses the same environment. Installing this scoring environment does not install or validate a foundation-model checkpoint.

The package versions above are the project-owned scoring baseline. Do not merge incompatible model CUDA stacks into this environment unless the selected scorer explicitly requires them.
