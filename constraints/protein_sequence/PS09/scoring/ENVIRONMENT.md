# PS09 scoring environment

## Supported Python

- Python 3.10 or 3.11
- Standard-library modules are not listed as dependencies.

## Required third-party packages

- `numpy==1.26.4`
- `scikit-learn==1.7.2`
- `safetensors==0.6.2` when using Safetensors input

## Create the environment

```bash
python -m venv .venv-ps09
# Linux/macOS
source .venv-ps09/bin/activate
# Windows PowerShell
# .venv-ps09\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install numpy==1.26.4 scikit-learn==1.7.2 safetensors==0.6.2
```

## Verify the entry point

Run the official scoring command shown in the constraint README. Installing this scoring environment does not install or validate a foundation-model checkpoint.

The package versions above are the project-owned scoring baseline. Do not merge incompatible model CUDA stacks into this environment unless the selected scorer explicitly requires them.
