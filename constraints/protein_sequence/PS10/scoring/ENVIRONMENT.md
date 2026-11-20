# PS10 scoring environment

## Validated scoring baseline

- Python 3.10.20
- Standard-library modules are not listed as dependencies.

## Required third-party packages

- `numpy==1.26.4`
- `scipy==1.14.1`
- `scikit-learn==1.7.2`

## Create the environment

```bash
python -m venv .venv-ps10
# Linux/macOS
source .venv-ps10/bin/activate
# Windows PowerShell
# .venv-ps10\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install numpy==1.26.4 scipy==1.14.1 scikit-learn==1.7.2
```

## Verify the entry point

Run the official scoring command shown in the constraint README. Installing this scoring environment does not install or validate a foundation-model checkpoint.

The package versions above are the project-owned scoring baseline and are also recorded in [`ENVIRONMENT.json`](ENVIRONMENT.json). All model representations, including the paper cases, are scored in this environment. Do not let a model's CUDA environment silently choose the PS10 optimizer implementation.
