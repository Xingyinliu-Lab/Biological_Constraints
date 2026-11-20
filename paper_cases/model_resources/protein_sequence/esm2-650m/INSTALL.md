# Install and load ESM-2 650M

This model-owned guide records how to recreate the ESM-2 runtime and load the checkpoint used by the formal executions. Constraint-specific RefSet inputs, mapping, pooling and output files are documented only in the matching `results/<domain>/<constraint_id>/models/<model_id>/forward/IMPLEMENTATION.md` manual.

## Executed environment

The recorded environment used Python 3.12.13, PyTorch 2.11.0 and Transformers 4.57.6. The complete observed package inventory is available in [`ENVIRONMENT.json`](ENVIRONMENT.json) and [`ENVIRONMENT.md`](ENVIRONMENT.md).

### Windows PowerShell

```powershell
py -3.12 -m venv .venv-esm2
.\.venv-esm2\Scripts\python.exe -m pip install --upgrade pip
.\.venv-esm2\Scripts\python.exe -m pip install numpy==2.5.3 scipy==1.18.1 transformers==4.57.6 torch==2.11.0 fair-esm
```

Verify with the same interpreter:

```powershell
.\.venv-esm2\Scripts\python.exe -c "import esm, numpy, scipy, torch, transformers; print(torch.__version__, torch.cuda.is_available(), transformers.__version__)"
```

### Linux or macOS

```bash
python3.12 -m venv .venv-esm2
source .venv-esm2/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy==2.5.3 scipy==1.18.1 transformers==4.57.6 torch==2.11.0 fair-esm
python -c "import esm, numpy, scipy, torch, transformers; print(torch.__version__, torch.cuda.is_available(), transformers.__version__)"
```

For a CUDA-specific PyTorch wheel, use the installation selector supplied by the PyTorch project for the local driver/runtime, then confirm the installed PyTorch version and CUDA availability. The repository records the environment used for the formal runs; it does not maintain or replace the upstream model software.

## Executed checkpoint

The formal executions used:

```text
esm2_t33_650M_UR50D.pt
SHA256 ea9d0522b335a8778dea6535a65301f10208dece28cd5865482b0b1fc446168c
```

The recorded loading interface is:

```python
import esm

model, alphabet = esm.pretrained.load_model_and_alphabet_local(checkpoint_path)
model.eval()
```

Verify the checkpoint before claiming reproduction of a registered result.

### Windows PowerShell

```powershell
Get-FileHash -Algorithm SHA256 .\esm2_t33_650M_UR50D.pt
```

### Linux or macOS

```bash
sha256sum esm2_t33_650M_UR50D.pt
```

The exact upstream source-code commit was not recorded in the formal execution evidence. This limitation is stated in [`checkpoint_identity.json`](checkpoint_identity.json); do not invent a revision. Equivalent underlying model parameters may be stored in a different serialization format without producing identical file bytes, but reproduction claims must still document the exact asset and extraction semantics used.

## Responsibility boundary

- [`ENVIRONMENT.md`](ENVIRONMENT.md) records the human-readable executed environment.
- [`ENVIRONMENT.json`](ENVIRONMENT.json) records the machine-readable package inventory.
- [`checkpoint_identity.json`](checkpoint_identity.json) binds the official source, loader and executed checkpoint identity.
- Constraint-specific forward instructions remain in `results/<domain>/<constraint>/models/<model>/forward/IMPLEMENTATION.md`.

Model installation does not define a constraint's representation interface and does not run its scorer.
