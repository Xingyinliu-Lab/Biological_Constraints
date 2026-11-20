# Comparator environment

This comparator owns its special runtime dependency:

- TM-align executable
- NumPy

Use a compatible TM-align executable. The constraint scoring environment alone does not install this external tool.

## Operational environment guide

This environment belongs to this comparator only. Install it separately from model-forward environments.

### Installation boundary

Use the official upstream installation instructions for external scientific software; this repository does not republish or replace upstream binaries.

Python-only dependencies can be installed in a dedicated virtual environment. External command-line tools must be discoverable on `PATH`, or supplied through the comparator's documented path argument. Do not hard-code an author-server path.

### Verification

Before a formal run:

1. record the Python version and installed package versions;
2. record the external executable version where applicable;
3. run the comparator program with `--help` when it exposes a CLI;
4. confirm that the executable accepts two coordinate-file arguments;
5. run on the complete frozen RefSet and preserve all formal IDs;
6. reject non-finite numeric output or incomplete population coverage.

The public inputs, arguments and result schema are documented in [`README.md`](README.md).

Generated caches, temporary alignments and intermediate score tables belong in a separate working directory.
