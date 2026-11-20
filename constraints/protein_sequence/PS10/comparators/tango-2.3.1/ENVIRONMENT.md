# Comparator environment

This comparator owns its special runtime dependency:

- TANGO 2.3.1 executable
- NumPy
- scikit-learn

Use the executable/version bound by the comparator's formal method or execution record. The constraint scoring environment alone does not install this external tool.

## Operational environment guide

This environment belongs to this comparator only. Install it separately from model-forward environments unless the recorded execution explicitly used the same environment.

### Installation boundary

The dependency names and special executable recorded above are the required public specification. Use the official upstream installation instructions for external scientific software; this repository does not republish or replace upstream binaries. Match the version named by the comparator method.

Python-only dependencies can be installed in a dedicated virtual environment. External command-line tools must be discoverable on `PATH`, or supplied through the comparator's documented path argument. Do not hard-code an author-server path.

### Verification

Before a formal run:

1. record the Python version and installed package versions;
2. record the external executable version where applicable;
3. run the comparator program with `--help` when it exposes a CLI;
4. run any retained method-specific tests or smoke checks;
5. run on the complete frozen RefSet and preserve all formal IDs;
6. reject non-finite numeric output or incomplete population coverage.

The full method, inputs, arguments and registered result are documented in [`README.md`](README.md).

Included Python entry/support files: `run_tango_ps10.py`, `score_ps10_tango.py`.

Generated caches, temporary alignments and intermediate score tables belong in a separate working directory. The canonical comparator result is selected by `registry/results.tsv`.
