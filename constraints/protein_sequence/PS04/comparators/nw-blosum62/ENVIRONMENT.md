# Comparator environment

This comparator owns its special runtime dependency:

- EMBOSS `needle` executable
- Python standard library
- NumPy where imported by the formal evaluator

Use EMBOSS Needle 6.6.0.0. The constraint scoring environment alone does not install this external tool.

## Operational environment guide

This environment belongs to this comparator only. Install it separately from model-forward environments unless the recorded execution explicitly used the same environment.

### Installation boundary

Use the official upstream installation instructions for external scientific software; this repository does not republish or replace upstream binaries.

Python-only dependencies can be installed in a dedicated virtual environment. External command-line tools must be discoverable on `PATH`, or supplied through the comparator's documented path argument. Do not hard-code an author-server path.

### Verification

Before a formal run:

1. confirm the external executable version;
2. run the comparator program with `--help`;
3. run on the complete released RefSet and preserve all formal IDs;
4. reject non-finite numeric output or incomplete population coverage.

The method, inputs, arguments and formal F result are documented in [`README.md`](README.md).

Generated caches, temporary alignments and intermediate score tables belong in a separate working directory. The committed result is `results/FORMAL_RESULT.json`.
