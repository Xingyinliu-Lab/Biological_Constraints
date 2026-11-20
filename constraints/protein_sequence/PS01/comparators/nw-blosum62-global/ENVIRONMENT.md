# Comparator environment

This comparator owns its special runtime dependency:

- EMBOSS `needle` 6.6.0.0 executable
- Python standard library
- NumPy where imported by the formal evaluator

Use the executable/version bound by the comparator's formal method or execution record. The constraint scoring environment alone does not install this external tool.

## Operational environment guide

This environment belongs to this comparator only. Install it separately from model-forward environments unless the recorded execution explicitly used the same environment.

### Installation boundary

The dependency names and special executable recorded above are the required public specification. Use the official upstream installation instructions for external scientific software; this repository does not republish or replace upstream binaries. Match a frozen version when one is named by the comparator method or result provenance.

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

The 2026-10-19 full-population validation used launcher SHA256
`d369c810d4c60532638db7c3c51903cd5f37a75405a3e26174ee2a4f5d380590`
and underlying `_needle` binary SHA256
`22d449c9f3ff7dd3d7e888c644f5ab75651563b47efd8de0db289d65b8828b39`.

Included Python entry/support files: `ps01_broad16k_nw_formal_wrapper_v2_0_1.py`, `ps01_formal_evaluator_v1_0_1.py`, `ps01_nw_blosum62_comparator.py`, `score_ps01_broad16k_nw_parallel.py`.

Generated caches, temporary alignments and intermediate score tables belong in a separate working directory. The canonical comparator result is [`results/FORMAL_RESULT.json`](results/FORMAL_RESULT.json); its ID-level audit outputs are [`results/INDIVIDUAL_SCORES.tsv.gz`](results/INDIVIDUAL_SCORES.tsv.gz) and [`results/RAW_MEASUREMENTS.tsv.gz`](results/RAW_MEASUREMENTS.tsv.gz).
