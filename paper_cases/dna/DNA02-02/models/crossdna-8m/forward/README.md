# crossdna-8m Forward for DNA02-02

This directory preserves the historical paper Forward implementation and its current runnable reference form. To reproduce the paper execution, load the checkpoint once and perform the two eval-mode encodings without reloading it between passes. The optional `--reload-between-passes` mode is a separate robustness run and must not be identified as reproduction of the archived paper representations.

It reads the released 54,020-object FASTA through the normalized object table and writes primary and repeat NPZ representations satisfying the DNA02-02 public contract.

Current runnable reference: `scripts/encode_dna01_crossdna.py` SHA256 `8d909dab38fd499cf4cefcfe6790666544e3fc6e06dc2f515525a40b50b5afd5`.

For the executable sequence and frozen parameters, follow [../runplan/COMMANDS.md](../runplan/COMMANDS.md). The original machine-readable RunPlan is [../runplan/RUNPLAN.json](../runplan/RUNPLAN.json).
