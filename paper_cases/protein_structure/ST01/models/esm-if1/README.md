# ST01 / esm-if1

This paper case records the esm-if1 Forward implementation and its use of the unique public ST01 scorer. Scientific definitions, RefSet assets and scoring code remain in [`constraints/protein_structure/ST01`](../../../../../constraints/protein_structure/ST01/README.md).

- Forward implementation: [`forward/README.md`](forward/README.md)
- Public scoring command: [`scoring_call/COMMANDS.md`](scoring_call/COMMANDS.md)
- Reconstructed public output: [`output/SCORE_SUMMARY.json`](output/SCORE_SUMMARY.json)

The model embedding is not duplicated in the repository. Recreate it from the released RefSet using `forward/scripts/`, then score the generated NPZ with the public command; no historical embedding-file SHA256 is asserted.
