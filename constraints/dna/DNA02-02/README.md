# DNA02-02 — single-edit adjacency

DNA02-02 tests whether a representation preserves the local geometry of two successive single-nucleotide edits. Each frozen ladder contains three equal-length sequences satisfying

\[
H(x_0,x_1)=1,\qquad H(x_1,x_2)=1,\qquad H(x_0,x_2)=2.
\]

The released Constraint is self-contained from model Forward output onward. Users run their own model; this package validates the submitted representations, calibrates the frozen repeat tolerance, computes all ladder contrasts and hierarchical scores, and reports F and F_95_CI and coverage.

## Formal scoring interface

DNA02-02 uses one formal score throughout the public interface. At the evaluation-unit level, an ordering win, tie and loss have `F_credit` values of `1`, `0` and `-1`. These are averaged through ladders, first-edit queries, loci and outer clusters to final `F`; the 20,000 outer-cluster bootstrap is reported as `F_95_CI`.

## Contract 1: RefSet to model Forward

Forward every record in [refset/DNA02_02_SEQUENCE_OBJECTS.fasta](refset/DNA02_02_SEQUENCE_OBJECTS.fasta). The FASTA identifier is the stable object ID; do not rename, omit or add objects. [refset/SEQUENCE_OBJECT.tsv](refset/SEQUENCE_OBJECT.tsv) gives the expected length and sequence SHA256 for each object.

Run the same fixed checkpoint and extraction protocol twice as separately executed Forward passes. Do not copy one output file to represent the repeat. The model must be in evaluation mode; retain the checkpoint identity, inference configuration, random-state controls, tokenization, windowing and pooling for both passes. Reloading the checkpoint between the two passes is permitted as an additional loading-path robustness check, but is not required for formal scoring. The representation must describe the complete sequence object; its dimension may be model-specific.

## Contract 2: model Forward to scorer

Submit a primary and repeat NPZ. Each file must contain exactly:

| Key | Shape | Requirement |
|---|---|---|
| `object_ids` | `(54020,)` | Unique strings exactly matching the frozen FASTA IDs; row order may differ between files. |
| `vectors` | `(54020, d)` | Finite real numeric sequence-level representations with `d > 0`; complex and boolean arrays are rejected. |

The two files must have the same object universe and representation dimension. The repeat is used only by the frozen stratified tau-calibration rule; users do not choose tau.

## Run the public scorer

Prepare the small scoring environment in [scoring/ENVIRONMENT.md](scoring/ENVIRONMENT.md), then run from the repository root:

```bash
python constraints/dna/DNA02-02/scoring/score_dna0202.py \
  --embeddings work/DNA02-02/my_model/primary.npz \
  --repeat-embeddings work/DNA02-02/my_model/repeat.npz \
  --model-id my_model \
  --output work/DNA02-02/my_model/score
```

The scorer automatically uses the released ladders, folds, object table, FASTA and tau controls. Advanced path overrides are available through `--help`; any result produced with alternative assets must retain those input identities and must not be represented as the released RefSet result.

## Outputs

| File | Contents |
|---|---|
| `INDIVIDUAL_SCORES.tsv` | One row per ladder: four raw cosine measurements, left/right contrasts, tie-aware endpoint credits and formal ladder credit. |
| `QUERY_SCORE_LEDGER.tsv` | First-edit-query aggregation with locus and cluster IDs. |
| `LOCUS_SCORE_LEDGER.tsv` | Parent-locus aggregation and outer-cluster mapping. |
| `CLUSTER_SCORE_LEDGER.tsv` | Outer-cluster scores used for the formal macro-average and bootstrap. |
| `TAU_EVIDENCE.json` | Per-length repeat-calibration evidence and selected frozen-rule tau. |
| `FORMAL_RESULT.json` | F and F_95_CI, coverage, QC, method settings and input identities. |
| `OUTPUT_MANIFEST.json` | Output file sizes and SHA256 values. |

`formal_scoring_protocol=true` means the frozen DNA02-02 computation was used. It does not automatically register a third-party result as a paper result.

## Comparator and paper cases

The built-in deterministic Hamming Level-0 reference, its code and canonical outputs are in [comparators/hamming-level0/](comparators/hamming-level0/). Paper model Forward implementations and their executions belong under `paper_cases/dna/DNA02-02/`; they call this same public scorer and are not dependencies of this Constraint.
