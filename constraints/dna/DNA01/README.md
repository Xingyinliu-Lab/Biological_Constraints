# DNA01 — Reverse-complement consistency

DNA01 measures whether a frozen DNA representation exposes the reverse-complement relation through one shared, low-complexity affine operator. A user performs model inference; this release validates the submitted representations and performs every DNA01-specific calculation from affine fitting through the final confidence interval.

## Scientific definition

For a forward DNA sequence `x`, the paired object is `RC(x)`. DNA01 fits a shared affine map on four outer folds and evaluates it on the held-out dependency-cluster fold: `z(RC(x)) ≈ z(x)T + b`.

For each held-out view, `e_i` is the prediction residual and `a_i` is the median distance from its observed RC target to eight frozen, same-view RC anchors. The per-view score is `F_i = a_i/(a_i+e_i)`.

The formal `F` is the equal-weight macro mean over views, then parent loci, then dependency clusters. A 20,000-draw dependency-cluster bootstrap gives the 95% percentile interval. At a single view, `F_i=0.5` means `e_i=a_i`; it is not a chance baseline.

The current public release is self-contained in this directory through its RefSet, scoring and comparator assets.

## Contract 1: RefSet to model Forward input

Use these files from [`refset/`](refset/):

- `DNA01_OBJECT_TABLE.tsv` defines the 2,920 formal objects. `object_id` is the required external identifier; `orientation` is `FWD` or `RC`; `length_bp` and `sequence_sha256` bind the sequence.
- `records.fasta` contains exactly those 2,920 sequences and is the ready-to-use model Forward input. FASTA identifiers are the formal `object_id` values.
- Do not alter sequences, reverse-complement them again, omit objects, add model-dependent filters, or replace the released fold/dependency assignments.

The model may use its own tokenizer, checkpoint, context handling and extraction implementation. Its output for each formal object must be one fixed-dimensional, sequence-level vector. Any token-to-nucleotide mapping and pooling needed to construct that vector belongs to the model Forward and must be applied consistently to FWD and RC objects.

DNA01 also requires a deterministic repeat Forward over the same 2,920 objects. This means repeating inference and representation extraction, not retraining the model. The frozen gate requires identical dtype, shape and value bit patterns after alignment by `object_id`; signed zero and dtype changes are therefore not treated as bitwise identical.

## Contract 2: Model Forward to scorer input

Save primary and repeat outputs as separate NumPy NPZ files. Each file must contain exactly:

| Key | Shape | Requirement |
|---|---|---|
| `object_ids` | `(2920,)` | Unique strings; exact `object_id` universe from the object table. Row order may differ between files. |
| `vectors` | `(2920, d)` | Finite numeric values; `d > 0`; one sequence-level vector per ID. |

Example export:

```python
from pathlib import Path
import numpy as np

out = Path("work/DNA01/my_model")
out.mkdir(parents=True, exist_ok=True)
np.savez_compressed(out / "primary.npz", object_ids=np.asarray(object_ids), vectors=np.asarray(vectors))
np.savez_compressed(out / "repeat.npz", object_ids=np.asarray(repeat_ids), vectors=np.asarray(repeat_vectors))
```

The scorer rejects missing, duplicate or extra IDs; non-finite values; zero-dimensional arrays; inconsistent representation dimensions; and any primary/repeat mismatch.

## Run the public scorer

Prepare the lightweight scoring environment in [`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md), then run from the repository root:

```bash
python constraints/dna/DNA01/scoring/score_dna01.py \
  --embeddings work/DNA01/my_model/primary.npz \
  --repeat-embeddings work/DNA01/my_model/repeat.npz \
  --model-key my_model \
  --output work/DNA01/my_model/score
```

The six RefSet tables default to the files shipped beside the scorer. Advanced users may point the corresponding path options at another internally compatible DNA01 RefSet; the scorer validates the references and grouping required to execute the frozen calculation and records the actual file identities in its result. RefSet biological relationships are validated during RefSet materialization, not re-audited at every model score. It does not claim that an alternative RefSet is the official release.

Alpha `0.001`, five folds, the hierarchy, dynamic-range gate, bootstrap unit, 20,000 draws and seed `20261013` are frozen and cannot be overridden from the CLI.

## Outputs and interpretation

| File | Meaning |
|---|---|
| `FORMAL_RESULT.json` | Primary `F`, 95% CI, coverage/QC, diagnostics, fixed numeric settings and input identities. |
| `INDIVIDUAL_SCORES.tsv` | Canonical ID-level output: relation/view identity, dependency grouping, eight observed target-to-anchor distances, raw residual and target scale, formal `F_credit`, specificity and cycle diagnostic. |
| `PARENT_F_LEDGER.tsv` | Equal-weight per-parent aggregation on the formal F scale. |
| `DEPENDENCY_CLUSTER_F_LEDGER.tsv` | Equal-weight dependency-cluster scores used for the final macro mean. |
| `F_BOOTSTRAP.npy` | The 20,000 formal dependency-cluster bootstrap means on the F scale. |
| `OUTPUT_MANIFEST.json` | Sizes and hashes of generated result files. |
| `SHA256SUMS` | Output integrity list. |
| `RUN_COMPLETE` | Written only after atomic successful publication. |

Top-level `status=PASS` means the run completed and the inputs and calculations passed validation. It does not by itself mean that a primary score was reportable. `score_status=REPORTED` means both informative-view and informative-cluster fractions passed the frozen 0.8 gate. `score_status=INSUFFICIENT_DYNAMIC_RANGE` means execution succeeded, but the representation cannot support a formal primary score. `Q_spec` and cycle residuals are diagnostics; they do not replace or multiply `F`.

The scorer fits a dense `d × d` affine operator five times. Memory grows approximately as `d²`, and solve time approximately as `d³`; very high-dimensional representations may therefore require substantial RAM and CPU time. This is the frozen DNA01 method, not an optional probe chosen by the user.

Common input failures are reported before fitting: missing/extra/duplicate IDs, non-real or non-finite vectors, inconsistent dimensions, and non-bitwise-identical repeat output. If a run reports insufficient dynamic range, inspect `coverage` and the `evaluation_status` / `unscorable_reason` columns rather than treating it as a software failure.

## Built-in reference

The exact 4-mer reverse-complement permutation is the registered analytic white-box reference. Its method, code, provenance, QC and canonical `F=1`, 95% CI `[1,1]` result are self-contained under [`comparators/wb-4mer-rc-permutation/`](comparators/wb-4mer-rc-permutation/).
