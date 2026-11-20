# DNA03 — Core-local representation

DNA03 measures whether a model's nucleotide-coordinate representation preserves the identity of a fixed DNA core when that core is observed in different surrounding contexts. The released result is the core-local formal score F.

This directory is a self-contained measurement package. A user performs model Forward only; DNA03 validates the returned representations, calibrates the repeat tolerance, computes relation-level F credits, aggregates by dependency cluster and transition, and reports F with its 95% confidence interval.

DNA03 has one formal readout: **core-local nucleotide-coordinate representation**.

## 2. RefSet to model Forward input

Forward every record in `refset/DNA03_FORWARD_VIEWS.fasta`.

The FASTA contains exactly 1,400 formal view objects. The token after `>` is the immutable `view_id`. Do not rename, omit, add or deduplicate records. `refset/CORE_VIEW_MEMBERSHIP.tsv` supplies the core interval for every view through `core_offset_start_0based` and `core_offset_end_exclusive`.

For each view, extract one vector from the frozen core interval. Token-to-nucleotide mapping, window stitching and pooling needed to form that core-local vector are part of the model Forward implementation. They must preserve the published coordinates.

For every reported evaluation, retain model and tokenizer identity, token-to-nucleotide/span rule, windowing procedure, treatment of tokens crossing the core boundary, and the final core-vector pooling/readout rule. The public scorer verifies released `view_id` coverage and numerical inputs; it does not load or independently validate a third-party tokenizer.

DNA03 also requires an independent repeat of the same Forward and extraction procedure. This is repeat inference, not repeat training.

## 3. Forward output contract

Save both Forward runs as NPZ files with exactly these arrays:

| Key | Shape | Requirement |
|---|---|---|
| `view_ids` | `(1400,)` | unique strings matching the complete FASTA ID universe |
| `embeddings` | `(1400, d)` | finite numeric core-local vectors, with `d > 0` |

Row order is arbitrary when each row remains paired with its `view_id`. Both NPZ files must contain the same IDs and representation dimension. Extra or missing IDs, duplicate IDs, NaN/Inf values and zero-dimensional embeddings are rejected.

## 4. Run the official scorer

Prepare the small scoring environment described in [`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md), then run:

```bash
python constraints/dna/DNA03/scoring/score_dna03.py \
  --embeddings work/DNA03/my_model/primary_embeddings.npz \
  --repeat-embeddings work/DNA03/my_model/repeat_embeddings.npz \
  --model-id my_model \
  --output work/DNA03/my_model/score
```

The scorer automatically uses the released queries, negative relations and 128 frozen tau controls. Bootstrap settings are fixed at 20,000 dependency-cluster draws with seed `20261014`; they are not user-selectable formal parameters.

## 5. Outputs

| File | Contents |
|---|---|
| `INDIVIDUAL_SCORES.tsv` | 35,840 relation-level rows with stable IDs, involved views, raw positive/negative cosine distances, signed margin, tie flag, calibrated tau and `F_credit` |
| `QUERY_SCORE_LEDGER.tsv` | 1,120 query-level F scores before higher-level aggregation |
| `CLUSTER_TRANSITION_SCORE_LEDGER.tsv` | dependency-cluster × transition F scores used for the formal macro-average and bootstrap |
| `TRANSITION_SCORE_LEDGER.tsv` | transition-level F scores whose macro-average is formal F |
| `TAU_EVIDENCE.json` | primary/repeat control distances and contrasts, their absolute contrast differences, and selected tau |
| `FORMAL_RESULT.json` | F, its 95% CI, transition F scores, coverage, QC and input identities |
| `OUTPUT_MANIFEST.json` | byte sizes and SHA256 values for the formal outputs |

Formal F is evaluated end-to-end: relation credits are `1` (win), `0` (tie), or `-1` (loss), and every higher aggregation is the prescribed mean of its direct F children. A result produced by this formal protocol is not automatically a registered paper result; registration additionally belongs to the corresponding paper case.

## 6. Comparator and paper cases

The built-in exact-core-containment Level-0 reference, its executable method and canonical results are under `comparators/exact-core-containment-level0/`.

The five paper model implementations and registered outputs are under `paper_cases/dna/DNA03/`. They call this same public scorer; DNA03 does not depend on those cases.
