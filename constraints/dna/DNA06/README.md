# DNA06 — Codon-local translation-consequence readability

DNA06 measures whether a model's local representation distinguishes ordinary missense changes from synonymous changes at matched coding sites. The released Constraint has one formal readout: the representation extracted from the three nucleotide positions of the affected codon.

The user performs model Forward. DNA06 validates the submitted rows, constructs WT/ALT difference features, runs the nested cross-fitted ridge classifier, computes held-out variant predictions and matched-site F credits, aggregates through parent CDS and dependency clusters, and reports F with its 95% confidence interval.

## 2. RefSet to model Forward input

`refset/records.fasta` contains the 11,935 formal sequence objects. Preserve each FASTA object ID. [`refset/FORWARD_TARGETS.tsv`](refset/FORWARD_TARGETS.tsv) is the direct Forward map: every required output row, its source object and sequence identity, and its canonical CDS codon interval. It contains 3,920 WT rows and 11,506 ALT rows, for 15,426 rows total. Coordinates are 0-based and half-open: `[codon_start_0based, codon_end_exclusive)`.

For every target row, extract a representation from the three nucleotide positions in its released codon interval. Token-to-nucleotide mapping, window stitching and pooling belong to model Forward. The Constraint defines biological coordinates, not a model's internal token indices.

For every reported evaluation, retain mapping evidence describing the model and tokenizer identity, token-to-nucleotide/span rule, windowing, treatment of boundary-crossing tokens and final codon-local pooling/readout. The public scorer verifies released row IDs and numeric inputs, but does not load or independently validate a third-party tokenizer.

## 3. Forward output contract

Save one NPZ containing exactly:

| Key | Shape | Requirement |
|---|---|---|
| `row_ids` | `(15426,)` | unique strings covering all 3,920 WT site rows and 11,506 ALT variant rows |
| `embeddings` | `(15426, d)` | finite real integer or floating-point codon-local vectors with `d > 0`; complex, boolean, string and object arrays are rejected |

Required row IDs are:

- WT: `<site_id>|WT|<parent_object_id>`
- ALT: `<site_id>|ALT|<variant_id>|<variant_object_id>`

Order is arbitrary when IDs and rows remain paired. Missing, duplicate or extra rows; non-finite values; and zero-dimensional embeddings are rejected.

The machine-readable contract is [`scoring/INPUT_CONTRACT.json`](scoring/INPUT_CONTRACT.json). It does not require a model-specific token index convention; it requires each submitted row to correspond to the released biological target in `FORWARD_TARGETS.tsv`.

## 4. Run the official scorer

Install the scoring-only dependencies in [`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md), then run:

```bash
python constraints/dna/DNA06/scoring/score_dna06.py \
  --embeddings work/DNA06/my_model/EMBEDDINGS.npz \
  --model-id my_model \
  --output-dir work/DNA06/my_model/score
```

All folds, the nine-value lambda grid, tie tolerance, aggregation levels, 20,000 bootstrap draws and seed `20261014` are fixed by the release.

## 5. Outputs

| File | Contents |
|---|---|
| `INDIVIDUAL_SCORES.tsv` | 11,506 held-out variant predictions with variant, site, parent, dependency cluster, fold, truth label and decision score |
| `PAIR_CREDIT_LEDGER.tsv` | matched missense–synonymous margins, tie status and relation `F_credit` |
| `SITE_SCORE_LEDGER.tsv` | 3,920 matched-site F scores |
| `PARENT_SCORE_LEDGER.tsv` | 560 parent-CDS F contributions before dependency-cluster aggregation |
| `CLUSTER_SCORE_LEDGER.tsv` | 554 dependency-cluster F scores used for aggregation and bootstrap |
| `FORMAL_RESULT.json` | F, its 95% CI, selected lambda per outer fold, coverage and representation identity |
| `OUTPUT_MANIFEST.json` | byte sizes and SHA256 identities of formal outputs |

Formal F is evaluated end-to-end: relation credits are `1` (win), `0` (tie), or `-1` (loss), and every higher aggregation is the prescribed mean of its direct F children. A score produced by the formal protocol is not automatically registered as a paper result.

## 6. Comparator and paper cases

The Standard Genetic Code Table 1 reference and raw-codon ridge reference, including executable methods and canonical outputs, are contained in `comparators/standard-genetic-code-table1/`.

The five registered model implementations and public-scorer outputs are under `paper_cases/dna/DNA06/`. DNA06 does not depend on those cases.
