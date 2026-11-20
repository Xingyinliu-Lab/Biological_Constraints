# RNA06 — paired-partner specificity

RNA06 measures whether nucleotide-addressable model representations distinguish each true paired partner from frozen paired nonpartners. The released formal result uses the frozen primary paired-nonpartner population.

## 1. RefSet to model Forward input

Encode the 105 sequences in [`refset/rna_sequences.fasta`](refset/rna_sequences.fasta), preserving the sequence SHA256 identifiers and zero-based nucleotide coordinates. The full primary population is fixed by:

- 3,176 directed endpoints in `RNA06_FORMAL_DIRECTED_ENDPOINT_LEDGER.tsv`;
- 14,948 paired-nonpartner candidates in `RNA06_FORMAL_PAIRED_NONPARTNER_CANDIDATE_LEDGER.tsv`;
- 1,588 reciprocal edge bundles in `RNA06_FORMAL_RECIPROCAL_EDGE_BUNDLE_LEDGER.tsv`.
- frozen positive-pair identities in `RNA06_FROZEN_PAIR_STATE_MAPS.tsv`.

Do not replace this population with the smaller Vienna-compatible comparator subset.

## 2. Model Forward output contract

Write one NPZ containing exactly:

| Key | Shape | Requirement |
|---|---|---|
| `sequence_sha256` | `(n,)` | released sequence IDs |
| `position0` | `(n,)` | zero-based nucleotide coordinates |
| `vectors` | `(n, d)` | finite numeric matrix with `d > 0` |

Every position referenced as an anchor, true partner or paired nonpartner must be present exactly once. Row order is arbitrary; IDs and vector rows must remain aligned. Historical sequence-keyed PyTorch matrices are accepted only as a compatibility input for registered paper cases.

### Position mapping evidence

`position0` is a biological nucleotide coordinate, not a model token index. For every anchor, true partner and paired nonpartner target, users must map the released `(sequence_sha256, position0)` to their model tokens and retain:

- model, checkpoint and tokenizer identities;
- token-to-nucleotide spans, including k-mer and overlapping-token behavior;
- BOS/EOS and other special-token handling;
- windowing, overlap reconciliation and restoration of original coordinates;
- the pooling or readout rule when a nucleotide maps to multiple tokens or a token spans multiple nucleotides.

Anchor, partner and nonpartner vectors must each be extracted from the sequence and coordinate named by the corresponding frozen ledger row; they must not be inferred by row order or reused across different targets. The public scorer validates submitted IDs, coordinate coverage, input structure and numerical values against the released scoring tables. Biological sequence, pair-state and dependency correctness is established when the RefSet is materialized and frozen; it is not re-audited during each scoring run. The scorer does not load or independently verify third-party tokenizers. Mapping evidence is Forward provenance, not an additional scorer input.

## 3. Run the official scorer

Install the lightweight scoring environment in [`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md), then run:

```bash
python constraints/rna/RNA06/scoring/score_rna06.py \
  --vectors work/RNA06/my_model/REPRESENTATIONS.npz \
  --model my_model \
  --output work/RNA06/my_model/score
```

The scorer automatically loads the full formal RefSet. It performs the frozen sum-plus-Hadamard feature construction, leave-one-cluster-out fitting, four-fold inner lambda selection, reciprocal hierarchy aggregation and 20,000-draw cluster bootstrap. Users do not choose formal scoring parameters.

## 4. Formal outputs

| File | Contents |
|---|---|
| `HELDOUT_PREDICTIONS.tsv` | true/candidate sample IDs, raw held-out scores and hierarchy IDs |
| `CANDIDATE_SCORES.tsv` | one row per frozen nonpartner candidate with true/nonpartner scores, margin and `F_credit` |
| `INDIVIDUAL_SCORES.tsv` | one `individual_F` and raw true/nonpartner measurements per directed endpoint |
| `EDGE_SCORE_LEDGER.tsv` | reciprocal `edge_F` scores |
| `PARENT_SCORE_LEDGER.tsv` | parent-level `parent_F` scores |
| `CLUSTER_SCORE_LEDGER.tsv` | `cluster_F` scores used for F and its CI |
| `FORMAL_RESULT.json` | F, F 95% CI, coverage, selected lambdas and method identity |
| `BOOTSTRAP_MEANS.npy` | frozen F bootstrap distribution |
| `OUTPUT_MANIFEST.json` | output sizes and SHA256 values |

F is the sole formal reporting score. `formal_result=true` means the released protocol was used; it does not add a third-party result to the paper registry.

## 5. Comparator and paper cases

The ViennaRNA partner-probability reference, its compatible scope, code and registered result are under [`comparators/viennarna-partner-probability/`](comparators/viennarna-partner-probability/). The five paper-model Forward implementations and outputs are under `paper_cases/rna/RNA06/`; they call this same scorer.
