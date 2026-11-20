# PS02A — higher-order protein homology

PS02A is a frozen public scientific measurement package. A user performs model
Forward on the released protein sequences and submits one whole-sequence vector
per formal endpoint. This package performs every PS02A-specific operation after
Forward: bundle matching, safe cosine measurements, relation credits, frozen
hierarchical aggregation, anchor-clan bootstrap, sensitivity analysis and
formal output generation.

Paper models and third-party models use the same public scorer. Checkpoint and
historical execution assertions belong to `paper_cases/`; they do not change the
score calculation.

## 1. Frozen scientific definition

Each of the 10,072 formal bundles contains an anchor, a positive endpoint from
the higher-order homologous relation and a matched negative endpoint. Primary
matching uses anchor-relative sequence identity and bilateral alignment
coverage. Source-side length, composition and taxonomy diagnostics do not alter
formal membership. The complete anchor clan is the uncertainty resampling unit.

The materialized RefSet release is `PS02A_ASSAY_REFSET_RELEASE_v1.0.1_16K_20261006`; historical smaller populations are not public scoring alternatives.

## 2. Contract A: RefSet to model Forward

Forward all 13,720 records in
[`refset/endpoint_sequences.faa`](refset/endpoint_sequences.faa). The FASTA
record token after `>` is the formal object ID. Preserve it exactly.

The bundle and dependency files are scorer inputs; users do not reconstruct or
edit them:

- `FORMAL_BUNDLE_LEDGER.tsv.gz`: 10,072 anchor/positive/negative bundles;
- `POINT_DEPENDENCY_MEMBERSHIP.tsv.gz`: frozen point-estimator grouping;
- `UNCERTAINTY_ANCHOR_CLAN_MEMBERSHIP.tsv.gz`: 158 bootstrap units.

Do not filter records for context length, tokenization success or embedding
quality. A model-specific windowing or pooling implementation belongs to the
user's Forward, must be applied consistently, and must produce exactly one
whole-sequence vector for every formal object.

## 3. Contract B: model Forward to scorer input

The preferred public input is an NPZ with exactly two arrays:

| Key | Shape | Requirement |
|---|---:|---|
| `object_ids` | `(13720,)` | unique strings; exact released Forward FASTA universe |
| `vectors` | `(13720, d)` | finite numeric matrix; `d > 0` and fixed across rows |

Row order is arbitrary because IDs control alignment. Missing, duplicate or
extra IDs relative to the released 13,720-record FASTA are rejected. The frozen
10,072 bundles use 13,573 of those objects; the remaining released records are
retained because the frozen Forward FASTA and registered paper representations
share that identity. They do not enter a score silently. Finite zero-norm
vectors are accepted and use the frozen total safe-cosine value `0.0`.

For the registered paper cases, the same scorer also accepts a safetensors file
plus an `embedding_index.tsv` containing `occurrence_id` and `embedding_key`.
This is a compatibility input, not a second scoring method.

## 4. Run the official scorer

Install the small scoring environment described in
[`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md), then run from the repository
root:

```bash
python constraints/protein_sequence/PS02A/scoring/score_ps02a.py \
  --representations work/PS02A/my_model/representations.npz \
  --model-id my_model \
  --output-dir work/PS02A/my_model/score
```

For an executed safetensors paper case:

```bash
python constraints/protein_sequence/PS02A/scoring/score_ps02a.py \
  --representations work/PS02A/my_model/embeddings.safetensors \
  --representation-index work/PS02A/my_model/embedding_index.tsv \
  --model-id my_model \
  --output-dir work/PS02A/my_model/score
```

The formal CLI does not expose scientific choices. It always uses
`delta_num=0`, 20,000 PCG64 draws and seed `20261006`.

## 5. Outputs

| File | Contents |
|---|---|
| `RAW_MEASUREMENTS.tsv.gz` | bundle ID and anchor-positive/anchor-negative cosine measurements |
| `INDIVIDUAL_SCORES.tsv.gz` | bundle/object identities, raw measurements, difference, F credit and dependency IDs |
| `DEPENDENCY_AGGREGATION_LEDGER.tsv.gz` | frozen dependency-group F means, child counts, equal weights and group contributions |
| `SCORE_SUMMARY.json` | F, F 95% CI, coverage, frozen parameters, sensitivity diagnostics and input identities |
| `OUTPUT_MANIFEST.json` | bytes and SHA256 for every formal output |

F credit is `1`, `0` or `-1` according to the frozen comparison of the
positive and negative measurements. The frozen hierarchy and confidence
interval operate directly on F; resampling uses complete anchor-clan blocks.
`formal_result=true` means the formal PS02A method was used; it does not add a
third-party result to the paper registry.

## 6. Comparator and paper cases

The built-in NW–BLOSUM62 comparator, its complete execution chain and canonical
reference result live under
[`comparators/nw-blosum62-primary-comparator/`](comparators/nw-blosum62-primary-comparator/).
It uses this same RefSet and public aggregation rules.

The five registered model implementations and their public-scorer outputs live under `paper_cases/protein_sequence/PS02A/`.
Nothing in this Constraint depends on those cases.
