# PS02B — higher-order protein homology, topology/HSF relation

PS02B is a frozen public scientific measurement package. Users encode the
released protein-domain sequences with their own model. This package then
performs all constraint-specific work: relation matching, total safe-cosine
measurements, frozen relation credits, hierarchical aggregation, dependency
bootstrap and formal output generation.

The five paper models and any third-party model use the same public scorer.
Paper-model Forward implementations are recorded separately in `paper_cases/`;
they never select the scoring rule.

## 1. Frozen release

The active population is the single large release
`PS02B_ASSAY_REFSET_RELEASE_v1.0.0_L_20261006`. It contains 8,208 formal
anchor/positive/negative relations, 2,764 anchors, 3,905 scored objects and 351
frozen dependency/bootstrap groups. Historical smaller populations are not
public alternatives.

Each relation asks whether the model ranks the designated higher-order
homology positive above its model-blind matched negative under the topology/HSF
scheme frozen in the Scientific Word and relation ledger. Matching diagnostics
are inputs to population construction, not user-tunable score parameters.

## 2. Contract A: RefSet to model Forward

Forward every one of the 3,937 records in
[`refset/endpoint_sequences.faa`](refset/endpoint_sequences.faa), preserving the
FASTA record token as `domain_id`. The relation ledger uses 3,905 of these
objects. The additional 32 records remain part of the frozen Forward FASTA and
registered representation identity but do not silently enter a score.

Users do not rebuild or edit:

- `FORMAL_RELATION_LEDGER.tsv.gz` — relation IDs and frozen roles;
- `FORMAL_OBJECT_LEDGER.tsv.gz` — object, sequence, topology/HSF and parent data;
- `POINT_DEPENDENCY_MEMBERSHIP.tsv.gz` — point-estimator grouping;
- `UNCERTAINTY_GROUP_MEMBERSHIP.tsv.gz` — 351 resampling groups.

Model tokenization, windowing and pooling belong to Forward and must be applied
consistently. The submitted output must contain one whole-domain vector for
every released FASTA record.

## 3. Contract B: Forward output to scorer

Preferred NPZ input:

| Key | Shape | Requirement |
|---|---:|---|
| `object_ids` | `(3937,)` | unique strings exactly matching the released FASTA |
| `vectors` | `(3937, d)` | finite numeric matrix with fixed `d > 0` |

Row order is arbitrary because IDs control alignment. Missing, duplicate or
extra IDs are rejected. Finite zero-norm vectors use the frozen total
safe-cosine value `0.0`.

Registered paper cases may instead supply safetensors plus an
`embedding_index.tsv` containing unique `domain_id` and `embedding_key` values.
Both containers feed the same numeric method.

## 4. Run the official scorer

Prepare the scoring environment from
[`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md), then run:

```bash
python constraints/protein_sequence/PS02B/scoring/score_ps02b.py \
  --representations work/PS02B/my_model/representations.npz \
  --model-id my_model \
  --output-dir work/PS02B/my_model/score
```

For safetensors:

```bash
python constraints/protein_sequence/PS02B/scoring/score_ps02b.py \
  --representations work/PS02B/my_model/embeddings.safetensors \
  --representation-index work/PS02B/my_model/embedding_index.tsv \
  --model-id my_model \
  --output-dir work/PS02B/my_model/score
```

The formal entry fixes `delta_num=0`, PCG64, 20,000 draws and seed `20261006`.
Users provide model outputs, not scoring decisions.

## 5. Outputs

| File | Contents |
|---|---|
| `RAW_MEASUREMENTS.tsv.gz` | relation ID and anchor-positive/anchor-negative similarities |
| `INDIVIDUAL_SCORES.tsv.gz` | relation/object IDs, raw values, difference, credit and dependency identities |
| `SCORE_SUMMARY.json` | formal F, F 95% CI, coverage, frozen parameters and input identities |
| `OUTPUT_MANIFEST.json` | bytes and SHA256 for each public output |

Relation F credit is `1`, `0` or `-1`. Aggregate F follows the frozen hierarchy.
The CI resamples the 351 frozen groups. A formal output is
comparable under this release but is not automatically a paper-registered
result.

## 6. Comparator and paper cases

The complete NW–BLOSUM62 reference method, execution chain and canonical
outputs belong to
[`comparators/nw-blosum62-primary-comparator/`](comparators/nw-blosum62-primary-comparator/).
The five paper implementations and model results belong to
`paper_cases/protein_sequence/PS02B/`. This Constraint has no dependency on
those cases.
