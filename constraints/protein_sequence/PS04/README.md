# PS04 — cross-superfamily topology relation

PS04 is a frozen public measurement package. A user performs model Forward on
the released protein-domain sequences and submits one whole-domain vector per
object. PS04 then computes all frozen relation measurements, hierarchical
credits, dependency-aware aggregation, sensitivity analysis and uncertainty.

Paper models and external models use the same public entry.

## 1. Frozen release

The active release is `PS04_ASSAY_REFSET_RELEASE_v1.0.1_L_4K_20261007`.
Its formal population contains 4,096 materialized objects, 3,641 endpoint
sequences used for Forward, 2,790 complete anchor bundles, 1,750,186 leaf
relations and 2,728 frozen bootstrap units. This is the only released large
population; historical development populations are not public alternatives.

Each leaf asks whether an anchor is closer to its designated positive than to
a frozen cross-superfamily negative under the hierarchy defined by the
Scientific Word, relation ledger and v0.1.2 numeric binding.

## 2. Contract A: RefSet to model Forward

Forward every record in
[`refset/FORMAL_ENDPOINTS.faa`](refset/FORMAL_ENDPOINTS.faa), preserving the
FASTA token as `object_id`. The file contains exactly 3,641 unique sequences.
Do not regenerate relations, substitute the 4,096-object materialization
ledger for the Forward FASTA, or filter objects according to model behavior.

The user controls the model, tokenizer, context handling and representation
extraction. The submitted representation must nevertheless be one consistent
whole-domain vector for each released FASTA object.

## 3. Contract B: Forward output to scorer

Preferred NPZ input:

| Key | Shape | Requirement |
|---|---:|---|
| `object_ids` | `(3641,)` | unique strings exactly matching `FORMAL_ENDPOINTS.faa` |
| `vectors` | `(3641, d)` | finite real numeric matrix with fixed `d > 0`; complex values are rejected |

Row order is arbitrary because IDs control alignment. Missing, duplicate or
extra objects are rejected. A finite zero-norm vector is accepted and produces
the total cosine value zero. Registered paper cases may instead use a
safetensors file plus an index containing unique `domain_id` and
`embedding_key` columns.

## 4. Run the official scorer

Prepare [`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md), then run:

```bash
python constraints/protein_sequence/PS04/scoring/score_ps04.py \
  --representations work/PS04/my_model/representations.npz \
  --model-id my_model \
  --output-dir work/PS04/my_model/score
```

For a registered safetensors Forward:

```bash
python constraints/protein_sequence/PS04/scoring/score_ps04.py \
  --representations work/PS04/my_model/embeddings.safetensors \
  --representation-index work/PS04/my_model/embedding_index.tsv \
  --model-id my_model \
  --output-dir work/PS04/my_model/score
```

The scorer fixes `delta_num=0`, PCG64, 20,000 draws, seed `20261007`, the
2,728-unit resampling hierarchy and the `anchor_bundle_reuse_count > 100`
sensitivity rule. Users provide model outputs, not scoring decisions.

The scorer validates the submitted RefSet schema, complete relation-to-leaf
mapping, object coverage and dependency structure. Actual input identities are
recorded in the formal output; byte identity alone is not an admission requirement.

## 5. Outputs and interpretation

| File | Contents |
|---|---|
| `RAW_MEASUREMENTS.tsv.gz` | every formal leaf identity and positive/negative similarity |
| `INDIVIDUAL_SCORES.tsv.gz` | stable evaluation ID, raw similarities, margin, F credit and bootstrap identity |
| `AGGREGATION_TRACE.jsonl.gz` | hierarchy node ID, child IDs, weights and aggregate F for tracing each formal contribution |
| `FORMAL_RESULT.json` | F, 95% CI, full coverage, hierarchy identity and reuse sensitivity |
| `OUTPUT_MANIFEST.json` | byte size and SHA256 of each public output |

Each relation receives F credit `1`, `0` or `-1` at frozen tolerance zero. The
formal F score is obtained through the frozen multi-level hierarchy.
`formal_result=true` means that the public formal
protocol was used. It does not add an external model to the paper registry.

## 6. Comparator and paper cases

The complete built-in comparators, their method records and canonical outputs
belong to [`comparators/`](comparators/README.md). Paper Forward implementations and
registered model results belong to `paper_cases/protein_sequence/PS04/`.
Neither is a dependency of this Constraint scorer.
