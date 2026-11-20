# PS05 — within-topology structural-neighbour relation

PS05 is a frozen public measurement package. Users Forward the released 321
protein-domain sequences and provide one whole-domain representation per ID.
PS05 performs every subsequent relation measurement, frozen hierarchical
aggregation, nested bootstrap, reuse sensitivity and formal output step.

## 1. Frozen release

The active release is `PS05_ASSAY_REFSET_RELEASE_v1.0.1_20261008`. It contains
321 Forward objects and 1,628 frozen leaf relations. Each relation compares an
anchor-positive pair with a frozen matched negative inside the released
within-topology population. Relation membership and all biological hierarchy
fields are fixed by the RefSet and planned-population identity.

## 2. Contract A: RefSet to model Forward

Forward every record in
[`refset/PS05_ENDPOINT_SEQUENCES.faa`](refset/PS05_ENDPOINT_SEQUENCES.faa),
preserving the FASTA token as `object_id`. Do not rebuild or filter the
population using model outputs. Tokenization, context handling and model
pooling belong to Forward; the submitted result must be one consistent
whole-domain vector for each of the 321 objects.

## 3. Contract B: Forward output to scorer

Preferred NPZ input:

| Key | Shape | Requirement |
|---|---:|---|
| `object_ids` | `(321,)` | unique strings exactly matching the released FASTA |
| `vectors` | `(321, d)` | finite, real, nonzero whole-domain vectors with fixed `d > 0` |

Row order is arbitrary. Missing, duplicate or extra IDs, complex values,
NaN/Inf and zero-norm vectors are rejected. Registered paper cases may use safetensors plus an index
with unique `domain_id` and `embedding_key` columns.

The scorer records SHA256 values for representation provenance, but it does
not require a submitted representation or index to match a paper-case hash.
This allows the same frozen RefSet and scoring method to evaluate newly
generated representations and user-supplied models. Hash checks retained for
the released RefSet and scoring contracts protect the fixed measurement
assets rather than restrict the model being evaluated.

## 4. Run the official scorer

After preparing [`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md):

```bash
python constraints/protein_sequence/PS05/scoring/score_ps05.py \
  --representations work/PS05/my_model/representations.npz \
  --model-id my_model \
  --output-dir work/PS05/my_model/score
```

For indexed safetensors, also pass
`--representation-index work/PS05/my_model/embedding_index.tsv`.

The formal entry fixes `delta_num=0`, PCG64, seed `20261008`, 20,000 draws,
the released hierarchy and the frozen endpoint-reuse sensitivity rule.

## 5. Outputs

| File | Contents |
|---|---|
| `RAW_MEASUREMENTS.tsv` | all 1,628 relation identities and positive/negative similarities |
| `INDIVIDUAL_SCORES.tsv.gz` | stable evaluation ID, raw measurements, margin and relation credit |
| `AGGREGATION_TRACE.jsonl` | frozen hierarchy nodes, child values, weights and aggregate values |
| `SCORE_SUMMARY.json` | F, F 95% CI, coverage, hierarchy, reuse sensitivity and identities |
| `OUTPUT_MANIFEST.json` | byte size and SHA256 of public outputs |

The formal reporting score is F; wins, ties and losses contribute `1`, `0` and
`-1`, respectively. Formal protocol use does not itself
register a third-party result as a paper result.

## 6. Comparators and paper cases

Comparator methods, complete executions and canonical results belong to
[`comparators/`](comparators/). Paper model Forward implementations and results
belong to `paper_cases/protein_sequence/PS05/`. All call this same public
scoring package.
