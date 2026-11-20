# ST01 — low-sequence-identity structural similarity

ST01 measures whether a structure-model representation ranks a frozen structurally similar target above a matched structurally dissimilar alternative when both pairs satisfy the low-sequence-identity policy. The released Constraint is a self-contained measurement package: the user performs model Forward; ST01 validates the representation and performs all relation scoring, hierarchy aggregation and uncertainty calculation.

## Formal release

- Numeric Method: [`scoring/ST01_NUMERIC_METHOD_v0.1.1_FROZEN_20261010.json`](scoring/ST01_NUMERIC_METHOD_v0.1.1_FROZEN_20261010.json)
- Input binding: [`scoring/ST01_NUMERIC_INPUT_BINDING.json`](scoring/ST01_NUMERIC_INPUT_BINDING.json)
- Formal population: 28 coordinate objects, 16 bidirectional positive edges and 32 directed contrasts.

## Contract 1 — RefSet to model Forward

Encode every object in [`refset/objects.tsv`](refset/objects.tsv). `object_id` is the stable public identifier; `coordinate_locator` identifies its PDB under [`refset/coordinates/`](refset/coordinates/), and `coordinate_sha256` verifies the released coordinate bytes. All 28 coordinates listed in the frozen coordinate manifest are included and must be encoded without changing the formal object set.

The user chooses the model, checkpoint and representation layer. Model-specific preprocessing must preserve the mapping from each global structure representation to the original `object_id`.

## Contract 2 — model Forward to ST01 scorer

Save a NumPy NPZ containing exactly:

| Array | Shape | Requirement |
|---|---|---|
| `object_ids` | `(28,)` | Unique strings covering every object referenced by `relations.tsv`; order is arbitrary when it matches the embedding rows. |
| `embeddings` | `(28, d)` | Finite numeric global structure representations, `d > 0`, with no zero-norm row. |

The scorer also accepts `--pair-scores` for built-in or independently implemented reference methods. That TSV must contain `anchor_object_id`, `candidate_object_id` and finite numeric `score` for every pair required by the frozen relations.

## Run the official scorer

Prepare the standalone environment described in [`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md), then run:

```bash
python constraints/protein_structure/ST01/scoring/score_st01_direct.py \
  --embeddings work/ST01/my_model/embeddings.npz \
  --output-dir work/ST01/my_model/score
```

The bundled relations, Numeric Method and input binding are selected automatically. Advanced users may point the corresponding CLI options at another compatible release; the output records the assets actually used.

## Outputs

| File | Contents |
|---|---|
| `INDIVIDUAL_SCORES.tsv` | One row per directed contrast: IDs, true/negative cosine scores, margin and tie-aware `F_credit`. |
| `ENDPOINT_SCORES.tsv` | Positive-edge endpoint F before anchor aggregation. |
| `ANCHOR_SCORES.tsv` | Anchor-level F within dependency clusters. |
| `CLUSTER_SCORES.tsv` | Source-dependency-cluster F entering the macro mean and bootstrap. |
| `SCORE_SUMMARY.json` | F, F 95% CI, counts and input identities. |
| `F_BOOTSTRAP.npy` | Frozen PCG64 bootstrap draws on the F scale. |
| `OUTPUT_MANIFEST.json` | Byte size and SHA256 for each output file. |

ST01 uses the frozen hierarchy `contrast → endpoint → anchor → source dependency cluster → unweighted cluster macro F`. A higher F indicates more frequent preservation of the frozen structural ordering. Use of the released protocol does not automatically register a third-party result as a paper result.

## Built-in reference and paper cases

The [`geometry-comparator`](comparators/geometry-comparator/README.md) is the frozen deterministic Level-0 geometry reference and includes its method, executable, ID-level scores and formal result. Paper-model Forward implementations and their calls to this same scorer are recorded separately under [`paper_cases/protein_structure/ST01`](../../../paper_cases/protein_structure/ST01/README.md).
