# PS08: progressive protein-edit adjacency

PS08 evaluates whether a model representation preserves the local ordering of a frozen two-step protein-edit ladder:

```text
x0 → x1 → x2
```

For both endpoints, the one-edit neighbor `x1` should be more similar than the two-edit endpoint. This directory is the PS08 formal public evaluation unit: it contains the frozen scientific authority, ready-to-use RefSet, representation contract and model-independent scorer.

## 1. Scientific definition

PS08 asks whether the one-edit neighbor is represented as more similar than the two-edit endpoint from both ends of the frozen path. A high score supports local edit-adjacency readability under the declared representation; it does not by itself establish mutation fitness, functional consequence, structural accuracy or evolutionary plausibility.

The two substitutions affect distinct amino-acid positions and the frozen sequence relations satisfy:

```text
H(x0,x1) = 1,  H(x1,x2) = 1,  H(x0,x2) = 2
```

The fixed numerical rules are:

- cosine similarity on one vector per frozen protein object;
- two endpoint comparisons per ladder;
- tie-aware credit with frozen numerical tolerance `delta_num = 0`;
- aggregation from leaf to first edit, source parent and source unit;
- 20,000 source-unit bootstrap draws with seed `20261002`;
- `F = 2R - 1`.

## 2. Frozen authority and ready-to-use RefSet

The released RefSet and public scoring implementation define the PS08 measurement. The public scoring path does not require a paper-model checkpoint, RunPlan or Forward manifest. The model-independent input contract is defined in Section 3 and enforced by `scoring/representations_to_leaf.py`.

Use these released files directly:

| File | Purpose |
|---|---|
| [`refset/ENDPOINT_SEQUENCES.faa.gz`](refset/ENDPOINT_SEQUENCES.faa.gz) | All 12,288 formal endpoint protein sequences; FASTA IDs are the model-output object IDs |
| [`refset/PLANNED_POPULATION_MANIFEST.json`](refset/PLANNED_POPULATION_MANIFEST.json) | The 4,096 frozen ladders and their x0/x1/x2, parent and dependency identities |
| [`refset/POPULATION_STATUS_SUPERSESSION.json`](refset/POPULATION_STATUS_SUPERSESSION.json) | Authority record promoting the byte-identical population to formal release |

Formal population identity:

```text
release: PS08_ASSAY_REFSET_RELEASE_v1.0.0_20261008
population identity: 2ee0bbde63d2ed1bcec5b0c12c00b8a10385cc98523a9aa8966ae1a4fd5a9fb3
endpoint FASTA SHA256: d5d8f39575f5b071278a223875570e32e9b4e3a321b5abba416e28f074445ea4
```

The historical filename `PLANNED_POPULATION_MANIFEST.json` is retained for byte identity. Its embedded development status is superseded by `POPULATION_STATUS_SUPERSESSION.json`. Results claimed against the official PS08 release must use this population. The scorer may process another schema-compatible population for research use, but records its actual semantic and byte identities so it cannot be confused with the registered release.

## 3. Representation input protocol

Generate one fixed-dimensional vector for every FASTA record and save it in a compressed NumPy archive. The scorer accepts representations from pretrained models, trained learners or probes; these are different experimental designs and must be reported and interpreted accordingly. A result presented as direct pretrained-representation readability must not use ladder roles, edit labels or evaluation feedback during representation construction.

```text
representations.npz
├── object_ids     string array, shape (12288,)
└── embeddings     real numeric array, shape (12288, d)
```

Requirements:

- every formal object ID appears exactly once;
- no additional object IDs are allowed;
- `d` is positive and identical for every object;
- all values are finite real numbers; complex and Boolean arrays are rejected;
- finite zero-norm vectors are permitted and use the frozen safe-cosine value `0` against every vector;
- the same checkpoint, layer, token/residue mapping, long-sequence rule and pooling rule are used for every object;
- row order is arbitrary when `object_ids` and `embeddings` remain aligned.

Minimal writer:

```python
import numpy as np
from pathlib import Path

output = Path("work/PS08/my_model/representations.npz")
output.parent.mkdir(parents=True, exist_ok=True)
np.savez_compressed(
    output,
    object_ids=np.asarray(object_ids, dtype="U128"),
    embeddings=np.asarray(vectors, dtype=np.float32),
)
```

The repository does not prescribe how an external model produces its vectors. Model installation and inference belong to that model; PS08 fixes only the biological input population and representation-to-score interface.

## 4. Run the public scoring pipeline

Create the standalone scoring environment described in [`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md), then run the two native PS08 stages from the repository root.

```bash
python constraints/protein_sequence/PS08/scoring/representations_to_leaf.py \
  --population constraints/protein_sequence/PS08/refset/PLANNED_POPULATION_MANIFEST.json \
  --objects-fasta constraints/protein_sequence/PS08/refset/ENDPOINT_SEQUENCES.faa.gz \
  --representations work/PS08/my_model/representations.npz \
  --output work/PS08/my_model/FORMAL_LEAF_SCORES.tsv

python constraints/protein_sequence/PS08/scoring/score_ps08_formal.py \
  --leaf-ledger work/PS08/my_model/FORMAL_LEAF_SCORES.tsv \
  --population constraints/protein_sequence/PS08/refset/PLANNED_POPULATION_MANIFEST.json \
  --objects-fasta constraints/protein_sequence/PS08/refset/ENDPOINT_SEQUENCES.faa.gz \
  --outdir work/PS08/my_model/score
```

The programs perform the complete model-independent evaluation:

1. recompute the population semantic identity, require its declaration to match the actual members, and validate endpoint sequence hashes plus the Hamming 1/1/2 relation;
2. validate the NPZ schema, complete object universe and finite numeric values;
3. generate the runtime leaf-score ledger using frozen safe cosine;
4. validate every leaf identity against the same frozen population;
5. compute F and the frozen source-unit bootstrap interval.

Use a new work directory for each evaluation so files from different models are not mixed.

## 5. Results and interpretation

The canonical outputs are:

```text
work/PS08/my_model/score/INDIVIDUAL_SCORES.tsv
work/PS08/my_model/score/SCORE_SUMMARY.json
work/PS08/my_model/score/OUTPUT_MANIFEST.json
```

`INDIVIDUAL_SCORES.tsv` contains all 4,096 evaluation IDs, the four raw similarities, left/right differences, F credits, dependency identities, hierarchy F scores and each leaf's contribution to F.

`SCORE_SUMMARY.json` contains:

| Field | Meaning |
|---|---|
| `F` | Formal signed score in `[-1,1]`; `0` is the tie reference |
| `F_95_CI` | Frozen source-unit bootstrap 95% confidence interval for F |
| `counts` | Leaf, first-edit, source-parent and source-unit coverage |
| `seed`, `draws`, `bootstrap` | Frozen uncertainty calculation identity |
| `population_identity_sha256` | Declared population identity retained for compatibility |
| `population_validation` | Actual file SHA256, semantic identity and registered-release match status |
| `individual_scores` | ID-level output path, row count and SHA256 |

Runtime intermediates such as `FORMAL_LEAF_SCORES.tsv` and the bootstrap array are generated inside the work directory. They are not required public inputs.

For a reportable third-party result, retain a run note beside the score containing the model ID, checkpoint revision or weight SHA256, selected layer/readout, token/residue mapping, long-sequence rule, pooling rule, representation-file SHA256, PS08 release ID and the exact commands above.

### Validity requirements

A result is valid only when:

- all 12,288 formal objects were encoded exactly once;
- the population schema, object set, relation identities and dependency hierarchy pass validation;
- no NaN or Inf values occur; zero-norm vectors follow the frozen safe-cosine rule and are not silently removed;
- the same declared representation protocol was applied to every object;
- coverage is 4,096 leaves, first edits, source parents and source units;
- the scorer reports seed `20261002`, 20,000 draws and `source_unit with replacement`;
- `delta_num` remains `0`;
- the released converter and formal scorer are used without modifying `delta_num`, population membership or aggregation.

## 6. Comparator reference

The exact Hamming method under [`comparators/hamming-oracle/`](comparators/hamming-oracle/) is the frozen formal Level-0 white-box reference, deterministic ceiling and positive control. Truth construction and the white-box comparator intentionally coincide for this mechanical relation. The PS08 comparator requirement is CLOSED; no additional NW, BLAST, HMMER, HHsearch or other traditional comparator is required.
