# ST02

Use this page to evaluate a model against the ST02 constraint. The commands below connect the released RefSet, representation interface and public F-only scorer.

## Public scoring interface

- Official public entry: `constraints/protein_structure/ST02/scoring/score_st02.py`
- Scoring environment: [`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md)
- All compatible representations use the same frozen RefSet, numeric method, aggregation and uncertainty procedure.

## Model representation requirements

Provide an NPZ containing exactly unique `object_ids` and a finite two-dimensional `embeddings` matrix of shape `(n, d)`, where `d > 0`. IDs must exactly cover every object used by the frozen ST02 relation table; row order is arbitrary when IDs and rows remain aligned.

The Forward population is the 509 rows in `refset/objects.tsv`. Each row gives
the stable `record_id`, the released PDB file in `coordinate_locator`, its
SHA256, and the structural classification metadata. Run the model on every
listed structure and preserve `record_id` as the corresponding `object_ids`
value. The public scorer performs cosine contrasts, hierarchy aggregation and
uncertainty calculation; users do not construct relation ledgers themselves.

For position- or residue-aware models, document the model/checkpoint,
the mapping from model residue or token indices to frozen
canonical residue IDs, special-token or padding handling, windowing/overlap,
and the final readout. This evidence documents model-specific extraction; it
is not an additional scorer input.

```python
import numpy as np
np.savez_compressed(
    "work/ST02/my_model/embeddings.npz",
    object_ids=np.asarray(object_ids, dtype="U128"),
    embeddings=np.asarray(vectors, dtype=np.float32),
)
```

## Run the official scorer

```bash
python constraints/protein_structure/ST02/scoring/score_st02.py \
  --refset-dir constraints/protein_structure/ST02/refset \
  --embeddings work/ST02/my_model/embeddings.npz \
  --model-id my_model \
  --output-dir work/ST02/my_model/score
```

This entrypoint validates the complete object universe and uses the frozen ST02 hierarchy and 20,000-draw dependency-cluster bootstrap. It does not require a paper-model execution binding. Its output records the supplied model identity and is not automatically added to the registered paper results.

`score/SCORE_SUMMARY.json` reports final `F`, `F_95_CI`, topology and
source-superfamily/topology stratified F scores with their intervals, aggregation
counts, contrast/dependency-cluster coverage, and input identities.
`INDIVIDUAL_SCORES.tsv` preserves all 4,584 contrast IDs, positive and negative
cosines, their raw difference and the formal relation credit. Endpoint, target
superfamily, anchor and cluster ledgers plus the frozen bootstrap draws and an
output manifest are written beside it.

## Formal method files

Frozen scoring implementation:

- `constraints/protein_structure/ST02/scoring/st02_numeric_core.py`
- `constraints/protein_structure/ST02/scoring/ST02_NUMERIC_METHOD_v0.1.3_F_ONLY.json`

## Released RefSet

Ready-to-use RefSet assets:

- `constraints/protein_structure/ST02/refset/objects.tsv`
- `constraints/protein_structure/ST02/refset/coordinates/` (509 PDB files)
- `constraints/protein_structure/ST02/refset/relations.tsv`
