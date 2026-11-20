# PS04 public scoring package

`score_ps04.py` is the only user-facing entry. It performs the full chain:

1. validate complete NPZ or indexed safetensors representations;
2. generate every frozen leaf-level positive and negative measurement;
3. assign relation credits at frozen `delta_num=0`;
4. apply the fixed endpoint/HSF/anchor/source/topology/architecture hierarchy;
5. run the frozen 20,000-draw source-unit bootstrap and reuse sensitivity;
6. save ID-level and aggregate formal outputs.

Representations must be finite real numeric vectors. Complex arrays are
rejected before conversion, so their imaginary components can never be
silently discarded.

The representation measurement uses the historically executed PS04 operator:
stored vectors are cast to float64, normalized with float64 L2 norms, and
compared by float64 dot product. A zero-norm vector maps to the all-zero vector.
The earlier float32 safe-cosine artifact was a repeatability preflight, not the
formal leaf-measurement operator.

Internal modules are not alternative scorers:

| File | Role |
|---|---|
| `representations_to_leaf.py` | model representations to raw leaf measurements |
| `ps04_scoring_core.py` | frozen relation credit and hierarchy |
| `ps04_anchor_source_unit_bootstrap.py` | frozen uncertainty calculation |
| `ps04_formal_population_gate.py` | complete-population validation |
| `score_ps04_formal.py` | leaf measurements to ID-level and aggregate outputs |
| `score_ps04.py` | single public orchestration entry |

`FORMAL_BINDING.json` fixes the released scoring parameters. Admission is based
on complete semantic-population validation, not byte-for-byte identity to a
specific local copy of the RefSet.
