# PS05 public scoring package

`score_ps05.py` is the single public entry. It validates complete model output,
computes all frozen leaf measurements, writes ID-level scores, applies the
released hierarchy, runs the 20,000-draw nested bootstrap and reports reuse
sensitivity plus aggregate F/F 95% CI.

Both NPZ and indexed safetensors inputs must contain finite real numeric
vectors. Complex arrays are rejected before conversion; imaginary components
are never silently discarded.

Internal modules are consecutive stages, not alternative scorers:

| File | Role |
|---|---|
| `representations_to_leaf.py` | representations to frozen leaf measurements |
| `score_ps05_formal.py` | validation, relation credits, hierarchy and outputs |
| `bootstrap_ps05_hsf_nested.py` | frozen nested bootstrap |
| `score_ps05.py` | single public orchestration entry |

Users cannot override formal scoring decisions. Model/checkpoint provenance is
recorded by paper cases or by the external user and does not change PS05 math.

The scorer validates submitted leaf identities against the supplied planned
population and applies the frozen PS05 contracts and numeric rules. It records
both the actual relation-ledger SHA256 and the registered release SHA256. A
byte mismatch is provenance, not a universal admission gate for an otherwise
schema-valid, internally consistent compatible RefSet.
