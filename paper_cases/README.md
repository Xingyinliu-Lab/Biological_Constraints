
# Paper implementation cases

`paper_cases/` contains the actual Foundation Model implementations and evaluation records reported by the paper. It does not define a Biological Constraint and contains no independent scoring algorithm.

Each case records the model resource, the released Forward implementation or execution manual, the handoff to the unique public Constraint scorer and the current registered result. Large Forward representations need not be distributed; where they are absent, full reproduction regenerates them from the recorded checkpoint and case-specific Forward method. Comparator methods and results remain exclusively under `constraints/<domain>/<constraint_id>/comparators/`.

The current repository-level index covers all 110 model × Constraint results.

## Current entries

- [`CASE_INDEX.tsv`](CASE_INDEX.tsv) lists every registered paper case, its model, current result path, formal `F` and 95% confidence interval.
- Domain and Constraint directories provide the corresponding model Forward and scoring handoffs.
