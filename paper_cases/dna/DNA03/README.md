# DNA03 paper cases

These five cases record the paper models' actual Forward implementations and their results from the public DNA03 v0.2.6 F-only scorer. They do not contain a second scoring algorithm.

All cases encode the formal 1,400 views, provide primary and repeat core-local representations, and call `constraints/dna/DNA03/scoring/score_dna03.py`. Model environments and checkpoint identities are shared through `paper_cases/model_resources/dna/`.

| Model | Registered result |
|---|---|
| `caduceus-ps` | [`FORMAL_RESULT.json`](models/caduceus-ps/output/DNA03_v0.2.6_f_only/FORMAL_RESULT.json) |
| `crossdna-8m` | [`FORMAL_RESULT.json`](models/crossdna-8m/output/DNA03_v0.2.6_f_only/FORMAL_RESULT.json) |
| `dnabert-2` | [`FORMAL_RESULT.json`](models/dnabert-2/output/DNA03_v0.2.6_f_only/FORMAL_RESULT.json) |
| `janusdna` | [`FORMAL_RESULT.json`](models/janusdna/output/DNA03_v0.2.6_f_only/FORMAL_RESULT.json) |
| `nt-v2` | [`FORMAL_RESULT.json`](models/nt-v2/output/DNA03_v0.2.6_f_only/FORMAL_RESULT.json) |
