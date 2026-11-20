# ST04

Use this page to evaluate a model against the ST04 constraint. The commands below connect the released RefSet, representation interface and public F-only scorer.

## Public scoring interface

- Official public entry: `constraints/protein_structure/ST04/scoring/score_st04.py`
- Scoring environment: [`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md)
- All compatible representations use the same frozen RefSet, numeric method, aggregation and uncertainty procedure.

## Model representation requirements

Provide an NPZ containing unique `residue_ids` and a finite two-dimensional `embeddings` matrix of shape `(n, d)`, where `d > 0`. Each ID must identify the exact structure/residue endpoint expected by the selected relation table; the file must provide complete coverage of that table and preserve ID-to-row alignment.

Residue IDs use the released structure-object and canonical polymer position,
for example `1byiA00:71`. The scorer is defined for residue-addressable
structure representations that can return every anchor, true-target and
negative endpoint in the relation ledger. A whole-structure pooled vector does
not satisfy this interface.

The Forward population is the 586 admitted structures in `refset/objects.tsv`.
Use each row's released `coordinate_locator`; map model residues back to the
canonical polymer positions used by the relation ledger, and emit one vector
for every required residue endpoint. Users do not construct the 152,180
relations or fit the nested learner themselves.

Document the Forward mapping with every submission: model/checkpoint and
tokenizer or structure-tokenizer identity, the residue-to-model-index mapping,
special-token or padding treatment, any windowing or cropping rule, and the
residue-level readout used for each emitted vector. This documentation is supplied
by the model user; the public scorer validates the submitted endpoint IDs and
numeric representation, not a third-party model's internal tokenization.

## Run the scorer

### Run the official scorer

Provide one finite residue-addressable vector for every frozen endpoint in an
NPZ containing `residue_ids` and `embeddings`. The public scorer retains the
frozen nested cross-fit, lambda selection, hierarchy, and bootstrap, but does
not require a paper-model execution binding:

```bash
python constraints/protein_structure/ST04/scoring/score_st04.py \
  --relations constraints/protein_structure/ST04/refset/relations.tsv \
  --refset-release constraints/protein_structure/ST04/refset/ST04_ASSAY_REFSET_RELEASE_v0.1.1_FROZEN_20261012.json \
  --numeric-method constraints/protein_structure/ST04/scoring/ST04_NUMERIC_METHOD_v0.1.4_F_ONLY.json \
  --split-release constraints/protein_structure/ST04/refset/ST04_SPLIT_RELEASE_v0.1.1_FROZEN_20261012.json \
  --embeddings my_model_residue_embeddings.npz \
  --model-id my_model \
  --out-dir work/ST04/my_model
```

The primary output is `EVALUATION_RESULT.json`, containing `F`, `F_95_CI`, population counts, bootstrap identity
and the SHA256 identities of all scored input assets. Selected lambda values are
recorded in `LAMBDA_SELECTION.json`. The directory also contains release QC,
coverage, hierarchy-score and held-out-prediction evidence.
`OOF_QUERY_SCORES.jsonl` is the formal ID-level output: each query record
contains its held-out positive score, negative scores, per-relation credits and
query credit. These values feed the frozen hierarchy and aggregate statistics.

The command above uses all 152,180 relations in the released RefSet. For the paper ESM-IF1 and ESM3-STRUCT interfaces, a documented structural-input limitation excludes seven relations: their 152,173-relation table is `paper_cases/protein_structure/ST04/scope/ESM_INPUT_LIMITED_RELATIONS.tsv`. These cases pass that table as `--relations` together with `--scope-record paper_cases/protein_structure/ST04/scope/ESM_INPUT_LIMITATION.json`; the scorer checks the scope record, relation-table hash and population counts. The other three paper models use all 152,180 relations.

## Formal method files

Frozen scoring implementation:

- `constraints/protein_structure/ST04/scoring/st04_numeric_core_v0_1_1.py`
- `constraints/protein_structure/ST04/scoring/st04_io_qc.py`
- `constraints/protein_structure/ST04/scoring/ST04_NUMERIC_METHOD_v0.1.4_F_ONLY.json`

## Released RefSet

Ready-to-use RefSet assets:

- `constraints/protein_structure/ST04/refset/dependency.tsv`
- `constraints/protein_structure/ST04/refset/FORMAL_COORDINATE_INPUT_MANIFEST.json`
- `constraints/protein_structure/ST04/refset/manifest.json`
- `constraints/protein_structure/ST04/refset/objects.tsv`
- `constraints/protein_structure/ST04/refset/objects/experimental_pdb/` (586 PDB files)
- `constraints/protein_structure/ST04/refset/PUBLIC_REFSET_VALIDATION.json`
- `constraints/protein_structure/ST04/refset/residue_evidence.tsv` (98,209 frozen residue mappings)
- `constraints/protein_structure/ST04/refset/ST04_ASSAY_REFSET_RELEASE_v0.1.1_FROZEN_20261012.json`
- `constraints/protein_structure/ST04/refset/ST04_SPLIT_RELEASE_v0.1.1_FROZEN_20261012.json`
- `constraints/protein_structure/ST04/refset/VALIDATION.json`
