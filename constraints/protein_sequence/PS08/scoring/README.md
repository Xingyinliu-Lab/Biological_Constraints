# PS08 scoring

## Environment

Create the standalone scoring environment described in [`ENVIRONMENT.md`](ENVIRONMENT.md). It is intentionally separate from model CUDA environments.

## Programs

### `representations_to_leaf.py`

Convert model representations into the frozen PS08 leaf-score ledger.

```bash
python constraints/protein_sequence/PS08/scoring/representations_to_leaf.py \
  --population constraints/protein_sequence/PS08/refset/PLANNED_POPULATION_MANIFEST.json \
  --objects-fasta constraints/protein_sequence/PS08/refset/ENDPOINT_SEQUENCES.faa.gz \
  --representations work/PS08/my_model/representations.npz \
  --output work/PS08/my_model/FORMAL_LEAF_SCORES.tsv
```

### `score_ps08_formal.py`

Compute the formal PS08 F score and confidence interval from the generated leaf ledger.

```bash
python constraints/protein_sequence/PS08/scoring/score_ps08_formal.py \
  --leaf-ledger work/PS08/my_model/FORMAL_LEAF_SCORES.tsv \
  --population constraints/protein_sequence/PS08/refset/PLANNED_POPULATION_MANIFEST.json \
  --objects-fasta constraints/protein_sequence/PS08/refset/ENDPOINT_SEQUENCES.faa.gz \
  --outdir work/PS08/my_model/score
```

| Argument | Required | Type/default | Meaning |
|---|---:|---|---|
| `--leaf-ledger` | yes | TSV/TSV.GZ path | Output from `representations_to_leaf.py`; all 4,096 frozen leaf identities and four relation-score columns are required. |
| `--population` | yes | JSON path | The same frozen population manifest used by the converter. |
| `--objects-fasta` | yes | FASTA.GZ/FASTA path | Endpoint sequences used to verify object hashes and the Hamming 1/1/2 relation. |
| `--outdir` | yes | directory | New output directory for ID-level scores, aggregate statistics, the bootstrap array and their manifest. |
| `--delta-num` | no | float | Frozen formal value `0` is used automatically; any nonzero override is rejected. |

## Inputs

Use the RefSet paths and representation schema documented in [`../README.md`](../README.md). Inputs must cover the complete ID universe, contain finite real numeric values and preserve frozen dependency units. The declared population identity must equal the semantic identity recomputed from `canonical_members`; registered-release status is determined from that recomputed identity rather than from a self-declared value.

## Outputs

Write outputs to `work/PS08/<model>/` for a new run.

| Output | Contents |
|---|---|
| `FORMAL_LEAF_SCORES.tsv` | converter output containing the four raw cosine measurements per leaf |
| `FORMAL_LEAF_SCORES.tsv.metadata.json` | actual population, representation and generated-ledger identities |
| `score/INDIVIDUAL_SCORES.tsv` | stable leaf and hierarchy IDs, four raw similarities, two margins, F credits, hierarchy F scores and overall F contribution |
| `score/SCORE_SUMMARY.json` | F, F 95% CI, coverage and frozen parameters |
| `score/BOOTSTRAP_F.float64le.bin` | frozen source-unit bootstrap F draws |
| `score/OUTPUT_MANIFEST.json` | size and SHA256 of every canonical scorer output |

The individual-score table and summary are canonical public outputs, not disposable intermediates. The scorer validates population schema, identities, endpoint closure and hierarchy consistency, and records the actual population-file and semantic identities.

## One public scoring pipeline
