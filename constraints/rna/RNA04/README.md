# RNA04 — Codon consequence readability

RNA04 measures whether a representation makes amino-acid-changing single-nucleotide variants more distinguishable from synonymous alternatives at the same coding site. The whole-CDS v0.2.2 F-only release uses one public scoring pipeline for every model.

## 1. Frozen RefSet and Forward input

Use [`refset/RNA04_unique_endpoints.fasta`](refset/RNA04_unique_endpoints.fasta). It contains the complete 11,935-sequence evaluation universe required by [`refset/COMPLETE_SITE_ALT_LEDGER.tsv`](refset/COMPLETE_SITE_ALT_LEDGER.tsv): 11,506 variant rows, 3,920 sites, 560 parents and 554 dependency clusters.

- Each FASTA identifier is a sequence SHA256 and is the stable model-input ID.
- Encode every record exactly once; do not filter by sequence length, tokenizer behavior or model output.
- Preserve the full coding-oriented sequence. Windowing and pooling belong to the model Forward and must yield one vector per sequence ID.

## 2. Forward output contract

Save one NPZ containing exactly:

| Key | Shape | Requirement |
|---|---|---|
| `sequence_sha256` | `(11935,)` | unique strings; exact RefSet coverage; order arbitrary |
| `vectors` | `(11935, d)` | finite numeric matrix; `d > 0`; row-aligned with IDs |

The scorer rejects missing, duplicate or extra IDs, complex/boolean/non-numeric arrays, non-finite values and invalid dimensions.

## 3. Run the public scorer

Install the small scoring environment in [`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md), then run from the repository root:

```bash
python constraints/rna/RNA04/scoring/score_rna04.py \
  --vectors work/RNA04/my_model/EMBEDDINGS.npz \
  --model-id my_model \
  --output work/RNA04/my_model/score
```

The bundled RefSet is selected automatically. A compatible alternative population must provide both `--ledger` and `--fasta`. Before scoring, the public program validates FASTA sequence identities, complete coding-RNA syntax, exact WT/ALT coverage, the declared single-nucleotide change and position, codons, translated consequences, unique site/ALT endpoints, site identities, and consistent parent-to-WT-sequence and parent-to-dependency-cluster mappings. It records the actual ledger, FASTA and semantic population identities; `registered_reference_population` states whether those semantics match the released population without imposing a byte-for-byte file gate.

The public scorer performs all RNA04-specific work after Forward: ALT-minus-WT feature construction, frozen five-fold outer cross-fit, four-fold inner regularization selection, hierarchical weighting, site/parent/cluster aggregation and the fixed 20,000-draw cluster bootstrap.

## 4. Outputs

| File | Contents |
|---|---|
| `HELDOUT_PREDICTIONS.tsv` | every variant row, object IDs and held-out missense probability |
| `PAIR_CREDIT_LEDGER.tsv` | every missense-versus-synonymous pair, raw scores, margin, tie state and `F_credit` |
| `INDIVIDUAL_SCORES.tsv` | every formal site ID, class counts, comparison count and `individual_F` |
| `PARENT_SCORE_LEDGER.tsv` | parent-level F values |
| `CLUSTER_SCORE_LEDGER.tsv` | dependency-cluster F values used for F and its interval |
| `FORMAL_RESULT.json` | F, F 95% CI, coverage, selected C values and input identities |
| `BOOTSTRAP_MEANS.npy` | fixed bootstrap distribution |
| `OUTPUT_MANIFEST.json` | output sizes and SHA256 values |

The formal output is F throughout: each missense-versus-synonymous comparison contributes `1`, `0` or `-1` for win, tie or loss, and the stated hierarchy averages these F values. `formal_result=true` means that the frozen RNA04 protocol was used, not that a third-party run is part of the paper registry.

## 5. Comparator and paper cases

- The built-in raw-CDS genetic-code reference, code and formal results are in [`comparators/raw-cds-genetic-code-comparator-current/`](comparators/raw-cds-genetic-code-comparator-current/).
- Registered paper-model Forward implementations and outputs are under `paper_cases/rna/RNA04/`. They call this same public scorer and do not define another scoring method.
- The executable formal definition is [`scoring/RNA04_NUMERIC_METHOD_v0.2.2_F_ONLY.json`](scoring/RNA04_NUMERIC_METHOD_v0.2.2_F_ONLY.json).
