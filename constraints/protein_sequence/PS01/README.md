# PS01 — protein family/clan relationship readability

PS01 measures whether whole-protein representations rank a same-family positive above a different-family negative drawn from the same Pfam clan. The only released population is the final Broad16K/L release; historical small populations are not formal outputs.

## 1. RefSet to model Forward input

Encode all 2,316 sequences in [`refset/FORMAL_SEQUENCES.faa.gz`](refset/FORMAL_SEQUENCES.faa.gz). FASTA record IDs are the stable object IDs used by the formal ledgers. Do not filter sequences because of model length or tokenizer behavior.

The frozen population contains 1,256 contrasts across 124 clans, 420 anchor families and 941 anchors. Object metadata and formal triples are in `FORMAL_OBJECT_LEDGER.tsv.gz` and `FORMAL_CONTRAST_LEDGER.tsv.gz`.

## 2. Model Forward output contract

Submit one NPZ containing exactly:

| Key | Shape | Requirement |
|---|---|---|
| `object_ids` | `(2316,)` | unique released FASTA IDs |
| `vectors` | `(2316, d)` | finite whole-protein vectors with `d > 0` and nonzero row norm |

Row order is arbitrary when IDs and vectors remain aligned. Pooling/tokenization belongs to the model Forward and must be documented by the user; PS01 does not accept precomputed contrast scores instead of representations.

The registered paper cases used `embeddings.safetensors` plus `embedding_index.tsv`. The same scorer accepts that executed format with `--representation-index`; it is a compatibility input to the same calculation, not a second scorer.

## 3. Run the public scorer

Install [`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md), then run:

```bash
python constraints/protein_sequence/PS01/scoring/score_ps01.py \
  --representations work/PS01/my_model/representations.npz \
  --model-id my_model \
  --output-dir work/PS01/my_model/score
```

The scorer computes both cosine similarities, tie-aware contrast credit, the frozen hierarchy and the fixed 20,000-draw clan bootstrap. Users do not choose formal scoring parameters.

## 4. Outputs

| File | Contents |
|---|---|
| `INDIVIDUAL_SCORES.tsv.gz` | 1,256 contrast IDs/objects, raw positive and negative cosine values, and formal `F_credit` |
| `F_BOOTSTRAP.npy` | 20,000 bootstrap draws on the formal F scale |
| `EVALUATION_RESULT.json` | F, F 95% CI, win/tie/loss counts, clan F values, coverage and input identities |
| `OUTPUT_MANIFEST.json` | output byte sizes and SHA256 values |

F is the formal score throughout the released pipeline. A win, tie and loss contribute `1`, `0` and `-1`, respectively; the frozen nested hierarchy and clan bootstrap operate directly on that F scale.

## 5. Comparator and paper cases

The NW+BLOSUM62 global-alignment reference, executable chain and registered result are under [`comparators/nw-blosum62-global/`](comparators/nw-blosum62-global/). Five paper-model cases are stored under `paper_cases/protein_sequence/PS01/`; they use this same public scorer.
