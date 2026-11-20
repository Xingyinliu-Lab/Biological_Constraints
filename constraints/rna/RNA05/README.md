# RNA05 — paired-state readability

RNA05 measures whether nucleotide-addressable representations distinguish paired from explicitly unpaired positions within frozen matched RNA blocks. The released population contains every complete 8-bin block; no runtime subset selection is allowed.

## 1. Frozen RefSet and Forward input

Use [`refset/RNA05_06_unique_inputs.fasta`](refset/RNA05_06_unique_inputs.fasta). It contains the 105 unique RNA sequences used by the 6,206 endpoint rows in [`refset/FORMAL_8BIN_ENDPOINT_LEDGER.tsv`](refset/FORMAL_8BIN_ENDPOINT_LEDGER.tsv).

- FASTA identifiers are sequence SHA256 values and are the stable model-input IDs.
- Preserve nucleotide coordinates exactly: `position0` is zero-based in the released sequence.
- Encode all required positions from all 105 sequences. Do not filter objects or positions using model behavior.

The formal population is 1,514 blocks, 103 parents and 20 dependency clusters. [`refset/RNA05_PRIMARY_MATCHED_BLOCKS_8BINS.tsv`](refset/RNA05_PRIMARY_MATCHED_BLOCKS_8BINS.tsv) records the block-level identities.

## 2. Forward output contract

The preferred public format is one sparse NPZ containing exactly:

| Key | Shape | Requirement |
|---|---|---|
| `sequence_sha256` | `(6206,)` | sequence IDs for formal endpoint rows |
| `position0` | `(6206,)` | integer zero-based coordinates |
| `vectors` | `(6206, d)` | finite numeric matrix; `d > 0` |

Each `(sequence_sha256, position0)` pair must be unique and the set must exactly match the formal endpoint ledger. Row order is arbitrary.

### Position mapping evidence

`position0` is a biological nucleotide coordinate, not a model token index. Users must map each released nucleotide position to their own tokenizer output and retain the following Forward provenance with the submitted representations:

- model, checkpoint and tokenizer identities;
- the token-to-nucleotide span rule, including k-mer or overlapping-token handling;
- BOS/EOS and other special-token handling;
- windowing, overlap reconciliation and restoration of original sequence coordinates;
- the pooling or readout rule used when one nucleotide maps to multiple tokens, or one token spans multiple nucleotides.

The mapping must produce one vector for every frozen `(sequence_sha256, position0)` target. The public scorer validates biological coordinates, exact target coverage and numerical inputs; it does not load or independently verify a third-party tokenizer. Mapping evidence is Forward provenance and is not an additional scorer input.

For registered historical paper cases, the scorer also accepts the executed sequence-keyed NPZ format: each NPZ key is one required sequence SHA256 and its value is a finite `(sequence_length, d)` matrix. This is a compatibility input to the same scorer, not a second scoring method.

## 3. Run the public scorer

Install [`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md), then run:

```bash
python constraints/rna/RNA05/scoring/score_rna05.py \
  --vectors work/RNA05/my_model/REPRESENTATIONS.npz \
  --model-id my_model \
  --output work/RNA05/my_model/score
```

The scorer automatically loads the formal endpoint ledger, block table and FASTA. Before fitting, it validates sequence hashes and alphabet, exact nucleotide coordinates and identities, unique endpoint coverage, paired-state labels, complete 8-bin membership, block counts, and parent/dependency-cluster consistency. Compatible alternative assets may be supplied with `--ledger`, `--blocks` and `--fasta`; their semantic population identity is recorded without imposing a byte-for-byte file gate. It then performs frozen leave-one-cluster-out fitting, four-fold inner lambda selection, hierarchical aggregation and the fixed 20,000-draw cluster bootstrap. Users do not choose scoring parameters.

## 4. Outputs

| File | Contents |
|---|---|
| `HELDOUT_PREDICTIONS.tsv` | all 6,206 endpoint IDs/coordinates and held-out paired probabilities |
| `INDIVIDUAL_SCORES.tsv` | 1,514 block IDs, raw paired/unpaired means and `F_credit` |
| `PARENT_SCORE_LEDGER.tsv` | 103 parent-level `parent_F` scores |
| `CLUSTER_SCORE_LEDGER.tsv` | 20 `cluster_F` scores used for F and its CI |
| `FORMAL_RESULT.json` | F, F 95% CI, coverage, selected lambdas and input identities |
| `BOOTSTRAP_MEANS.npy` | frozen F bootstrap distribution |
| `OUTPUT_MANIFEST.json` | output sizes and SHA256 values |

F is the sole formal reporting score. `formal_result=true` identifies use of the released protocol; it does not register a third-party result as a paper result.

## 5. Authority, comparator and paper cases

- Formal method: [`scoring/RNA05_NUMERIC_METHOD_v0.1.3_F_ONLY.json`](scoring/RNA05_NUMERIC_METHOD_v0.1.3_F_ONLY.json).
- The ViennaRNA pairing-propensity reference and its explicitly compatible evaluation scope are under [`comparators/viennarna-pairing-propensity/`](comparators/viennarna-pairing-propensity/).
- Registered paper-model Forward implementations and results are under `paper_cases/rna/RNA05/`; they call this same public scorer.
