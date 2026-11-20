# RNA10 | Exon Recognition

RNA10 measures whether whole-RNA representations separate single-nucleotide variants with disrupting exon-recognition effects from matched preserving variants under the released MFASS-derived reference population.

## Public inputs and formal method

`refset/` contains the released WT and mutant sequences, object maps and truth relations. `scoring/RNA10_NUMERIC_METHOD_v0.1.3_F_ONLY.json`, `scoring/RNA10_FORMAL_RUNPLAN.json` and `scoring/RNA10_FORMAL_SCORING_BINDING.json` define the one public evaluation protocol.

Encode every sequence in `refset/RNA10_UNIQUE_INPUTS.fasta` as one finite whole-sequence vector. The NPZ must contain exactly `sequence_sha256` and `vectors`; its IDs must cover the released WT and mutant sequence universe exactly.

## Score a new model

From the repository root:

```bash
python constraints/rna/RNA10/scoring/build_measurement_from_vectors.py \
  --truth-relations constraints/rna/RNA10/refset/truth_relations.jsonl \
  --objects-wt constraints/rna/RNA10/refset/objects_wt.jsonl \
  --objects-mutant constraints/rna/RNA10/refset/objects_mutant.jsonl \
  --vectors work/RNA10/<model>/REPRESENTATIONS.npz \
  --out work/RNA10/<model> \
  --model <model>

python constraints/rna/RNA10/scoring/score_rna10.py \
  --measurement-ledger work/RNA10/<model>/MEASUREMENT_LEDGER.json \
  --model-id <model> \
  --out work/RNA10/<model>/FORMAL_RESULT.json
```

The supplied scoring files are defaults and may be explicitly overridden only with compatible public replacements.

## F-only score

The frozen cross-fit protocol constructs `z_mutant - z_WT`, uses the fixed L2 logistic-ridge selection and leave-one-source-unit-out evaluation, then compares matched disrupting and preserving held-out scores. Each leaf receives `F_credit` of `1`, `0`, or `-1` for win, tie, or loss. Equal block and source-unit means produce F. `F_95_CI` is the 20,000-draw percentile interval over released uncertainty dependency units (seed 10010).

Raw held-out decision scores remain in the measurement ledger as measurement evidence; all public scoring outputs are F-only.

## Comparator

The deterministic positional k-mer ridge comparator is in [`comparators/raw-positional-kmer-weighted-ridge-comparator/`](comparators/raw-positional-kmer-weighted-ridge-comparator/).
