# RNA01-02 — whole-RNA edit-path adjacency

RNA01-02 tests whether a model's whole-sequence representation preserves the local topology of a two-step RNA edit path. Each formal ladder contains two one-edit edges and one two-edit diagonal. The released scorer checks both endpoints, aggregates by dependency-closed outer cluster, and reports the formal score `F` with a 95% cluster-bootstrap interval. Endpoint F credits are `1` (win), `0` (tie), and `-1` (loss).

The frozen **L16K large population** is the sole public and formal population.

## 1. RefSet to model Forward input

Use [`refset/RNA0102_L16K_EDIT_REFSET/OBJECT_LEDGER.json`](refset/RNA0102_L16K_EDIT_REFSET/OBJECT_LEDGER.json). It contains 49,152 formal RNA objects representing 49,149 unique sequence hashes. For each row, `object_id` is the stable evaluation ID, `sequence` is the exact A/C/G/U RNA input, and `sequence_sha256` is the scorer-facing row identity. Identical sequences shared by multiple formal objects are encoded once and joined by their common hash.

Encode every listed sequence exactly once. Do not change the sequence, object universe, edit relations, cluster membership, or large-population selection. Model loading, tokenization, windowing and pooling are the user's Forward responsibility; the final representation must be one whole-sequence vector per `sequence_sha256`.

The calibration panel is [`refset/RNA0102_calibration_ladders.json`](refset/RNA0102_calibration_ladders.json). Run this smaller panel twice through the same Forward path. These are repeat inferences, not repeat model training.

Materialize the two ready-to-Forward FASTA files from the released JSON assets:

```bash
python constraints/rna/RNA01-02/refset/materialize_forward_fastas.py \
  --object-ledger constraints/rna/RNA01-02/refset/RNA0102_L16K_EDIT_REFSET/OBJECT_LEDGER.json \
  --calibration-ladders constraints/rna/RNA01-02/refset/RNA0102_calibration_ladders.json \
  --l16k-fasta work/RNA01-02/forward_input/L16K_UNIQUE_SEQUENCES.fasta \
  --calibration-fasta work/RNA01-02/forward_input/CALIBRATION_UNIQUE_SEQUENCES.fasta \
  --manifest work/RNA01-02/forward_input/MANIFEST.json
```

The outputs contain 49,149 and 384 hash-addressed unique sequences, respectively. Use the L16K FASTA once and the calibration FASTA for two independent Forward passes.

## 2. Forward output contract

Submit three NPZ files: the L16K vectors and two independent calibration-repeat files. Each must contain exactly:

| key | shape | meaning |
|---|---:|---|
| `sequence_sha256` | `(n,)` | unique string IDs aligned to vector rows |
| `vectors` | `(n, d)` | finite real numeric whole-sequence vectors; `d > 0` |

The L16K file must cover exactly the 49,149 unique hashes in `OBJECT_LEDGER.json`. Each calibration file must cover exactly the 384 endpoint hashes referenced by the 128 calibration ladders. Row order is arbitrary when IDs and rows remain aligned. All three files must use the same representation dimension and extraction semantics. Finite zero vectors are accepted under the frozen safe-distance rule.

## 3. Run the official scorer

Prepare the standalone environment in [`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md), then run:

```bash
python constraints/rna/RNA01-02/scoring/score_rna0102_global_ci.py \
  --vectors work/RNA01-02/my_model/L16K_VECTORS.npz \
  --calibration-repeat-1 work/RNA01-02/my_model/CALIBRATION_REPEAT_1.npz \
  --calibration-repeat-2 work/RNA01-02/my_model/CALIBRATION_REPEAT_2.npz \
  --model-id my_model \
  --output work/RNA01-02/my_model/score
```

The scorer locates the RefSet automatically. It validates sequence content and SHA256 identities, edit coordinates, exact Hamming 1/1/2 ladders, population membership, dependency clusters, representation dtype and exact coverage. Population identity is recomputed from content rather than trusted from a declaration. A semantically valid alternative RefSet may be evaluated, but its actual identity is recorded and it is not labeled as the registered L16K population. The scorer then derives `tau_delta` from repeat complete-ladder contrasts using the executed frozen rule, computes every ladder measurement, aggregates outer clusters, and runs the fixed 20,000-draw bootstrap. Users do not select tau, population size, aggregation, bootstrap unit, draw count or thresholds.

## 4. Outputs

| file | contents |
|---|---|
| `INDIVIDUAL_SCORES.tsv` | one row per ladder: IDs, `d01`, `d12`, `d02`, endpoint contrasts and F credits |
| `CLUSTER_SCORE_LEDGER.tsv` | outer-cluster ID, leaf count and `cluster_F` |
| `CALIBRATION.json` | repeat-panel support, calibration policy and derived `tau_delta` |
| `FORMAL_RESULT.json` | `F`, `F_95_CI`, coverage, population/bootstrap identity and input identities |
| `BOOTSTRAP_MEANS.npy` | 20,000 formal cluster-bootstrap F means |
| `OUTPUT_MANIFEST.json` | SHA256 and byte size of each output |

`formal_result=true` means that the released protocol was used. It does not register a third-party result as a paper result. Different pooling, supervised probes or trained readouts may be evaluated, but their interpretation must state that representation design.

## 5. Formal measurement definition

- Public numerical method: [`scoring/RNA01_02_NUMERIC_METHOD_v0.2.6_F_ONLY.json`](scoring/RNA01_02_NUMERIC_METHOD_v0.2.6_F_ONLY.json)
- Exact-Hamming white-box reference: [`comparators/exact-hamming-comparator-current/README.md`](comparators/exact-hamming-comparator-current/README.md)

The public F-only method defines the two-repeat absolute-contrast-difference calibration and scoring procedure.

Paper-model Forward implementations and registered FM outputs belong in `paper_cases/rna/RNA01-02/`; they call this same scorer and are not dependencies of this Constraint.
