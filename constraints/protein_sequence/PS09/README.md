# PS09: disorder-to-order transition readability

PS09 asks whether model representations make experimentally disordered residues distinguishable from coordinate-supported ordered residues in matched protein blocks. The released method uses a frozen five-fold, block-state-weighted logistic readout and evaluates only held-out endpoints. It measures information accessible through this prespecified low-complexity readout; it is not a probe-free cosine assay.

The RefSet is `PS09_ASSAY_REFSET_RELEASE_v1.0.0_20261008`; the formal input binding is [`scoring/PS09_NUMERIC_INPUT_BINDING.json`](scoring/PS09_NUMERIC_INPUT_BINDING.json).

## 1. RefSet to model Forward input

Use these released assets:

- [`refset/SELECTED_PROTEINS.faa`](refset/SELECTED_PROTEINS.faa): 81 protein sequences to encode.
- [`refset/FORMAL_BLOCK_MEMBERSHIP.tsv.gz`](refset/FORMAL_BLOCK_MEMBERSHIP.tsv.gz): 15,708 formal residue endpoints with block, state, protein accession, one-based residue position and homology-cluster identity.
- [`refset/PS09_FROZEN_CLUSTER_FOLD_ASSIGNMENT.tsv`](refset/PS09_FROZEN_CLUSTER_FOLD_ASSIGNMENT.tsv): frozen five-fold assignment. The scorer consumes this file; the user must not refit or alter folds.

Encode every sequence in the FASTA with the model being evaluated. For each membership row, extract the representation corresponding to `uniprot_position_1based`. Token-to-residue mapping, window stitching and model-layer selection belong to the model Forward and must be applied consistently. Do not change the released sequences, endpoint population, state labels, block membership or fold assignment.

For every reported evaluation, retain Forward mapping evidence: model and tokenizer identity, token-to-residue/span rule, windowing procedure, treatment of tokens crossing the target residue, and final residue-vector readout/pooling rule. The public scorer validates released endpoint IDs and numerical inputs; it does not load or independently validate a third-party tokenizer.

The scorer endpoint ID for a membership row is:

```text
<block_id>|<state>|<uniprot_position_1based>
```

If the same biological residue participates in more than one formal block, export the same Forward-derived residue vector under each required formal endpoint ID.

## 2. Model Forward output contract

The public scorer accepts either NPZ or Safetensors.

NPZ must contain exactly:

```text
object_ids  string array, shape (15708,)
vectors     real numeric array, shape (15708, d), d > 0
```

Safetensors must contain exactly the 15,708 formal endpoint IDs as keys, each mapped to one numeric vector of the same positive dimension.

For both formats:

- IDs must be unique and exactly cover the released membership; row/key order is irrelevant.
- All values must be finite.
- Complex, Boolean, object and non-numeric arrays are rejected; the scorer never discards an imaginary component by conversion.
- No extra or missing endpoint is allowed.
- The user supplies representations only. Scaling, frozen cross-fit, class readout, pair F credits, cluster aggregation and bootstrap are performed by PS09.

## 3. Run the official scorer

Install the small scoring environment described in [`scoring/ENVIRONMENT.md`](scoring/ENVIRONMENT.md), then run one command from the repository root:

```bash
python constraints/protein_sequence/PS09/scoring/score_ps09.py \
  --representations work/PS09/my_model/representations.npz \
  --model-id my_model \
  --output-dir work/PS09/my_model/score
```

Use the same command for a Safetensors input by changing the `--representations` path. The scorer automatically loads the released RefSet and all frozen scientific parameters. Before fitting, it verifies FASTA accession and sequence closure, one-based residue range and amino-acid identity, matched-block endpoint counts, MMseqs2 cluster membership, and fold/cluster/block-count consistency. `--refset` exists only for an explicitly supplied compatible release; the result records the assets actually used.

## 4. Official outputs

The output directory contains:

- `HELDOUT_ENDPOINT_SCORES.tsv`: held-out disorder probability for every formal endpoint and its frozen fold.
- `INDIVIDUAL_SCORES.tsv`: endpoint ID, block, state, cluster, fold and held-out raw score.
- `PAIR_CREDIT_LEDGER.tsv.gz`: every disorder-versus-order endpoint comparison, raw scores, margin, tie and F_credit.
- `BLOCK_SCORES.tsv`: block-level F computed from disorder versus ordered endpoints.
- `CLUSTER_SCORES.tsv`: block-mean F for each homology cluster; these are the formal bootstrap units.
- `SCORE_SUMMARY.json`: F, F 95% CI, coverage, bootstrap settings, model ID and input identities.
- `SCORE_SUMMARY.json.bootstrap.float64le.bin`: the 20,000 frozen bootstrap draws.
- `OUTPUT_MANIFEST.json`: file sizes and SHA256 identities for the run.

The fixed formal settings are five frozen folds, `LogisticRegression(C=1.0, solver="lbfgs", max_iter=1000)`, training-only block-state-balanced weights and scaling, homology-cluster aggregation, PCG64 seed `20261008`, 20,000 draws and linear 2.5/97.5 percentiles. Users do not choose these scoring decisions.

`formal_result=true` means the released PS09 protocol completed successfully. It does not register a third-party run as a paper result.

## 5. Interpretation

Within each block, pair outcomes are reported as `F_credit` (win = 1, tie = 0, loss = -1), equivalent to `F = 2AUC - 1`. F is the formal mean cluster-level readability score; chance-level ranking corresponds to `F=0`. The confidence interval resamples homology clusters, not individual residues. Coverage must show all 15,708 endpoints, 1,030 blocks and 78 clusters before a result is treated as complete.

The paper model implementations and registered FM outputs are maintained outside this independent Constraint package. Every paper model and every third-party model enters this same public scorer through the Forward output contract above.

## 6. Built-in comparator

The IUPred2A comparator, its execution method and canonical ID-level and aggregate outputs are in [`comparators/iupred2a/`](comparators/iupred2a/). It uses the same PS09 aggregation scorer after producing its endpoint scores.
