# PS09 public scoring

[`score_ps09.py`](score_ps09.py) is the only user-facing PS09 scorer. It accepts model representations, validates the complete biological RefSet closure, performs the frozen five-fold readout and writes endpoint-, block- and cluster-level F plus F 95% CI.

```bash
python constraints/protein_sequence/PS09/scoring/score_ps09.py \
  --representations work/PS09/my_model/representations.npz \
  --model-id my_model \
  --output-dir work/PS09/my_model/score
```

See [`../README.md`](../README.md) for the complete Forward input/output contract and output schemas, and [`ENVIRONMENT.md`](ENVIRONMENT.md) for dependencies.

## Internal modules

- `representations_to_leaf.py` implements the frozen representation validation and five-fold held-out readout.
- `ps09_refset.py` validates FASTA residues and coordinates plus membership, block, cluster and fold closure.
- `score_ps09_formal.py` implements the frozen endpoint-to-block-to-cluster aggregation and bootstrap.

They are consecutive modules of the one public PS09 pipeline, not alternative scorers. Paper cases and third-party models use the same `score_ps09.py` entry. Paper-specific checkpoint and Forward assertions remain outside this package and do not alter scoring.

The formal parameters are not exposed as user-selectable CLI options. Supplying a compatible alternate RefSet is allowed through `--refset`; the output records the actual asset identities and is comparable to the released paper results only when the same formal population and method are used.

For each fold, block-state weights must sum to the number of complete training
blocks derived from that RefSet. The official population still yields the
frozen historical value 824 in every fold; compatible populations are checked
against their own complete training-block population rather than that number.
