# RNA06 ViennaRNA partner-probability comparator

This comparator derives base-pair partner probabilities from the released RNA06 sequences using ViennaRNA and the frozen Turner-2004 parameter file `RNA_Turner2004_ViennaRNA_2.4.17.par`. It is evaluated on the explicitly compatible subset: 1,274 of 1,588 formal reciprocal edges; 314 edges are reported as out of scope.

`vienna_partition_partner_probability.py` produces the raw comparator measurement. `RUN_COMPARATOR.py` is the portable released entry point that applies the frozen compatibility rules and creates the 12,182 candidate-level measurements. It then calls the same public RNA06 aggregation implementation used for model representations; no comparator-specific scorer is maintained. The registered result is F `0.9473389148952038`, with its F 95% CI in `FORMAL_RESULT.json`.

The comparator-specific compatible ledgers are under `refset/`; they do not replace the full formal model population in the parent Constraint RefSet.

## Run

Install the versions listed in `ENVIRONMENT.md`, then run from the repository root:

```bash
python constraints/rna/RNA06/comparators/viennarna-partner-probability/RUN_COMPARATOR.py \
  --output work/RNA06/viennarna-partner-probability
```

The command defaults to the released RNA06 RefSet, including its frozen pair-state map. It writes the compatibility audit, candidate-, endpoint-, edge-, parent- and cluster-level F ledgers, regenerated Turner parameter artifact, formal aggregate result and output manifest.
