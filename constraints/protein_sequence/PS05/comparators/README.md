# PS05 comparator

PS05 publishes NW+BLOSUM62 as its single formal comparator. It starts from the released PS05 RefSet, produces the standard leaf-measurement ledger, and delegates all relation credits, hierarchy, bootstrap and F/F 95% CI calculation to the single public PS05 scorer.

| Comparator | Method | Canonical F | Reproduction guide |
|---|---|---:|---|
| NW+BLOSUM62 | global Needleman-Wunsch alignment | 0.288054 | [`nw-blosum62/README.md`](nw-blosum62/README.md) |

The comparator directory contains the public implementation, its formal result,
and output identities. Full measurement ledgers are regenerated locally from
the released RefSet when needed.
