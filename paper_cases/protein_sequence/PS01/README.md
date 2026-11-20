# PS01 paper implementation cases

Five registered protein-model Forward implementations use the same 2,316-object contract and the unique public PS01 scorer. Large representation files are not duplicated in Git; comparator execution and results remain only in the Constraint release.

| Model | F | 95% CI | Public parity |
|---|---:|---|---|
| carp-640m | 0.41595903488472596 | [0.3188714553921924, 0.510513907175773] | PASS |
| esm2-650m | 0.34465591927753514 | [0.2392500075849009, 0.44596360747174413] | PASS |
| esmc-600m | 0.25669077799334133 | [0.15690225473131578, 0.35634534116676975] | PASS |
| prott5 | 0.6316691748791403 | [0.5426318920774566, 0.7140438806265421] | PASS |
| esm3-seq | 0.6257036401844871 | [0.5433230168289212, 0.701826496462873] | PASS |

Each case retains the Forward code, exact public scoring command, and current ID-level plus aggregate public F outputs. Superseded results are not retained.
