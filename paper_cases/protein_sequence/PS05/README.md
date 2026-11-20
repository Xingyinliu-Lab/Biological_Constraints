# PS05 paper implementation cases

Five paper models use the same 321-object Forward contract and the single
public PS05 scorer. Large representation files are not distributed; they are
regenerated from the released FASTA and the recorded model checkpoint.

The four non-ESM3 cases use the parameterized public adapter
[`shared_forward/ps05_forward_public.py`](shared_forward/ps05_forward_public.py).
It implements the model loading, registered representation layer and
residue-mean pooling recorded by the formal PS05 model audit. ESM3-seq retains
its separately executed sequence-track adapter under its case directory. Each
model directory documents the exact Forward command, scorer handoff and current
formal result. Comparator implementations and results belong to the PS05
Constraint release.

Paper checkpoint and Forward-output hashes are provenance records for the
reported executions. The public scorer does not use those hashes as an
admission gate: regenerated or user-supplied representations are accepted
when they satisfy the released object-ID, schema and numerical-validity
requirements.

Models: `esm2-650m`, `prott5`, `esmc-600m`, `carp-640m`, `esm3-seq`.
