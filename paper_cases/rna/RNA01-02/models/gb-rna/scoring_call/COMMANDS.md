# RNA01-02 / gb-rna public scoring call

Run from the repository root after Forward produces the three contract files:

```bash
python constraints/rna/RNA01-02/scoring/score_rna0102_global_ci.py \
  --vectors work/RNA01-02/gb-rna/L16K_VECTORS.npz \
  --calibration-repeat-1 work/RNA01-02/gb-rna/CALIBRATION_REPEAT_1.npz \
  --calibration-repeat-2 work/RNA01-02/gb-rna/CALIBRATION_REPEAT_2.npz \
  --model-id AIDO-GB.RNA-650M \
  --output work/RNA01-02/gb-rna/score
```

Canonical F-only outputs are stored under `output/RNA01-02_v0.2.6_F_ONLY/`.
