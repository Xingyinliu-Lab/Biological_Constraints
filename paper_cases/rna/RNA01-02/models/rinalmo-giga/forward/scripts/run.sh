#!/usr/bin/env bash
set -euo pipefail
exec > "/home/mik/fd/m2_rna/jobs/rna_rinalmo_official_closure_20261005/rna0102/rinalmo_official.log" 2>&1
CUDA_VISIBLE_DEVICES=0 PYTHONPATH="/home/mik/fd/m2_rna/jobs/rna09_closure_20261005/rinalmo_official/repo" "/home/mik/fd/m2_rna/jobs/rna09_closure_20261005/rinalmo_official/venv/bin/python" "/home/mik/fd/m2_rna/jobs/rna_rinalmo_official_closure_20261005/rna0102/encode_official_rinalmo_global.py" --fasta "/home/mik/fd/m2_rna/jobs/rna0102_four_model_20261002/RNA0102_L16K_unique_endpoints.fasta" --repo "/home/mik/fd/m2_rna/jobs/rna09_closure_20261005/rinalmo_official/repo" --out "/home/mik/fd/m2_rna/jobs/rna_rinalmo_official_closure_20261005/rna0102/rinalmo_official.npz" --device cuda --batch-size 1
touch "/home/mik/fd/m2_rna/jobs/rna_rinalmo_official_closure_20261005/rna0102/RUN_COMPLETE"
