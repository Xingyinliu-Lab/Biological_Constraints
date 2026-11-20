set -euo pipefail
R04=/home/mik/fd/m2_rna/jobs/rna04_four_model_20261002/forwards
J=/home/mik/fd/m2_rna/jobs/rna0102_four_model_20261002
W=$J/rna0102_overlap_remediation
for x in rnafm gb_rna rinalmo ernie_rna; do while [ ! -s "$R04/$x.npz" ]; do sleep 60; done; done
python3 "$W/prepare_rna0102_boundary_overlap.py" --fasta "$J/RNA0102_L16K_unique_endpoints.fasta" --out "$W/boundary_21.fasta" --summary "$W/boundary_summary.json"
nohup env CUDA_VISIBLE_DEVICES=0 PYTHONPATH="$J/rnafm_source" /home/mik/fd/m2_rna/envs/rnafm/bin/python "$W/encode_rnafm_boundary_overlap.py" --fasta "$W/boundary_21.fasta" --weights /home/mik/fd/m2_rna/models/RNA-FM/RNA-FM_pretrained.pth --source-root "$J/rnafm_source" --output "$W/rnafm_boundary_overlap.npz" > "$W/rnafm_boundary.log" 2>&1 & p0=$!
nohup env CUDA_VISIBLE_DEVICES=1 /home/mik/fd/m2_rna/envs/ernie_rna/bin/python "$W/encode_ernie_boundary_overlap.py" --fasta "$W/boundary_21.fasta" --model-dir /home/mik/fd/m2_rna/models/ERNIE-RNA-HF --output "$W/ernie_boundary_overlap.npz" > "$W/ernie_boundary.log" 2>&1 & p1=$!
wait $p0; wait $p1
/home/mik/fd/m2_rna/envs/rnafm/bin/python "$W/splice_rna0102_boundary_vectors.py" --base "$J/forwards/rnafm_l16k.npz" --patch "$W/rnafm_boundary_overlap.npz" --out "$W/rnafm_l16k_overlap.npz"
/home/mik/fd/m2_rna/envs/rnafm/bin/python "$W/splice_rna0102_boundary_vectors.py" --base "$J/forwards/ernie_l16k.npz" --patch "$W/ernie_boundary_overlap.npz" --out "$W/ernie_l16k_overlap.npz"
touch "$W/COMPLETE"