set -euo pipefail
J=/home/mik/fd/m2_rna/jobs/orthrus_4track_rna_20261004
M=/home/mik/fd/m2_rna/models/orthrus-4-track
C=/home/mik/miniconda3/bin/conda
run(){ LD_LIBRARY_PATH="$M/cuda11/nvidia/cuda_runtime/lib:/home/mik/miniconda3/envs/orthrus4track/lib/python3.10/site-packages/torch/lib:${LD_LIBRARY_PATH:-}" CUDA_VISIBLE_DEVICES=1 HF_HUB_OFFLINE=1 HF_HOME="$M/hf_cache" "$C" run -n orthrus4track python "$J/encode_orthrus_4track.py" "$@" --model-id "$M/repository" --device cuda; }
run --fasta /home/mik/fd/m2_rna/jobs/rna04_four_model_20261002/RNA04_unique_endpoints.fasta --output "$J/forwards/rna04_global_gpu1.npz" --mode global --batch-size 32
run --fasta /home/mik/fd/m2_rna/jobs/rna07_four_model_20261003/RNA07_V020_UNIQUE_INPUTS.fasta --output "$J/forwards/rna07_global_gpu1.npz" --mode global --batch-size 8
run --fasta /home/mik/fd/m2_rna/jobs/rna09_four_model_v021_20261004/RNA09_PANEL_ENDPOINTS.fasta --output "$J/forwards/rna09_global_gpu1.npz" --mode global --batch-size 8
