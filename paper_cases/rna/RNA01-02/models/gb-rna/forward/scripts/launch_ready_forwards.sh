#!/usr/bin/env bash
set -euo pipefail
J=/home/mik/fd/m2_rna/jobs/rna0102_four_model_20261002
F=$J/RNA0102_L16K_unique_endpoints.fasta
if [ ! -e $J/forwards/rnafm_l16k.npz ]; then
 nohup env CUDA_VISIBLE_DEVICES=0 PYTHONPATH=$J/rnafm_source /home/mik/fd/m2_rna/envs/rnafm/bin/python $J/encode_rnafm_global.py --fasta $F --weights /home/mik/fd/m2_rna/models/RNA-FM/RNA-FM_pretrained.pth --source-root $J/rnafm_source --output $J/forwards/rnafm_l16k.npz --batch-size 4 > $J/forwards/rnafm.log 2>&1 & echo $! > $J/forwards/rnafm.pid
fi
if [ ! -e $J/forwards/gb_rna_l16k.npz ]; then
 nohup env CUDA_VISIBLE_DEVICES=1 /home/mik/fd/m2_rna/envs/gb_rna_650m/bin/python $J/encode_gb_rna_global.py --fasta $F --model-dir /home/mik/fd/m2_rna/models/GB.RNA-650M --output $J/forwards/gb_rna_l16k.npz --batch-size 2 > $J/forwards/gb_rna.log 2>&1 & echo $! > $J/forwards/gb_rna.pid
fi
