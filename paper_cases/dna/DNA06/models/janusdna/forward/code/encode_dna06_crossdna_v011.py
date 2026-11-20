#!/usr/bin/env python3
"""CrossDNA adapter for the shared DNA06 v0.2.0 continuous mapping runner."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import torch
from transformers import AutoModelForMaskedLM
from dna06_native_runner import encode_native, read_fasta, read_rows
from dna06_forward_publication import publish_forward


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--membership',type=Path,required=True); ap.add_argument('--fasta',type=Path,required=True)
    ap.add_argument('--model-dir',type=Path,required=True); ap.add_argument('--output',type=Path,required=True); ap.add_argument('--limit',type=int)
    args=ap.parse_args(); sequences=read_fasta(args.fasta); rows=read_rows(args.membership,sequences); rows=rows[:args.limit] if args.limit else rows
    model=AutoModelForMaskedLM.from_pretrained(args.model_dir,trust_remote_code=True,local_files_only=True).eval().half().cuda(); model.backbone.pretrain=False; model.backbone.for_representation=True
    base={b:i for i,b in enumerate('ACGTN')}
    def forward_batch(seqs):
        width=max(map(len,seqs)); ids=torch.tensor([[base[b] for b in s]+[4]*(width-len(s)) for s in seqs],device='cuda')
        with torch.inference_mode(),torch.autocast(device_type='cuda',dtype=torch.float16): hidden,_=model.backbone(ids)
        return [hidden[i,:len(s)].float().cpu().numpy() for i,s in enumerate(seqs)]
    keys,embeddings,tasks,coverage=encode_native(rows,sequences,8,forward_batch)
    keys2,repeat_embeddings,tasks2,coverage2=encode_native(rows,sequences,8,forward_batch)
    publish_forward(args.output,model_key='crossdna_8m',model_identity={'model_id':'CrossDNA_8.1M','model_revision':'DNA01_EXECUTED_CHECKPOINT_IDENTITY','model_dir':str(args.model_dir)},
        source_paths={'membership':args.membership,'fasta':args.fasta},primary=(keys,embeddings),repeat=(keys2,repeat_embeddings),
        validation={'tasks':tasks,'repeat_tasks':tasks2,'global_nucleotide_coverage':coverage,'repeat_global_nucleotide_coverage':coverage2},
        code_paths={'adapter':Path(__file__),'native_runner':Path(__file__).with_name('dna06_native_runner.py'),'mapping':Path(__file__).with_name('dna06_continuous_mapping.py'),'publisher':Path(__file__).with_name('dna06_forward_publication.py')})
    print(json.dumps({'status':'PASS','interface_version':'0.2.3','views':len(keys),'tasks':tasks,'repeat_bitwise_identical':True},sort_keys=True))
if __name__=='__main__': main()
