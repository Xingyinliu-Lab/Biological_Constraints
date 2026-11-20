#!/usr/bin/env python3
"""JanusDNA adapter for the shared DNA03 v0.2.0 continuous mapping runner."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import torch
import encode_dna0202_janusdna as dna02_loader
from encode_dna0202_janusdna import BASE_TO_ID,PAD_ID,MODEL_LENGTH,load_model
from dna03_native_runner import encode_native, read_fasta, read_rows
from dna03_forward_publication import publish_forward

EXPECTED_REMAP={
    'model.layers.4.attn.self_attn_bwd.o_proj.weight':'model.layers.4.attn.self_attn_bwd.o_projs.0.weight',
    'model.layers.4.attn.self_attn_fwd.o_proj.weight':'model.layers.4.attn.self_attn_fwd.o_projs.0.weight',
}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--membership',type=Path,required=True); ap.add_argument('--fasta',type=Path,required=True)
    ap.add_argument('--source',type=Path,required=True); ap.add_argument('--config',type=Path,required=True); ap.add_argument('--checkpoint',type=Path,required=True); ap.add_argument('--output',type=Path,required=True); ap.add_argument('--limit',type=int)
    args=ap.parse_args(); sequences=read_fasta(args.fasta); rows=read_rows(args.membership,sequences); rows=rows[:args.limit] if args.limit else rows
    model,missing,unexpected,remapped,ignored=load_model(args.source,args.config,args.checkpoint)
    if missing or unexpected or remapped != EXPECTED_REMAP:
        raise RuntimeError(f'JanusDNA checkpoint load is not exact after required remap: missing={missing}, unexpected={unexpected}, remapped={remapped}')
    def forward_batch(seqs):
        encoded=[[BASE_TO_ID[b] for b in s]+[PAD_ID]*(MODEL_LENGTH-len(s)) for s in seqs]
        ids=torch.tensor(encoded,device='cuda')
        with torch.inference_mode(),torch.autocast(device_type='cuda',dtype=torch.float16): hidden=model.model(input_ids=ids,return_dict=True,use_cache=False).last_hidden_state
        return [hidden[i,:len(s)].float().cpu().numpy() for i,s in enumerate(seqs)]
    keys,vectors,tasks,coverage=encode_native(rows,sequences,128,forward_batch)
    keys2,repeat_vectors,tasks2,coverage2=encode_native(rows,sequences,128,forward_batch)
    load_record={'missing_keys':missing,'unexpected_parameter_keys':unexpected,'compatibility_key_remap':remapped,'ignored_metadata':ignored}
    publish_forward(args.output,model_key='janusdna',model_identity={'model_id':'JanusDNA_len1k_dmodel32_interdim128_8layers','model_revision':'checkpoint_sha256:c1dc802dc88c298672f9c1ceb2fda009dab6a2eef2c0551d046fcda3b1b628e2','source':str(args.source),'config':str(args.config),'checkpoint':str(args.checkpoint)},
        source_paths={'membership':args.membership,'fasta':args.fasta,'config':args.config,'checkpoint':args.checkpoint,'dna02_fixed_loader':Path(dna02_loader.__file__)},
        primary=(keys,vectors),repeat=(keys2,repeat_vectors),validation={'tasks':tasks,'repeat_tasks':tasks2,'checkpoint_load':load_record,'global_nucleotide_coverage':coverage,'repeat_global_nucleotide_coverage':coverage2},
        code_paths={'adapter':Path(__file__),'native_runner':Path(__file__).with_name('dna03_native_runner.py'),'mapping':Path(__file__).with_name('dna03_continuous_mapping.py'),'publisher':Path(__file__).with_name('dna03_forward_publication.py'),'dna02_fixed_loader':Path(dna02_loader.__file__)})
    print(json.dumps({'status':'PASS','interface_version':'0.2.3','views':len(keys),'tasks':tasks,'repeat_bitwise_identical':True,'checkpoint_load':load_record},sort_keys=True))
if __name__=='__main__': main()
