from __future__ import annotations
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np,torch
from gb_rna import RNABertModel,RNABertTokenizer,get_vocab_filepath,validate_sequences

def fasta(p):
 r=[];i=None;s=[]
 for x in Path(p).read_text().splitlines():
  if x.startswith('>'):
   if i is not None:r.append((i,''.join(s)))
   i=x[1:].strip();s=[]
  else:s.append(x.strip().upper())
 if i is not None:r.append((i,''.join(s)))
 return r

def main():
 a=argparse.ArgumentParser();a.add_argument('--fasta',required=True);a.add_argument('--model-dir',required=True);a.add_argument('--output',required=True);a.add_argument('--batch-size',type=int,default=2);a.add_argument('--limit',type=int);a.add_argument('--device',default='cuda');ns=a.parse_args()
 rows=fasta(ns.fasta);rows=rows[:ns.limit] if ns.limit else rows
 tok=RNABertTokenizer(get_vocab_filepath(),version='v2');model=RNABertModel.from_pretrained(ns.model_dir).to(ns.device).to(torch.bfloat16).eval();vals=[];ids=[];t=time.time(); valid_ids=torch.tensor([5,6,7,8,9,10],device=ns.device)
 with torch.inference_mode():
  for i in range(0,len(rows),ns.batch_size):
   b=rows[i:i+ns.batch_size];seqs=validate_sequences([x[1] for x in b]);inp=tok(seqs,max_length=1024,padding=True,truncation=False,return_tensors='pt');inp={k:v.to(ns.device) for k,v in inp.items()}; h=model(input_ids=inp['input_ids'],attention_mask=inp['attention_mask'],output_hidden_states=False).last_hidden_state.float();mask=(inp['input_ids'].unsqueeze(-1)==valid_ids).any(-1); denom=mask.sum(1,keepdim=True); 
   if (denom==0).any():raise ValueError('no valid nucleotide token')
   z=(h*mask.unsqueeze(-1)).sum(1)/denom
   vals.extend(z.cpu().tolist());ids.extend(x[0] for x in b)
   if i%(ns.batch_size*50)==0 or i+len(b)==len(rows):print(json.dumps({'done':i+len(b),'total':len(rows)}),flush=True)
 z=np.stack(vals).astype(np.float32)
 if not np.isfinite(z).all():raise ValueError('nonfinite')
 Path(ns.output).parent.mkdir(parents=True,exist_ok=True)
 np.savez(ns.output,sequence_sha256=np.array(ids),vectors=z)
 Path(str(ns.output)+'.metrics.json').write_text(json.dumps({'model_id':'AIDO/GB.RNA-650M','records':len(ids),'embedding_dim':int(z.shape[1]),'pooling':'mean valid A/G/C/T/U/N token states; special tokens excluded','elapsed_seconds':time.time()-t,'weights_sha256':hashlib.sha256((Path(ns.model_dir)/'model.safetensors').read_bytes()).hexdigest()},sort_keys=True)+'\n')
if __name__=='__main__':main()

