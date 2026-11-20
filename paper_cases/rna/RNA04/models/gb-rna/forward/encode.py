from __future__ import annotations
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np,torch
from gb_rna import RNABertModel,RNABertTokenizer,get_vocab_filepath,validate_sequences

def fasta(p):
 out=[]; k=None; s=[]
 for x in Path(p).read_text().splitlines():
  if x.startswith('>'):
   if k is not None: out.append((k,''.join(s)))
   k=x[1:].strip();s=[]
  else:s.append(x.strip().upper())
 if k is not None:out.append((k,''.join(s)))
 return out

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--fasta',required=True);ap.add_argument('--model-dir',required=True);ap.add_argument('--output',required=True);ap.add_argument('--batch-size',type=int,default=2);ap.add_argument('--device',default='cuda');a=ap.parse_args();Path(a.output).parent.mkdir(parents=True,exist_ok=True)
 rows=fasta(a.fasta);tok=RNABertTokenizer(get_vocab_filepath(),version='v2');model=RNABertModel.from_pretrained(a.model_dir).to(a.device).to(torch.bfloat16).eval();valid=torch.tensor([5,6,7,8,9,10],device=a.device);ids=[];out=[];windowed=0;t=time.time()
 def encode(chunks):
  inp=tok(validate_sequences(chunks),max_length=1024,padding=True,truncation=False,return_tensors='pt');inp={k:v.to(a.device) for k,v in inp.items()};h=model(input_ids=inp['input_ids'],attention_mask=inp['attention_mask']).last_hidden_state.float();mask=(inp['input_ids'].unsqueeze(-1)==valid).any(-1)
  counts=mask.sum(1);expected=torch.tensor([len(sequence) for sequence in chunks],device=a.device,dtype=counts.dtype)
  if not torch.equal(counts,expected):raise ValueError('GB-RNA tokenizer did not preserve one nucleotide token per input nucleotide')
  return [(h[i][mask[i]].sum(0),int(mask[i].sum())) for i in range(len(chunks))]
 with torch.inference_mode():
  for i,(ident,seq) in enumerate(rows):
   chunks=[seq[q:q+1022] for q in range(0,len(seq),1022)];windowed+=len(chunks)>1;sums=[]
   for j in range(0,len(chunks),a.batch_size):sums.extend(encode(chunks[j:j+a.batch_size]))
   total=sum((x for x,n in sums),torch.zeros_like(sums[0][0]));n=sum(n for x,n in sums);out.append((total/n).cpu().tolist());ids.append(ident)
   if i%50==0 or i+1==len(rows):print(json.dumps({'done':i+1,'total':len(rows),'windowed_nodes':windowed}),flush=True)
 z=np.stack(out).astype(np.float32)
 if not np.isfinite(z).all():raise ValueError('nonfinite')
 np.savez(a.output,sequence_sha256=np.array(ids),vectors=z)
 Path(str(a.output)+'.metrics.json').write_text(json.dumps({'model_id':'AIDO/GB.RNA-650M','records':len(ids),'embedding_dim':int(z.shape[1]),'pooling':'mean valid nucleotide states across 1022-nt windows; special tokens excluded','windowed_boundary_nodes':windowed,'window_nt':1022,'elapsed_seconds':time.time()-t,'weights_sha256':hashlib.sha256((Path(a.model_dir)/'model.safetensors').read_bytes()).hexdigest()},sort_keys=True)+'\n')
if __name__=='__main__':main()
