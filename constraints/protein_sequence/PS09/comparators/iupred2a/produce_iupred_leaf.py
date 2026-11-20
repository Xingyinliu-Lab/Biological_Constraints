import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
from pathlib import Path

ap=argparse.ArgumentParser();ap.add_argument('--refset-dir',type=Path,required=True);ap.add_argument('--iupred-dir',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args()
rel=a.refset_dir;run=a.output_dir;app=a.iupred_dir;run.mkdir(parents=True,exist_ok=True);(run/'tmp').mkdir(parents=True,exist_ok=True)
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for x in iter(lambda:f.read(1<<20),b''):
            h.update(x)
    return h.hexdigest()
seq={}; name=None; buf=[]
for line in (rel/'SELECTED_PROTEINS.faa').read_text().splitlines():
    if line.startswith('>'):
        if name: seq[name]=''.join(buf)
        name=line[1:].split()[0]; buf=[]
    else:
        buf.append(line.strip())
if name: seq[name]=''.join(buf)
profiles={}; rows=[]
for i,(acc,s) in enumerate(sorted(seq.items()),1):
    f=run/'tmp'/f'{acc}.faa'
    f.write_text(f'>{acc}\n{s}\n')
    out=subprocess.check_output(['python3',str(app/'iupred2a.py'),'-d',str(app),str(f),'long'],text=True)
    q=[]
    for ln in out.splitlines():
        m=re.match(r'^(\d+)\s+([A-Z])\s+([0-9.eE+-]+)$',ln)
        if m: q.append((int(m.group(1)),m.group(2),float(m.group(3))))
    if len(q)!=len(s) or any(pos!=j+1 or aa!=s[j] or not 0<=v<=1 for j,(pos,aa,v) in enumerate(q)):
        raise RuntimeError(f'IUPred profile invariant failed: {acc}')
    profiles[acc]={pos:v for pos,aa,v in q}
    rows.extend((acc,pos,aa,v) for pos,aa,v in q)
with gzip.open(run/'IUPRED2A_LONG_PROFILE.tsv.gz','wt',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['uniprot_accession','position_1based','amino_acid','iupred_long_score']); w.writerows(rows)
leaf=[]
with gzip.open(rel/'FORMAL_BLOCK_MEMBERSHIP.tsv.gz','rt') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        acc=r['uniprot_accession']; pos=int(r['uniprot_position_1based'])
        leaf.append({'endpoint_id':f"{r['block_id']}|{r['state']}|{pos}",'score':profiles[acc][pos]})
with open(run/'FORMAL_LEAF_SCORES.tsv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['endpoint_id','score'],delimiter='\t'); w.writeheader(); w.writerows(leaf)
binding={'method':'IUPred2A','mode':'long','source_repo':'https://github.com/doszilab/iupred2a.git','commit':subprocess.check_output(['git','-C',str(app.parent),'rev-parse','HEAD'],text=True).strip(),'implementation_path':str(app/'iupred2a.py'),'implementation_sha256':sha(app/'iupred2a.py'),'input_fasta_sha256':sha(rel/'SELECTED_PROTEINS.faa'),'membership_sha256':sha(rel/'FORMAL_BLOCK_MEMBERSHIP.tsv.gz'),'profile_sha256':sha(run/'IUPRED2A_LONG_PROFILE.tsv.gz'),'leaf_sha256':sha(run/'FORMAL_LEAF_SCORES.tsv'),'profile_count':len(profiles),'profile_rows':len(rows),'leaf_rows':len(leaf),'query_time_database_or_retrieval':False}
(run/'IUPRED_EXECUTION_BINDING.json').write_text(json.dumps(binding,indent=2)+'\n')
print(json.dumps(binding,sort_keys=True))
