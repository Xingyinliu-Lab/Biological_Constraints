import argparse,csv,gzip,hashlib,json,subprocess,sys
from pathlib import Path
FORMAL_FASTA_SHA256='d5d8f39575f5b071278a223875570e32e9b4e3a321b5abba416e28f074445ea4'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def fasta(p):
 d={};k=None;b=[]
 with gzip.open(p,'rt') as f:
  for z in f:
   z=z.strip()
   if z.startswith('>'):
    if k:d[k]=''.join(b)
    k=z[1:];b=[]
   else:b.append(z)
  if k:d[k]=''.join(b)
 return d
def hd(a,b):
 if len(a)!=len(b):raise ValueError('length')
 return sum(x!=y for x,y in zip(a,b))
p=argparse.ArgumentParser();p.add_argument('--model');p.add_argument('--development-leaf',type=Path);p.add_argument('--fasta',type=Path,required=True);p.add_argument('--population',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True);p.add_argument('--oracle',action='store_true');a=p.parse_args();a.outdir.mkdir(parents=True,exist_ok=False)
scorer=Path(__file__).resolve().parents[2]/'scoring'/'score_ps08_formal.py'
if not scorer.is_file():raise FileNotFoundError('PS08 public formal scorer not found: '+str(scorer))
F=['leaf_id','first_edit_id','source_parent_id','source_unit_id','x0_object_id','x1_object_id','x2_object_id','left_near_score','left_far_score','right_near_score','right_far_score']
if a.oracle:
 actual_fasta_sha256=sha(a.fasta)
 seq=fasta(a.fasta);pop=json.load(open(a.population))['canonical_members'];rows=[]
 for r in pop:
  x0,x1,x2=(seq[r[f'x{i}_object_id']] for i in range(3));v=(hd(x0,x1),hd(x0,x2),hd(x2,x1),hd(x2,x0))
  if v!=(1,2,1,2):raise ValueError('oracle invariant')
  rows.append({**{k:r[k] for k in F[:7]},'left_near_score':-v[0],'left_far_score':-v[1],'right_near_score':-v[2],'right_far_score':-v[3]})
 src={'kind':'EXACT_HAMMING_ORACLE','endpoint_fasta_sha256':actual_fasta_sha256,'registered_endpoint_fasta_sha256_match':actual_fasta_sha256==FORMAL_FASTA_SHA256}
else:
 with gzip.open(a.development_leaf,'rt',newline='') as f:q=list(csv.DictReader(f,delimiter='\t'))
 rows=[{**{k:r[k] for k in F[:7]},'left_near_score':r['safe_cosine_01'],'left_far_score':r['safe_cosine_02'],'right_near_score':r['safe_cosine_21'],'right_far_score':r['safe_cosine_20']} for r in q];src={'kind':'EXISTING_FORWARD_REUSE','development_leaf_sha256':sha(a.development_leaf)}
leaf=a.outdir/'FORMAL_LEAF_SCORES.tsv'
with leaf.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=F,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
subprocess.run([sys.executable,str(scorer),'--leaf-ledger',str(leaf),'--population',str(a.population),'--objects-fasta',str(a.fasta),'--outdir',str(a.outdir/'scorer'),'--delta-num','0.0'],check=True)
r=json.load(open(a.outdir/'scorer/SCORE_SUMMARY.json'));r.update({'status':'FORMAL_COMPLETE','method':'EXACT_HAMMING_ORACLE' if a.oracle else a.model,'formal_leaf_sha256':sha(leaf),**src});(a.outdir/'FORMAL_RESULT.json').write_text(json.dumps(r,indent=2)+'\n')
