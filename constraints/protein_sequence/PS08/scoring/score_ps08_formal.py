#!/usr/bin/env python3
import argparse,csv,gzip,json,math
from collections import defaultdict
from pathlib import Path
import numpy as np
from ps08_contract import ID_FIELDS, load_population, sha256 as sha
SEED=20261002;DRAWS=20000;DELTA=0.0
SCORE=['left_near_score','left_far_score','right_near_score','right_far_score']
def rows(p):
 op=gzip.open if p.suffix=='.gz' else open
 with op(p,'rt',newline='') as f:
  out=list(csv.DictReader(f,delimiter='\t'))
 if not out or any(set(ID_FIELDS+tuple(SCORE))-set(r) for r in out):raise ValueError('formal leaf schema failure')
 if len({r['leaf_id'] for r in out})!=len(out):raise ValueError('duplicate leaf_id')
 for r in out:
  v=[]
  for k in SCORE:
   x=float(r[k]);
   if not math.isfinite(x):raise ValueError('nonfinite score')
   v.append(x)
  dl=v[0]-v[1];dr=v[2]-v[3];r['delta_left']=dl;r['delta_right']=dr;r['F_credit_left']=1. if dl>DELTA else -1. if dl<-DELTA else 0.;r['F_credit_right']=1. if dr>DELTA else -1. if dr<-DELTA else 0.;r['F_value']=(r['F_credit_left']+r['F_credit_right'])/2
 return out
def validate(p,fasta_path,r):
 expect,provenance=load_population(p,fasta_path);e={z['leaf_id']:z for z in expect}
 if len(r)!=len(expect) or set(e)!={z['leaf_id'] for z in r}:raise ValueError('population mismatch')
 for z in r:
  q=e[z['leaf_id']]
  for k in ID_FIELDS:
   if z[k]!=q[k]:raise ValueError('identity drift '+k)
 return provenance
def aggregate(r):
 first=defaultdict(list);meta={}
 for z in r:first[z['first_edit_id']].append(z['F_value']);meta[z['first_edit_id']]=(z['source_parent_id'],z['source_unit_id'])
 first_score={k:math.fsum(v)/len(v) for k,v in first.items()}
 parent=defaultdict(list);punit={}
 for k,v in first_score.items():p,u=meta[k];parent[p].append(v);punit[p]=u
 parent_score={k:math.fsum(v)/len(v) for k,v in parent.items()}
 unit=defaultdict(list)
 for p,v in parent_score.items():unit[punit[p]].append(v)
 unit_score={k:math.fsum(v)/len(v) for k,v in unit.items()}
 F=math.fsum(unit_score.values())/len(unit_score)
 for z in r:
  f=z['first_edit_id'];p=z['source_parent_id'];u=z['source_unit_id']
  z['first_edit_F']=first_score[f];z['source_parent_F']=parent_score[p];z['source_unit_F']=unit_score[u]
  z['overall_contribution_F']=z['F_value']/len(first[f])/len(parent[p])/len(unit[u])/len(unit_score)
 return unit_score,F,{'leaf':len(r),'first_edit':len(first),'source_parent':len(parent),'source_unit':len(unit),'overall':1}
def write_individual(path,r):
 fields=[*ID_FIELDS,*SCORE,'delta_left','delta_right','F_credit_left','F_credit_right','leaf_F','first_edit_F','source_parent_F','source_unit_F','overall_contribution_F']
 with path.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
  for z in r:w.writerow({**{k:z[k] for k in ID_FIELDS+tuple(SCORE)},'delta_left':repr(z['delta_left']),'delta_right':repr(z['delta_right']),'F_credit_left':repr(z['F_credit_left']),'F_credit_right':repr(z['F_credit_right']),'leaf_F':repr(z['F_value']),'first_edit_F':repr(z['first_edit_F']),'source_parent_F':repr(z['source_parent_F']),'source_unit_F':repr(z['source_unit_F']),'overall_contribution_F':repr(z['overall_contribution_F'])})
def main():
 a=argparse.ArgumentParser();a.add_argument('--leaf-ledger',type=Path,required=True);a.add_argument('--population',type=Path,required=True);a.add_argument('--objects-fasta',type=Path,required=True);a.add_argument('--outdir',type=Path,required=True);a.add_argument('--delta-num',type=float,default=0.0,help='Frozen formal value; nonzero values are rejected.');x=a.parse_args()
 if x.delta_num!=0.:raise ValueError('PS08 frozen delta_num is 0')
 r=rows(x.leaf_ledger);population=validate(x.population,x.objects_fasta,r);u,F,counts=aggregate(r);vals=np.array([u[k] for k in sorted(u)],dtype=np.float64);rng=np.random.Generator(np.random.PCG64(SEED));draws=np.empty(DRAWS,dtype=np.float64)
 for start in range(0,DRAWS,500):
  n=min(500,DRAWS-start);idx=rng.integers(0,len(vals),size=(n,len(vals)));draws[start:start+n]=vals[idx].mean(axis=1)
 ci=tuple(float(q) for q in np.quantile(draws,[.025,.975],method='linear'))
 x.outdir.mkdir(parents=True,exist_ok=True);bootstrap_path=x.outdir/'BOOTSTRAP_F.float64le.bin';individual_path=x.outdir/'INDIVIDUAL_SCORES.tsv';summary_path=x.outdir/'SCORE_SUMMARY.json';bootstrap_path.write_bytes(draws.tobytes());write_individual(individual_path,r);semantic=population['semantic_population_sha256'];res={'kind':'PS08_FORMAL_SCORER_RESULT','formal_result':True,'formal_interface':'F_ONLY','score_transform_version':'PS08_v0.1.2_F_ONLY','F':F,'F_95_CI':list(ci),'counts':counts,'seed':SEED,'draws':DRAWS,'delta_num':0.,'bootstrap':'source_unit with replacement','leaf_ledger_sha256':sha(x.leaf_ledger),'population_identity_sha256':semantic,'population_validation':population,'individual_scores':{'path':str(individual_path),'sha256':sha(individual_path),'rows':len(r)},'bootstrap_F_float64le_sha256':sha(bootstrap_path)};summary_path.write_text(json.dumps(res,indent=2)+'\n',encoding='utf-8');manifest={p.name:{'bytes':p.stat().st_size,'sha256':sha(p)} for p in (individual_path,bootstrap_path,summary_path)};(x.outdir/'OUTPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
