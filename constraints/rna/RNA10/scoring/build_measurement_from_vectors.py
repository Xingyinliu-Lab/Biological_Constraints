#!/usr/bin/env python3
"""RNA10 signed-delta LOO ridge measurement builder."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from rna10_measurement_core import build_measurement_records

def file_sha256(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def jsonl(path):
 return [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]
def unique_map(rows,key,value,label):
 out={}
 for row in rows:
  item_id=row[key]
  if item_id in out:raise ValueError(f'duplicate {label}: {item_id}')
  out[item_id]=row[value]
 return out
def main():
 root=Path(__file__).resolve().parent
 a=argparse.ArgumentParser();a.add_argument('--truth-relations',type=Path,required=True);a.add_argument('--objects-wt',type=Path,required=True);a.add_argument('--objects-mutant',type=Path,required=True);a.add_argument('--vectors',type=Path,required=True);a.add_argument('--runplan',type=Path,default=root/'RNA10_FORMAL_RUNPLAN.json');a.add_argument('--binding',type=Path,default=root/'RNA10_FORMAL_SCORING_BINDING.json');a.add_argument('--out',type=Path,required=True);a.add_argument('--model',required=True);x=a.parse_args();x.out.mkdir(parents=True,exist_ok=True)
 rows=jsonl(x.truth_relations);plan=json.loads(x.runplan.read_text());binding=json.loads(x.binding.read_text())
 wt=unique_map(jsonl(x.objects_wt),'wt_object_id','sequence_sha256','WT object ID')
 mut=unique_map(jsonl(x.objects_mutant),'mutation_id','sequence_sha256','mutation ID')
 for r in rows:
  r['wt_sequence_sha256']=wt[r['wt_object_id']];r['mut_sequence_sha256']=mut[r['mutation_id']]
 needed={r['wt_sequence_sha256'] for r in rows}|{r['mut_sequence_sha256'] for r in rows}
 with np.load(x.vectors,allow_pickle=False) as v:
  if set(v.files)!={'sequence_sha256','vectors'}:raise ValueError('vector NPZ must contain exactly sequence_sha256 and vectors')
  ids=[str(k) for k in v['sequence_sha256'].tolist()];raw_vectors=v['vectors']
  if np.issubdtype(raw_vectors.dtype,np.bool_) or not np.issubdtype(raw_vectors.dtype,np.number) or np.issubdtype(raw_vectors.dtype,np.complexfloating):raise ValueError('vectors must use a real numeric dtype')
  vectors=np.asarray(raw_vectors,dtype=np.float64)
 if vectors.ndim!=2 or vectors.shape[1]<=0 or vectors.shape[0]!=len(ids):raise ValueError('invalid vector matrix shape')
 if len(ids)!=len(set(ids)):raise ValueError('duplicate sequence_sha256 in vector NPZ')
 if set(ids)!=needed:raise ValueError('vector object universe mismatch')
 if not np.isfinite(vectors).all():raise ValueError('nonfinite representation')
 vm={k:z for k,z in zip(ids,vectors)}
 X=np.stack([vm[r['mut_sequence_sha256']]-vm[r['wt_sequence_sha256']] for r in rows]).astype(float)
 if not np.isfinite(X).all():raise ValueError('nonfinite representation delta')
 records,outer=build_measurement_records(rows,X,plan)
 runplan_sha=hashlib.sha256(json.dumps(plan,sort_keys=True).encode()).hexdigest()
 source_inputs={'vectors_sha256':file_sha256(x.vectors),'truth_relations_sha256':file_sha256(x.truth_relations),'objects_wt_sha256':file_sha256(x.objects_wt),'objects_mutant_sha256':file_sha256(x.objects_mutant),'runplan_file_sha256':file_sha256(x.runplan),'binding_file_sha256':file_sha256(x.binding)}
 out={'artifact_type':'C_RNA10_MEASUREMENT_LEDGER','model_id':x.model,'measurement_artifact_release_id':f'RNA10_{x.model}_REGENERATED_MEASUREMENT','runplan_sha256':runplan_sha,'population_identity_sha256':binding['population_identity_sha256'],'numeric_method_sha256':binding['numeric_method_sha256'],'split_release_sha256':binding['split_release_sha256'],'representation_dimension':int(X.shape[1]),'delta_definition':'z_mut-z_WT','label_encoding':{'PRESERVING':0,'DISRUPTING':1},'source_inputs':source_inputs,'records':records,'outer_LOO':outer}
 (x.out/'MEASUREMENT_LEDGER.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
