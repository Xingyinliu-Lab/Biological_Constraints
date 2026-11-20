#!/usr/bin/env python3
"""Single public PS05 entry from representations to formal outputs."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
REFSET=HERE.parent/'refset'
def sha(path:Path)->str:
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1<<20),b''): h.update(block)
 return h.hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--representations',type=Path,required=True);p.add_argument('--representation-index',type=Path);p.add_argument('--model-id',required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
 raw=a.output_dir/'RAW_MEASUREMENTS.tsv';individual=a.output_dir/'INDIVIDUAL_SCORES.tsv.gz'
 cmd=[sys.executable,str(HERE/'representations_to_leaf.py'),'--relations',str(REFSET/'FORMAL_RELATION_LEDGER.tsv'),'--planned-population',str(REFSET/'PLANNED_SCOREABLE_POPULATION.json'),'--objects-fasta',str(REFSET/'PS05_ENDPOINT_SEQUENCES.faa'),'--representations',str(a.representations),'--model-id',a.model_id,'--output',str(raw)]
 if a.representation_index:cmd.extend(['--representation-index',str(a.representation_index)])
 subprocess.run(cmd,check=True)
 subprocess.run([sys.executable,str(HERE/'score_ps05_formal.py'),'--contract-root',str(REFSET),'--contract-binding',str(HERE/'CONTRACT_BINDING.json'),'--planned-population',str(REFSET/'PLANNED_SCOREABLE_POPULATION.json'),'--leaf-ledger',str(raw),'--delta-num','0','--reuse-ledger',str(REFSET/'ENDPOINT_REUSE_LEDGER.tsv'),'--individual-scores',str(individual),'--output-dir',str(a.output_dir)],check=True)
 summary=a.output_dir/'SCORE_SUMMARY.json';data=json.loads(summary.read_text());data['model_id']=a.model_id;data['formal_interface']='F_ONLY';data['score_transform_version']='PS05_v1.0.2_F_ONLY';data['representation_sha256']=sha(a.representations);data['representation_index_sha256']=sha(a.representation_index) if a.representation_index else None;summary.write_text(json.dumps(data,indent=2)+'\n')
 outputs=(raw,individual,summary,a.output_dir/'AGGREGATION_TRACE.jsonl');manifest={'status':'PASS','files':[{'path':x.name,'bytes':x.stat().st_size,'sha256':sha(x)} for x in outputs]};(a.output_dir/'OUTPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
if __name__=='__main__':main()
