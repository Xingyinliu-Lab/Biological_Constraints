#!/usr/bin/env python3
"""Bind a frozen PS04 direct-alignment comparator to the frozen formal scorer."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys
from pathlib import Path

def sha(p: Path) -> str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def emboss_version(tool: Path) -> str:
    proc=subprocess.run([str(tool),'-version'],capture_output=True,text=True,check=True)
    version=(proc.stdout+'\n'+proc.stderr).strip()
    if '6.6.0' not in version:
        raise SystemExit(f'unsupported EMBOSS version: {version}')
    return version

def main():
    p=argparse.ArgumentParser()
    constraint_root = Path(__file__).resolve().parents[2]
    refset = constraint_root / 'refset'
    scoring = constraint_root / 'scoring'
    p.add_argument('--relations', type=Path, default=refset / 'FORMAL_RELATION_LEDGER.tsv.gz')
    p.add_argument('--fasta', type=Path, default=refset / 'FORMAL_ENDPOINTS.faa')
    p.add_argument('--uncertainty-membership', type=Path, default=refset / 'UNCERTAINTY_ANCHOR_SOURCE_UNIT_MEMBERSHIP.tsv.gz')
    p.add_argument('--endpoint-reuse-ledger', type=Path, default=refset / 'ENDPOINT_REUSE_LEDGER.tsv.gz')
    p.add_argument('--formal-binding', type=Path, default=scoring / 'FORMAL_BINDING.json')
    p.add_argument('--comparator', type=Path, default=Path(__file__).with_name('ps04_direct_alignment_comparator.py'))
    p.add_argument('--formal-scorer', type=Path, default=scoring / 'score_ps04_formal.py')
    p.add_argument('--tool', type=Path, required=True)
    p.add_argument('--outdir',type=Path,required=True); p.add_argument('--workers',type=int,default=64); args=p.parse_args()
    mode='NW'
    version=emboss_version(args.tool)
    args.outdir.mkdir(parents=True,exist_ok=False); measurement=args.outdir/'measurement'; result=args.outdir/'FORMAL_RESULT.json'
    subprocess.run([sys.executable,str(args.comparator),'--relations',str(args.relations),'--fasta',str(args.fasta),'--tool',str(args.tool),'--mode',mode,'--outdir',str(measurement),'--workers',str(args.workers)],check=True)
    subprocess.run([sys.executable,str(args.formal_scorer),'--binding',str(args.formal_binding),'--formal-relations',str(args.relations),'--leaf-ledger',str(measurement/'FORMAL_LEAF_SCORES.tsv.gz'),'--uncertainty-membership',str(args.uncertainty_membership),'--endpoint-reuse-ledger',str(args.endpoint_reuse_ledger),'--output-dir',str(args.outdir/'scorer')],check=True)
    scorer_result=args.outdir/'scorer'/'FORMAL_RESULT.json'; result.write_bytes(scorer_result.read_bytes())
    (args.outdir/'COMPARATOR_EXECUTION.json').write_text(json.dumps({'kind':'PS04_NW_BLOSUM62_COMPARATOR_EXECUTION','status':'FORMAL_EXECUTION_COMPLETE','method':'Needleman-Wunsch with EBLOSUM62','mode':mode,'measurement_ledger_sha256':sha(measurement/'FORMAL_LEAF_SCORES.tsv.gz'),'result_sha256':sha(result),'execution_driver_sha256':sha(Path(__file__).resolve()),'tool_sha256':sha(args.tool),'tool_version':version,'comparator_code_sha256':sha(args.comparator),'formal_scorer_sha256':sha(args.formal_scorer),'formal_binding_sha256':sha(args.formal_binding),'workers':args.workers},indent=2)+'\n')

if __name__=='__main__': main()
