#!/usr/bin/env python3
import argparse, importlib.util, json, math
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--scorer',type=Path,required=True);p.add_argument('--contract-root',type=Path,required=True);p.add_argument('--leaf',type=Path,required=True);p.add_argument('--score-dir',type=Path,required=True);p.add_argument('--delta',type=float,default=0.0);a=p.parse_args()
spec=importlib.util.spec_from_file_location('ps04score',a.scorer);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
binding=mod.validate_contract(a.contract_root); rows=mod.load_rows(a.leaf,a.delta); trace,raw=mod.score(rows,mod.digest(a.leaf))
summary=json.loads((a.score_dir/'SCORE_SUMMARY.json').read_text())
if not math.isclose(raw,summary['raw']['value'],abs_tol=1e-12):raise ValueError(f'RAW_SCORE_MISMATCH {raw} {summary["raw"]["value"]}')
anchors=[dict(x,meta=None) for x in trace['ANCHOR']]
(a.score_dir/'AGGREGATION_TRACE.jsonl').write_text('\n'.join(mod.canonical_json(x) for x in anchors)+'\n')
(a.score_dir/'ANCHOR_TRACE_REPORT.json').write_text(json.dumps({'kind':'PS04_ANCHOR_TRACE_FOR_FROZEN_DEPENDENCY_BOOTSTRAP','leaf_count':len(rows),'anchor_count':len(anchors),'raw_R':raw,'binding':binding},indent=2)+'\n')
print(json.dumps({'status':'PASS','leaf_count':len(rows),'anchors':len(anchors),'raw_R':raw}))