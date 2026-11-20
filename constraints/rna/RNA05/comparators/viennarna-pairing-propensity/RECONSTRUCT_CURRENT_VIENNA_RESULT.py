#!/usr/bin/env python3
"""Reconstruct the RNA05 Vienna comparator F from released raw measurements."""
import argparse, json, hashlib, math
from pathlib import Path

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--score-ledger', required=True)
    ap.add_argument('--compatibility-audit', required=True)
    ap.add_argument('--result', required=True)
    ap.add_argument('--out', required=True)
    a=ap.parse_args()
    score=json.loads(Path(a.score_ledger).read_text())
    audit=json.loads(Path(a.compatibility_audit).read_text())
    result=json.loads(Path(a.result).read_text())
    rows=score['records']
    clusters={}
    for r in rows:
        m=r['member']; clusters.setdefault(m['dependency_cluster_id'],{}).setdefault(m['parent_id'],{}).setdefault(m['matched_block_id'],[]).append(float(r['F_credit']))
    cluster_scores={}
    for cid,parents in clusters.items():
        pvals=[]
        for _pid,blocks in parents.items():
            pvals.append(sum(sum(v)/len(v) for v in blocks.values())/len(blocks))
        cluster_scores[cid]=sum(pvals)/len(pvals)
    F=sum(cluster_scores.values())/len(cluster_scores)
    out={'F':F,'clusters':len(cluster_scores),'leaves':len(rows),'compatible_blocks':audit['compatible_blocks'],'oos_blocks':audit['oos_blocks'],'score_ledger_sha256':sha(a.score_ledger),'compatibility_audit_sha256':sha(a.compatibility_audit),'reported_result_sha256':sha(a.result),'reported_F':result['result']['F'],'F_matches_reported':abs(F-result['result']['F'])<1e-12,'cluster_F':cluster_scores}
    Path(a.out).write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
if __name__=='__main__': main()
