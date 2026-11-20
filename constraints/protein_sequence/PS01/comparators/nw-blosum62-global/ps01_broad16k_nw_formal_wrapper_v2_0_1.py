"""Package-ready formal PS01 Broad16K NW comparator wrapper.
Execution-equivalent to the retained parallel runner, but resolves the frozen
comparator source next to this file and reads formal a/p/n ledger fields.
"""
import argparse, csv, gzip, hashlib, importlib.util, json, sys, tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
SOURCE = Path(__file__).with_name("ps01_nw_blosum62_comparator.py")
spec = importlib.util.spec_from_file_location("nw", SOURCE)
nw = importlib.util.module_from_spec(spec); spec.loader.exec_module(nw)
SCORING = Path(__file__).resolve().parents[2] / "scoring"
sys.path.insert(0, str(SCORING))
from ps01_scoring_core import evaluate
def sha(p):
 h=hashlib.sha256();
 with p.open("rb") as f:
  for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
 return h.hexdigest()
def read_fasta(p):
 d={}; ident=None; seq=[]
 with gzip.open(p,"rt") as f:
  for line in f:
   line=line.strip()
   if line.startswith(">"):
    if ident: d[ident]="".join(seq)
    ident=line[1:].split()[0]; seq=[]
   else: seq.append(line)
  if ident: d[ident]="".join(seq)
 return d
ARGS=None
def one(task):
 k,x=task
 with tempfile.TemporaryDirectory(prefix="ps01_nw_") as directory:
  tmpdir=Path(directory)
  if k=="self": return task,nw.needle_score(ARGS.needle,ARGS.seqs[x[0]],ARGS.seqs[x[0]],tmpdir)
  return task,nw.needle_score(ARGS.needle,ARGS.seqs[x[0]],ARGS.seqs[x[1]],tmpdir)
def main():
 global ARGS
 ap=argparse.ArgumentParser(); ap.add_argument("--contrasts",type=Path,required=True); ap.add_argument("--fasta",type=Path,required=True); ap.add_argument("--needle",type=Path,required=True); ap.add_argument("--outdir",type=Path,required=True); ap.add_argument("--workers",type=int,default=64); ARGS=ap.parse_args(); ARGS.outdir.mkdir(parents=True,exist_ok=True)
 with gzip.open(ARGS.contrasts,"rt",newline="") as f: rows=list(csv.DictReader(f,delimiter="\t"))
 seqs=read_fasta(ARGS.fasta); ids=sorted({r["a"] for r in rows}|{r["p"] for r in rows}|{r["n"] for r in rows}); ARGS.seqs=seqs
 if len(rows)!=1256 or len(ids)!=2316 or any(x not in seqs for x in ids): raise ValueError("unexpected Broad16K formal scope")
 pairs=sorted({tuple(sorted((r["a"],r["p"]))) for r in rows}|{tuple(sorted((r["a"],r["n"]))) for r in rows}); tasks=[("self",(x,)) for x in ids]+[("pair",x) for x in pairs]
 with ThreadPoolExecutor(max_workers=ARGS.workers) as pool: scores=dict(pool.map(one,tasks))
 selfs={x:scores[("self",(x,))] for x in ids}; pairscores={x:scores[("pair",x)] for x in pairs}
 raw_path=ARGS.outdir/"RAW_MEASUREMENTS.tsv.gz"
 with gzip.open(raw_path,"wt",newline="") as f:
  w=csv.writer(f,delimiter="\t",lineterminator="\n"); w.writerow(["id_a","id_b","needle_score","self_a","self_b","normalized_similarity"])
  for a,b in pairs: s=pairscores[(a,b)]; w.writerow([a,b,repr(s),repr(selfs[a]),repr(selfs[b]),repr(nw.norm_score(s,selfs[a],selfs[b]))])
 individual_path=ARGS.outdir/"INDIVIDUAL_SCORES.tsv.gz"
 with gzip.open(individual_path,"wt",newline="") as f:
  fields=["anchor","positive","clan_id","anchor_family","negative_family","negative","anchor_parent","positive_parent","negative_parent","positive_identity_count","positive_denominator","negative_identity_count","negative_denominator","length_ratio_pn","composition_l1_5_pn","same_superkingdom_pn","superkingdom_anchor","superkingdom_positive","superkingdom_negative","S_anchor_positive","S_anchor_negative","delta","F_credit"]
  w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader()
  for r in rows:
   a,p,n=r["a"],r["p"],r["n"]; sp=nw.norm_score(pairscores[tuple(sorted((a,p)))],selfs[a],selfs[p]); sn=nw.norm_score(pairscores[tuple(sorted((a,n)))],selfs[a],selfs[n]); d=sp-sn
   out=dict(r)
   out.update(anchor=a,positive=p,negative=n,
              positive_identity_count=r.get("identity_p_num",""),positive_denominator=r.get("identity_p_den",""),
              negative_identity_count=r.get("identity_n_num",""),negative_denominator=r.get("identity_n_den",""),
              S_anchor_positive=repr(sp),S_anchor_negative=repr(sn),delta=repr(d),F_credit=repr(1.0 if d>0 else 0.0 if d==0 else -1.0))
   w.writerow({key:out.get(key,"") for key in fields})
 result=evaluate(individual_path,return_bootstrap=True); F_bootstrap=result.pop("F_bootstrap")
 result.update({"formal_result":True,"formal_interface":"F_ONLY","score_transform_version":"PS01_v2.0.2_F_ONLY","comparator_id":"nw-blosum62-global","status":"PASS"})
 result_path=ARGS.outdir/"FORMAL_RESULT.json"; result_path.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
 import numpy as np; np.save(ARGS.outdir/"F_BOOTSTRAP.npy",np.asarray(F_bootstrap,dtype=np.float64),allow_pickle=False)
 summary={"status":"PASS","method":"Needleman-Wunsch + EBLOSUM62","formal_interface":"F_ONLY","score_transform_version":"PS01_v2.0.2_F_ONLY","row_count":len(rows),"unique_sequence_count":len(ids),"workers":ARGS.workers,"wins":result["leaf_wins"],"ties":result["leaf_ties"],"losses":result["leaf_losses"],"pair_score_count":len(pairs),"cells_sha256":sha(raw_path),"pair_score_cache_sha256":sha(raw_path),"contrast_score_ledger_sha256":sha(ARGS.contrasts),"fasta_sha256":sha(ARGS.fasta),"individual_scores_sha256":sha(individual_path),"formal_result_sha256":sha(result_path)}
 (ARGS.outdir/"COMPARATOR_LEDGER_SUMMARY.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n",encoding="utf-8")
 outputs=[{"path":path.name,"bytes":path.stat().st_size,"sha256":sha(path)} for path in (raw_path,individual_path,ARGS.outdir/"F_BOOTSTRAP.npy",result_path,ARGS.outdir/"COMPARATOR_LEDGER_SUMMARY.json")]
 (ARGS.outdir/"OUTPUT_MANIFEST.json").write_text(json.dumps({"status":"PASS","files":outputs},indent=2,sort_keys=True)+"\n",encoding="utf-8")
 print(json.dumps(summary,indent=2,sort_keys=True))
if __name__=="__main__": main()
