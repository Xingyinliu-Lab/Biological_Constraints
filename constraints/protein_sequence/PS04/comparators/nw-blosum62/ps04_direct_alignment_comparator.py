#!/usr/bin/env python3
"""PS04 raw-domain alignment measurement producer; no aggregation occurs here."""
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, math, os, re, sqlite3, subprocess, tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

COMMON = {'matrix':'EBLOSUM62','gapopen':'10.0','gapextend':'0.5'}
NW_EXTRA = {'endweight':'Y','endopen':'10.0','endextend':'0.5'}
PROTEIN = re.compile('[ACDEFGHIKLMNPQRSTVWYBXZJUO]+')

def sha(p: Path) -> str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def open_text(p: Path): return gzip.open(p,'rt',newline='') if p.suffix=='.gz' else p.open(newline='')

def fasta(p: Path) -> dict[str,str]:
    out={}; key=None; parts=[]
    with p.open() as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if key is not None:
                    if key in out: raise ValueError('duplicate FASTA id')
                    out[key]=''.join(parts).upper()
                key=line[1:].split()[0]; parts=[]
            else: parts.append(line)
    if key is not None:
        if key in out: raise ValueError('duplicate FASTA id')
        out[key]=''.join(parts).upper()
    if not out or any(not PROTEIN.fullmatch(x) for x in out.values()): raise ValueError('invalid endpoint FASTA')
    return out

def write_fasta(p: Path, s: str): p.write_text('>x\n'+s+'\n')

def parse_score(text: str) -> float:
    m=re.search(r'\(([-+0-9.eE]+)\)',text) or re.search(r'Score:\s*([-+0-9.eE]+)',text)
    if not m: raise ValueError('unparseable EMBOSS score')
    return float(m.group(1))

def align(binary: Path, mode: str, a: str, b: str) -> float:
    with tempfile.TemporaryDirectory(prefix='ps04_align_') as d:
        d=Path(d); fa,fb,out=d/'a.fa',d/'b.fa',d/'score.txt'; write_fasta(fa,a);write_fasta(fb,b)
        cmd=[str(binary),'-asequence',str(fa),'-bsequence',str(fb),'-datafile',COMMON['matrix'],'-gapopen',COMMON['gapopen'],'-gapextend',COMMON['gapextend']]
        if mode=='NW': cmd += ['-endweight',NW_EXTRA['endweight'],'-endopen',NW_EXTRA['endopen'],'-endextend',NW_EXTRA['endextend']]
        cmd += ['-aformat','score','-outfile',str(out),'-auto']
        env=os.environ.copy(); env['PATH']=str(binary.parent)+os.pathsep+env.get('PATH','')
        env.setdefault('EMBOSS_ACDROOT',str(binary.parent.parent/'share/EMBOSS/acd')); env.setdefault('EMBOSS_DATA',str(binary.parent.parent/'share/EMBOSS/data'))
        q=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env)
        if q.returncode: raise RuntimeError(q.stderr[:500])
        return parse_score(out.read_text())

def norm(raw: float, self_a: float, self_b: float) -> float:
    d=self_a+self_b
    if not math.isfinite(d) or d==0: raise ValueError('invalid self alignment score')
    return 2.0*raw/d

def relation_to_leaf(r: dict, ap: float, an: float) -> dict:
    return {'negative_endpoint_id':r['negative_id'],'negative_hsf_id':r['negative_hsf_id'],'negative_topology_id':r['negative_topology_id'],'positive_endpoint_id':r['positive_id'],'positive_hsf_id':r['positive_hsf_id'],'anchor_id':r['anchor_id'],'anchor_parent_id':r['anchor_parent_proxy'],'anchor_source_unit_id':r['anchor_source_unit_id'],'anchor_hsf_id':r['anchor_hsf_id'],'anchor_topology_id':r['anchor_topology_id'],'anchor_architecture_id':r['anchor_architecture_id'],'positive_score':repr(ap),'negative_score':repr(an)}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--relations',type=Path,required=True); ap.add_argument('--fasta',type=Path,required=True); ap.add_argument('--tool',type=Path,required=True); ap.add_argument('--mode',choices=('NW','SW'),required=True); ap.add_argument('--outdir',type=Path,required=True); ap.add_argument('--workers',type=int,default=1); args=ap.parse_args()
    if args.workers<1: raise ValueError('workers must be positive')
    args.outdir.mkdir(parents=True,exist_ok=False); seq=fasta(args.fasta); db=args.outdir/'pair_index.sqlite'
    con=sqlite3.connect(db); con.execute('CREATE TABLE pair (a TEXT NOT NULL,b TEXT NOT NULL,PRIMARY KEY(a,b))')
    endpoint_ids=set(); n=0
    with open_text(args.relations) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            n+=1; a,p,neg=r['anchor_id'],r['positive_id'],r['negative_id']; endpoint_ids.update((a,p,neg))
            con.execute('INSERT OR IGNORE INTO pair VALUES (?,?)',tuple(sorted((a,p))))
            con.execute('INSERT OR IGNORE INTO pair VALUES (?,?)',tuple(sorted((a,neg))))
    con.commit()
    if set(endpoint_ids)-set(seq): raise ValueError('formal relation endpoint missing from FASTA')
    pairs=list(con.execute('SELECT a,b FROM pair ORDER BY a,b'))
    def job(pair): return pair,align(args.tool,args.mode,seq[pair[0]],seq[pair[1]])
    self_scores={}
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        for identifier,value in ex.map(lambda x:(x,align(args.tool,args.mode,seq[x],seq[x])),sorted(endpoint_ids)): self_scores[identifier]=value
    raw_scores={}
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        for pair,value in ex.map(job,pairs): raw_scores[pair]=value
    cache=args.outdir/'PAIR_SCORE_CACHE.tsv.gz'
    with gzip.open(cache,'wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['id_a','id_b','alignment_score','self_a','self_b','normalized_similarity'],delimiter='\t',lineterminator='\n');w.writeheader()
        for a,b in pairs: w.writerow({'id_a':a,'id_b':b,'alignment_score':repr(raw_scores[(a,b)]),'self_a':repr(self_scores[a]),'self_b':repr(self_scores[b]),'normalized_similarity':repr(norm(raw_scores[(a,b)],self_scores[a],self_scores[b]))})
    leaf=args.outdir/'FORMAL_LEAF_SCORES.tsv.gz'; detail=args.outdir/'CONTRAST_SCORE_LEDGER.tsv.gz'; wins=ties=losses=0
    fields=['negative_endpoint_id','negative_hsf_id','negative_topology_id','positive_endpoint_id','positive_hsf_id','anchor_id','anchor_parent_id','anchor_source_unit_id','anchor_hsf_id','anchor_topology_id','anchor_architecture_id','positive_score','negative_score']
    with open_text(args.relations) as f, gzip.open(leaf,'wt',newline='') as lf, gzip.open(detail,'wt',newline='') as df:
        lw=csv.DictWriter(lf,fieldnames=fields,delimiter='\t',lineterminator='\n'); lw.writeheader(); dw=csv.DictWriter(df,fieldnames=['anchor_id','positive_id','negative_id','S_anchor_positive','S_anchor_negative','delta','F_credit'],delimiter='\t',lineterminator='\n'); dw.writeheader()
        for r in csv.DictReader(f,delimiter='\t'):
            a,p,neg=r['anchor_id'],r['positive_id'],r['negative_id']; sp=norm(raw_scores[tuple(sorted((a,p)))],self_scores[a],self_scores[p]); sn=norm(raw_scores[tuple(sorted((a,neg)))],self_scores[a],self_scores[neg]); delta=sp-sn; f_credit=1. if delta>0 else -1. if delta<0 else 0.
            wins+=f_credit==1.;ties+=f_credit==0.;losses+=f_credit==-1.; lw.writerow(relation_to_leaf(r,sp,sn));dw.writerow({'anchor_id':a,'positive_id':p,'negative_id':neg,'S_anchor_positive':repr(sp),'S_anchor_negative':repr(sn),'delta':repr(delta),'F_credit':repr(f_credit)})
    executed_parameters = dict(COMMON)
    if args.mode == 'NW': executed_parameters.update(NW_EXTRA)
    summary={'artifact_type':'PS04_DIRECT_ALIGNMENT_MEASUREMENT','status':'PASS','formal_result':False,'mode':args.mode,'relation_count':n,'endpoint_count':len(endpoint_ids),'unique_pair_count':len(pairs),'wins':wins,'ties':ties,'losses':losses,'tool_sha256':sha(args.tool),'relation_ledger_sha256':sha(args.relations),'endpoint_fasta_sha256':sha(args.fasta),'pair_score_cache_sha256':sha(cache),'leaf_ledger_sha256':sha(leaf),'detail_ledger_sha256':sha(detail),'fixed_method_prior':'EBLOSUM62 substitution matrix','terminal_gaps_penalized':args.mode=='NW','parameters':executed_parameters}
    (args.outdir/'COMPARATOR_LEDGER_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')

if __name__=='__main__': main()
