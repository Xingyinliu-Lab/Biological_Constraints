#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, math, os, re, subprocess, tempfile
from pathlib import Path

NEEDLE_PARAMS = {"matrix":"EBLOSUM62","gapopen":"10.0","gapextend":"0.5","endweight":"Y","endopen":"10.0","endextend":"0.5"}

def sha256_file(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024), b''): h.update(chunk)
    return h.hexdigest()

def open_text(path: Path): return gzip.open(path,'rt') if path.suffix=='.gz' else path.open('rt')

def read_fasta(path: Path) -> dict[str,str]:
    seqs={}; cur=None; buf=[]
    with open_text(path) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if cur is not None:
                    seq=''.join(buf).upper();
                    if cur in seqs: raise ValueError(f'duplicate fasta id: {cur}')
                    if not re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWYBXZJUO]+', seq): raise ValueError(f'non-protein characters for {cur}')
                    seqs[cur]=seq
                cur=line[1:].split()[0]; buf=[]
            else: buf.append(line)
    if cur is not None:
        seq=''.join(buf).upper()
        if cur in seqs: raise ValueError(f'duplicate fasta id: {cur}')
        if not re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWYBXZJUO]+', seq): raise ValueError(f'non-protein characters for {cur}')
        seqs[cur]=seq
    return seqs

def write_one_fasta(path: Path, ident: str, seq: str) -> None:
    with path.open('wt') as f:
        f.write(f'>{ident}\n')
        for i in range(0,len(seq),80): f.write(seq[i:i+80]+'\n')

def parse_needle_score(text: str) -> float:
    for line in text.splitlines():
        if 'Score:' in line:
            m=re.search(r'Score:\s*([-+0-9.eE]+)', line)
            if m: return float(m.group(1))
    # EMBOSS -aformat score may emit: "a b 20 (116.0)"
    m=re.search(r'\(([-+0-9.eE]+)\)', text)
    if m: return float(m.group(1))
    nums=re.findall(r'[-+]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[eE][-+]?\d+)?', text)
    if len(nums)==1: return float(nums[0])
    raise ValueError(f'cannot parse EMBOSS needle score from: {text[:200]!r}')

def needle_score(needle: Path, seq_a: str, seq_b: str, tmpdir: Path) -> float:
    fa=tmpdir/'a.fa'; fb=tmpdir/'b.fa'; out=tmpdir/'out.score'
    write_one_fasta(fa,'a',seq_a); write_one_fasta(fb,'b',seq_b)
    cmd=[str(needle),'-asequence',str(fa),'-bsequence',str(fb),'-datafile',NEEDLE_PARAMS['matrix'],
         '-gapopen',NEEDLE_PARAMS['gapopen'],'-gapextend',NEEDLE_PARAMS['gapextend'],
         '-endweight',NEEDLE_PARAMS['endweight'],'-endopen',NEEDLE_PARAMS['endopen'],'-endextend',NEEDLE_PARAMS['endextend'],
         '-aformat','score','-outfile',str(out),'-auto']
    env=os.environ.copy(); env['PATH']=str(needle.parent)+os.pathsep+env.get('PATH','')
    env.setdefault('EMBOSS_ACDROOT', str(needle.parent.parent/'share/EMBOSS/acd'))
    env.setdefault('EMBOSS_DATA', str(needle.parent.parent/'share/EMBOSS/data'))
    env.setdefault('PLPLOT_LIB', str(needle.parent.parent/'share/EMBOSS'))
    cp=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env)
    if cp.returncode!=0: raise RuntimeError(f'needle failed rc={cp.returncode}: {cp.stderr[:500]}')
    return parse_needle_score(out.read_text())

def norm_score(pair_score: float, self_a: float, self_b: float) -> float:
    den=self_a+self_b
    if den==0 or not math.isfinite(den): raise ValueError('invalid self-score denominator')
    return 2.0*pair_score/den

def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument('--cells', required=True, type=Path)
    ap.add_argument('--fasta', required=True, type=Path)
    ap.add_argument('--needle', required=True, type=Path)
    ap.add_argument('--outdir', required=True, type=Path)
    ap.add_argument('--limit', type=int, default=None, help='Preflight row limit; omit for full score ledger')
    args=ap.parse_args(); args.outdir.mkdir(parents=True, exist_ok=True)
    seqs=read_fasta(args.fasta)
    with open_text(args.cells) as f: rows=list(csv.DictReader(f, delimiter='\t'))
    if not rows: raise ValueError('empty cells ledger')
    required={'anchor','positive','negative'}
    if not required.issubset(rows[0].keys()): raise ValueError(f'missing columns {required}')
    if args.limit is not None: rows=rows[:args.limit]
    ids=sorted({r['anchor'] for r in rows}|{r['positive'] for r in rows}|{r['negative'] for r in rows})
    missing=[x for x in ids if x not in seqs]
    if missing: raise ValueError(f'{len(missing)} ids missing from fasta; first={missing[:5]}')
    pair_keys=set()
    for r in rows:
        a,p,n=r['anchor'],r['positive'],r['negative']
        pair_keys.add(tuple(sorted((a,p)))); pair_keys.add(tuple(sorted((a,n))))
    with tempfile.TemporaryDirectory(prefix='ps01_nw_') as td:
        tmp=Path(td)
        self_scores={i:needle_score(args.needle, seqs[i], seqs[i], tmp) for i in ids}
        pair_scores={k:needle_score(args.needle, seqs[k[0]], seqs[k[1]], tmp) for k in sorted(pair_keys)}
    pair_path=args.outdir/'PAIR_SCORE_CACHE.tsv.gz'
    with gzip.open(pair_path,'wt',newline='') as f:
        w=csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['id_a','id_b','needle_score','self_a','self_b','normalized_similarity'])
        for a,b in sorted(pair_scores):
            s=pair_scores[(a,b)]; w.writerow([a,b,repr(s),repr(self_scores[a]),repr(self_scores[b]),repr(norm_score(s,self_scores[a],self_scores[b]))])
    ledger_path=args.outdir/'CONTRAST_SCORE_LEDGER.tsv.gz'; wins=ties=losses=0
    with gzip.open(ledger_path,'wt',newline='') as f:
        cols=list(rows[0].keys())+['S_anchor_positive','S_anchor_negative','delta','F_credit']
        w=csv.DictWriter(f, fieldnames=cols, delimiter='\t', lineterminator='\n'); w.writeheader()
        for r in rows:
            a,p,n=r['anchor'],r['positive'],r['negative']
            sp=norm_score(pair_scores[tuple(sorted((a,p)))], self_scores[a], self_scores[p])
            sn=norm_score(pair_scores[tuple(sorted((a,n)))], self_scores[a], self_scores[n])
            delta=sp-sn
            if delta>0: F_credit=1.0; wins+=1
            elif delta<0: F_credit=-1.0; losses+=1
            else: F_credit=0.0; ties+=1
            out=dict(r); out.update(S_anchor_positive=repr(sp), S_anchor_negative=repr(sn), delta=repr(delta), F_credit=repr(F_credit)); w.writerow(out)
    summary={'artifact_type':'PS01_NW_BLOSUM62_COMPARATOR_SCORE_LEDGER_SUMMARY','status':'PASS',
             'formal_result': False if args.limit is not None else 'SCORE_LEDGER_ONLY_NOT_AGGREGATED',
             'cells_path':str(args.cells),'cells_sha256':sha256_file(args.cells),'fasta_path':str(args.fasta),'fasta_sha256':sha256_file(args.fasta),
             'needle_path':str(args.needle),'needle_params':NEEDLE_PARAMS,'row_count':len(rows),'unique_sequence_count':len(ids),
             'pair_score_count':len(pair_scores),'wins':wins,'ties':ties,'losses':losses,
             'pair_score_cache_sha256':sha256_file(pair_path),'contrast_score_ledger_sha256':sha256_file(ledger_path)}
    (args.outdir/'COMPARATOR_LEDGER_SUMMARY.json').write_text(json.dumps(summary, indent=2, sort_keys=True)+'\n')
if __name__=='__main__': main()
