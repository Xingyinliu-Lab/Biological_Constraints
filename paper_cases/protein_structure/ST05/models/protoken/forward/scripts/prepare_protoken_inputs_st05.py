#!/usr/bin/env python3
import argparse, hashlib, json, shutil
from pathlib import Path
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--input-dir',type=Path,required=True); ap.add_argument('--output-dir',type=Path,required=True); a=ap.parse_args(); a.output_dir.mkdir(parents=True,exist_ok=True); rows=[]
    for path in sorted(a.input_dir.glob('*.pdb')):
        lines=path.read_text(errors='strict').splitlines(keepends=True); affected=any(x.startswith(('ATOM  ','HETATM')) and len(x)>=27 and x[26]!=' ' for x in lines); out=a.output_dir/path.name
        if not affected: shutil.copyfile(path,out); rows.append({'object_id':path.stem,'input_sha256':sha(path),'output_sha256':sha(out),'normalized':False}); continue
        residue_map={}; chain_next={}; new=[]
        for line in lines:
            if line.startswith(('ATOM  ','HETATM')) and len(line)>=27:
                chain=line[21]; old=(chain,line[22:26],line[26])
                if old not in residue_map: chain_next[chain]=chain_next.get(chain,0)+1; residue_map[old]=chain_next[chain]
                line=line[:22]+f'{residue_map[old]:4d}'+' '+line[27:]
            new.append(line)
        out.write_text(''.join(new)); old=[x for x in lines if x.startswith(('ATOM  ','HETATM'))]; nxt=[x for x in new if x.startswith(('ATOM  ','HETATM'))]
        if len(old)!=len(nxt) or any(x[:22]!=y[:22] or x[27:]!=y[27:] for x,y in zip(old,nxt)): raise ValueError('normalization changed protected fields')
        rows.append({'object_id':path.stem,'input_sha256':sha(path),'output_sha256':sha(out),'normalized':True,'atom_records':len(old),'residues':len(residue_map)})
    rep={'status':'PASS','policy':'ST02 frozen insertion-code normalization; only resSeq/iCode columns change','input_file_count':len(rows),'normalized_file_count':sum(x['normalized'] for x in rows),'coordinate_atom_identity_check':'PASS','records':rows}; (a.output_dir/'INPUT_NORMALIZATION_AUDIT.json').write_text(json.dumps(rep,indent=2,sort_keys=True)+'\n'); print(json.dumps({k:rep[k] for k in ('status','input_file_count','normalized_file_count')},sort_keys=True))
if __name__=='__main__': main()
