from __future__ import annotations
import argparse,csv,hashlib,json,os,shutil,tempfile,sys
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scoring'))
from dna06_numeric_core import bootstrap_distribution,hierarchical,pair_credit_ledger,parent_scores,site_scores,site_scores_from_pair_ledger

CODE={'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W','CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R','ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
def sh(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''): h.update(b)
 return h.hexdigest()
def tsv(p): return list(csv.DictReader(open(p,encoding='utf-8',newline=''),delimiter='\t'))
def write_tsv(p,rows,fields):
 with open(p,'w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def fasta(p):
 d={};k=None;v=[]
 for line in open(p,encoding='utf-8'):
  line=line.strip()
  if line.startswith('>'):
   if k:d[k]=''.join(v)
   k=line[1:].split()[0];v=[]
  else:v.append(line.upper())
 if k:d[k]=''.join(v)
 return d
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--refset-dir',required=True);ap.add_argument('--output-dir',required=True);a=ap.parse_args();r=Path(a.refset_dir)
 sites=tsv(r/'codon_site.tsv');vars=tsv(r/'variant_endpoint.tsv');members=tsv(r/'site_variant_membership.tsv');parents=tsv(r/'parent_cds.tsv');objs=tsv(r/'sequence_object.tsv');seq=fasta(r/'records.fasta')
 if (len(sites),len(vars),len(parents),len(set(x['dependency_cluster_id'] for x in parents)))!=(3920,11506,560,554):raise ValueError('formal population mismatch')
 sm={x['site_id']:x for x in sites};pm={x['parent_cds_id']:x for x in parents};om={x['object_id']:x for x in objs};by_sha={x['sequence_sha256']:x for x in objs};mm={x['variant_id']:x for x in members};rows=[]
 if set(mm)!={x['variant_id'] for x in vars}:raise ValueError('incomplete membership')
 for v in vars:
  s=sm[v['site_id']];p=pm[v['parent_cds_id']];wt=by_sha.get(p['cds_sequence_sha256'])
  if wt is None:raise ValueError('WT identity')
  ws=seq[wt['object_id']];ms=seq[v['object_id']];diff=[i for i,(x,y) in enumerate(zip(ws,ms)) if x!=y]
  if len(ws)!=len(ms) or diff!=[int(s['codon_nt_position_0based'])+3*int(s['codon_index_0based'])]:raise ValueError('not declared one-SNV')
  q=3*int(s['codon_index_0based']);wc=ws[q:q+3];mc=ms[q:q+3];derived='SYNONYMOUS' if CODE[wc]==CODE[mc] else 'ORDINARY_MISSENSE';m=mm[v['variant_id']]
  if CODE[wc]=='*' or CODE[mc]=='*' or derived!=v['consequence_label'] or m['site_id']!=v['site_id'] or m['consequence_label']!=derived or m['endpoint_role']!='ALT':raise ValueError('Table 1 truth or membership mismatch')
  rows.append({'variant_id':v['variant_id'],'site_id':v['site_id'],'parent_cds_id':v['parent_cds_id'],'dependency_cluster_id':p['dependency_cluster_id'],'label':derived,'decision_score':1 if derived=='ORDINARY_MISSENSE' else 0})
 scoring_rows=[dict(x,score=x['decision_score']) for x in rows];pairs=pair_credit_ledger(scoring_rows,0.0);sr=site_scores_from_pair_ledger(pairs)
 if sr!=site_scores(scoring_rows,0.0):raise AssertionError('pair-to-site aggregation mismatch')
 pr=parent_scores(sr);F,cr=hierarchical(sr)
 if F!=1 or any(x!=1 for x in cr.values()):raise ValueError('oracle not exact')
 out=Path(a.output_dir);out.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix=out.name+'.stage.',dir=str(out.parent)))
 try:
  result={'constraint_id':'DNA06','method_version':'0.1.3','comparator_id':'DNA06_STANDARD_GENETIC_CODE_TABLE1','comparator_class':'ANALYTIC_BIOLOGICAL_KNOWLEDGE_REFERENCE','formal_result':True,'score_status':'REPORTED','formal_interface':'F_ONLY','score_transform_version':'DNA06_v0.1.3_F_ONLY','F':1.0,'F_95_CI':[1.0,1.0],'sites':3920,'variants':11506,'dependency_clusters':554,'prediction_source':'independent raw WT/ALT sequence translation; stored consequence and amino-acid fields not used for prediction'}
  write_tsv(stage/'INDIVIDUAL_SCORES.tsv',rows,list(rows[0]))
  write_tsv(stage/'PAIR_CREDIT_LEDGER.tsv',pairs,list(pairs[0]))
  write_tsv(stage/'SITE_SCORE_LEDGER.tsv',sr,list(sr[0]))
  write_tsv(stage/'PARENT_SCORE_LEDGER.tsv',pr,list(pr[0]))
  cluster_rows=[{'dependency_cluster_id':k,'cluster_F':format(v,'.17g')} for k,v in sorted(cr.items())]
  write_tsv(stage/'CLUSTER_SCORE_LEDGER.tsv',cluster_rows,list(cluster_rows[0]))
  (stage/'FORMAL_RESULT.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');np.save(stage/'F_BOOTSTRAP.npy',bootstrap_distribution(cr),allow_pickle=False);(stage/'INPUT_BINDING.json').write_text(json.dumps({x:sh(r/x) for x in ['codon_site.tsv','variant_endpoint.tsv','parent_cds.tsv','sequence_object.tsv','records.fasta']},indent=2,sort_keys=True)+'\n')
  files={p.name:{'bytes':p.stat().st_size,'sha256':sh(p)} for p in stage.iterdir()}
  (stage/'OUTPUT_MANIFEST.json').write_text(json.dumps(files,indent=2,sort_keys=True)+'\n');(stage/'RUN_COMPLETE').write_text('PASS\n')
  if out.exists():raise FileExistsError(out)
  os.replace(stage,out)
 except Exception:shutil.rmtree(stage,ignore_errors=True);raise
if __name__=='__main__':main()
