from __future__ import annotations
import argparse,csv,hashlib,json,os,shutil,sys,tempfile
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scoring'))
from dna06_numeric_core import bootstrap_distribution,fold_map,standardize_fit,ridge_logistic_fit,predict,inner_select,site_scores,hierarchical,pair_credit_ledger,site_scores_from_pair_ledger,parent_scores

BASES='ACGT'; CODONS=[a+b+c for a in BASES for b in BASES for c in BASES]; CI={x:i for i,x in enumerate(CODONS)}; BI={x:i for i,x in enumerate(BASES)}
def tsv(p):return list(csv.DictReader(open(p,encoding='utf-8',newline=''),delimiter='\t'))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def write_tsv(p,rows,fields):
 with open(p,'w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def feature(w,a,pos,ref,alt):
 x=np.zeros(139);x[CI[w]]=1;x[64+CI[a]]=1;x[128+pos]=1;x[131+BI[ref]]=1;x[135+BI[alt]]=1;return x
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--refset-dir',required=True);ap.add_argument('--output-dir',required=True);a=ap.parse_args();r=Path(a.refset_dir)
 sites={x['site_id']:x for x in tsv(r/'codon_site.tsv')};parents={x['parent_cds_id']:x for x in tsv(r/'parent_cds.tsv')};vars=tsv(r/'variant_endpoint.tsv')
 if len(sites)!=3920 or len(vars)!=11506 or len(parents)!=560 or len(set(x['dependency_cluster_id'] for x in parents.values()))!=554:raise ValueError('population mismatch')
 X=[];Y=[];C=[];S=[];P=[]
 for v in vars:
  s=sites[v['site_id']]
  if v['consequence_label'] not in {'SYNONYMOUS','ORDINARY_MISSENSE'}:raise ValueError('invalid consequence label')
  X.append(feature(s['wt_codon'],v['variant_codon'],int(s['codon_nt_position_0based']),v['ref_allele'],v['alt_allele']));Y.append(v['consequence_label']=='ORDINARY_MISSENSE');C.append(parents[v['parent_cds_id']]['dependency_cluster_id']);S.append(v['site_id']);P.append(v['parent_cds_id'])
 X=np.asarray(X);Y=np.asarray(Y,dtype=int);C=np.asarray(C);S=np.asarray(S);P=np.asarray(P);fm=fold_map(C,5,'DNA06_OUTER_V010');pred=np.empty(len(Y));sel={}
 for f in range(5):
  tr=np.array([fm[c]!=f for c in C]);va=~tr;lam,_=inner_select(X[tr],Y[tr],C[tr],S[tr],P[tr]);mu,sd=standardize_fit(X[tr]);b=ridge_logistic_fit((X[tr]-mu)/sd,Y[tr],lam);pred[va]=predict((X[va]-mu)/sd,b);sel[str(f)]=lam
 if not np.isfinite(pred).all():raise FloatingPointError('non-finite R_input prediction')
 rows=[{'variant_id':vars[i]['variant_id'],'site_id':S[i],'parent_cds_id':P[i],'dependency_cluster_id':C[i],'outer_fold':fm[C[i]],'label':'ORDINARY_MISSENSE' if Y[i] else 'SYNONYMOUS','decision_score':format(pred[i],'.17g')} for i in range(len(Y))];scoring_rows=[dict(x,score=x['decision_score']) for x in rows];pairs=pair_credit_ledger(scoring_rows);sr=site_scores_from_pair_ledger(pairs)
 if sr!=site_scores(scoring_rows):raise AssertionError('pair-to-site aggregation mismatch')
 pr=parent_scores(sr);F,cr=hierarchical(sr);boot=bootstrap_distribution(cr);ci=[float(np.quantile(boot,.025)),float(np.quantile(boot,.975))]
 out={'constraint_id':'DNA06','companion_id':'DNA06_RAW_CODON_RIDGE','status':'REPORTED','formal_interface':'F_ONLY','score_transform_version':'DNA06_v0.1.3_F_ONLY','F':F,'F_95_CI':ci,'selected_lambda_by_outer_fold':sel,'feature_dimension':139,'forbidden_label_or_amino_acid_features':True}
 output=Path(a.output_dir)
 if output.exists():raise FileExistsError(output)
 output.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix=output.name+'.stage.',dir=str(output.parent)))
 try:
  write_tsv(stage/'INDIVIDUAL_SCORES.tsv',rows,list(rows[0]));write_tsv(stage/'PAIR_CREDIT_LEDGER.tsv',pairs,list(pairs[0]));write_tsv(stage/'SITE_SCORE_LEDGER.tsv',sr,list(sr[0]));write_tsv(stage/'PARENT_SCORE_LEDGER.tsv',pr,list(pr[0]))
  cluster_rows=[{'dependency_cluster_id':k,'cluster_F':format(v,'.17g')} for k,v in sorted(cr.items())];write_tsv(stage/'CLUSTER_SCORE_LEDGER.tsv',cluster_rows,list(cluster_rows[0]));np.save(stage/'F_BOOTSTRAP.npy',boot,allow_pickle=False)
  (stage/'FORMAL_RESULT.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
  inputs={x:{'bytes':(r/x).stat().st_size,'sha256':sha(r/x)} for x in ['codon_site.tsv','variant_endpoint.tsv','parent_cds.tsv']};(stage/'INPUT_BINDING.json').write_text(json.dumps(inputs,indent=2,sort_keys=True)+'\n')
  files={p.name:{'bytes':p.stat().st_size,'sha256':sha(p)} for p in stage.iterdir()};(stage/'OUTPUT_MANIFEST.json').write_text(json.dumps(files,indent=2,sort_keys=True)+'\n');(stage/'RUN_COMPLETE').write_text('PASS\n')
  os.replace(stage,output)
 except Exception:shutil.rmtree(stage,ignore_errors=True);raise
if __name__=='__main__':main()
