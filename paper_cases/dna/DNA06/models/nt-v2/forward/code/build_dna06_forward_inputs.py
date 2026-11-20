from __future__ import annotations
import argparse,csv,hashlib
from pathlib import Path
def tsv(p):return list(csv.DictReader(open(p,encoding='utf-8',newline=''),delimiter='\t'))
def fasta(p):
 d={};k=None
 for line in open(p,encoding='utf-8'):
  line=line.strip()
  if line.startswith('>'):k=line[1:].split()[0];d[k]=''
  else:d[k]+=line.upper()
 return d
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--refset-dir',required=True);ap.add_argument('--output-dir',required=True);a=ap.parse_args();r=Path(a.refset_dir);o=Path(a.output_dir);o.mkdir(parents=True,exist_ok=False)
 sites={x['site_id']:x for x in tsv(r/'codon_site.tsv')};vars=tsv(r/'variant_endpoint.tsv');parents={x['parent_cds_id']:x for x in tsv(r/'parent_cds.tsv')};objs=tsv(r/'sequence_object.tsv');seq=fasta(r/'records.fasta');bysha={x['sequence_sha256']:x['object_id'] for x in objs};rows=[];outseq={}
 for sid,s in sites.items():
  pid=s['parent_cds_id'];oid=bysha[parents[pid]['cds_sequence_sha256']];rid=f'{sid}|WT|{oid}';outseq[rid]=seq[oid];start=int(s['cds_offset_start_0based']);rows.append({'view_id':rid,'target_length_nt':len(seq[oid]),'view_sequence_sha256':hashlib.sha256(seq[oid].encode()).hexdigest(),'core_offset_start_0based':start,'core_offset_end_exclusive':start+3})
 for v in vars:
  s=sites[v['site_id']];rid=f"{v['site_id']}|ALT|{v['variant_id']}|{v['object_id']}";outseq[rid]=seq[v['object_id']];start=int(s['cds_offset_start_0based']);rows.append({'view_id':rid,'target_length_nt':len(outseq[rid]),'view_sequence_sha256':hashlib.sha256(outseq[rid].encode()).hexdigest(),'core_offset_start_0based':start,'core_offset_end_exclusive':start+3})
 if len(rows)!=15426 or len(outseq)!=15426:raise ValueError('formal endpoint plan mismatch')
 with open(o/'DNA06_FORWARD_MEMBERSHIP.tsv','w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
 with open(o/'DNA06_FORWARD_INPUT.fasta','w',encoding='utf-8') as f:
  for rid in sorted(outseq):f.write(f'>{rid}\n{outseq[rid]}\n')
if __name__=='__main__':main()
