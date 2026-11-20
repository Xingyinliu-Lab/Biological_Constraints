import argparse, csv, hashlib, json, re, subprocess
from pathlib import Path

ap=argparse.ArgumentParser();ap.add_argument('--panel',type=Path,required=True);ap.add_argument('--rnapaln',default='RNApaln');ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args()
PANEL=a.panel;OUT=a.output_dir;OUT.mkdir(parents=True,exist_ok=True)
BASE=[a.rnapaln,'--mode=p']
REQUIRED=['--printAlignment']

def run_pair(a,b,print_alignment=True):
    command=BASE + (REQUIRED if print_alignment else [])
    inp=f'>a\n{a}\n>b\n{b}\n'
    proc=subprocess.run(command,input=inp,text=True,capture_output=True,timeout=180)
    numbers=[float(x) for x in re.findall(r'(?m)^(-?\d+(?:\.\d+)?)\s*$',proc.stdout)]
    if proc.returncode or not numbers:
        raise RuntimeError(f'RNApaln failure: rc={proc.returncode}; stderr={proc.stderr[-300:]}')
    return numbers[0],proc.stdout

# Regression oracle: the executable's score is materialized only through backtracking.
toy='GGGAAACCC'
sentinel,_=run_pair(toy,toy,False)
finite,_=run_pair(toy,toy,True)
if sentinel != -9999.0 or not (finite != -9999.0 and finite == finite):
    raise RuntimeError(f'RNApaln printAlignment semantic oracle failed: noflag={sentinel}, flag={finite}')
oracle={'without_printAlignment_score':sentinel,'with_printAlignment_score':finite,'status':'PASS'}
(OUT/'RNApaln_PRINTALIGNMENT_REGRESSION.json').write_text(json.dumps(oracle,indent=2)+'\n')

seq={}; ident=None; buf=[]
for line in (PANEL/'scoring_input_v01/RNA09_PANEL_ENDPOINTS.fasta').read_text().splitlines():
    if line.startswith('>'):
        if ident: seq[ident]=''.join(buf).upper().replace('T','U')
        ident=line[1:].strip();buf=[]
    else: buf.append(line.strip())
if ident: seq[ident]=''.join(buf).upper().replace('T','U')
obj={r['object_id']:r['sequence_sha256'] for r in csv.DictReader(open(PANEL/'scoring_input_v01/RNA09_ENDPOINT_OBJECT_MAP.tsv'),delimiter='\t')}
contrasts=list(csv.DictReader(open(PANEL/'FORMAL_CONTRAST_LEDGER.tsv'),delimiter='\t'))
# Endpoint roles and truth labels are discarded before method invocation.
pairs={(obj[r['anchor_id']],obj[r[field]]) for r in contrasts for field in ('positive_partner_id','negative_partner_id')}
needs={(x,x) for x in seq}|pairs|{(b,a) for a,b in pairs}
scores={}
for i,(a,b) in enumerate(sorted(needs)):
    score,_=run_pair(seq[a],seq[b],True)
    if score == -9999.0 or not (score == score): raise RuntimeError(f'non-materialized score for {a},{b}')
    scores[a+'\t'+b]=score
    print(json.dumps({'done':i+1,'total':len(needs)}),flush=True)
parameters={'comparator_id':'CMP-RNA09-RNAPALN-RAWSEQ-v0.1.0','rnapaln_version':subprocess.check_output([a.rnapaln,'--version'],text=True,stderr=subprocess.STDOUT).strip(),'base_parameters':BASE,'required_execution_flag':'--printAlignment','parser':'first full-line numeric score; alignment text ignored','pair_count':len(scores),'sequence_count':len(seq),'score_sha256':hashlib.sha256(json.dumps(scores,sort_keys=True).encode()).hexdigest()}
(OUT/'PRIMARY_PARAMETERS.json').write_text(json.dumps(parameters,indent=2)+'\n')
(OUT/'RAW_DIRECTIONAL_SCORES.json').write_text(json.dumps(scores,sort_keys=True,indent=2)+'\n')
