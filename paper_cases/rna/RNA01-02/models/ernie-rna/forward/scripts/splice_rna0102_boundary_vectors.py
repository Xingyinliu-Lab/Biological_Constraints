import argparse,numpy as np
p=argparse.ArgumentParser();p.add_argument('--base',required=True);p.add_argument('--patch',required=True);p.add_argument('--out',required=True);a=p.parse_args();b=np.load(a.base);q=np.load(a.patch);ids=list(b['sequence_sha256']);idx={x:i for i,x in enumerate(ids)};v=b['vectors'].copy()
for ident,z in zip(q['sequence_sha256'].tolist(),q['vectors']):
 if ident not in idx:raise ValueError('patch id absent from base')
 v[idx[ident]]=z
np.savez(a.out,sequence_sha256=b['sequence_sha256'],vectors=v)
print({'records':len(ids),'patched':len(q['sequence_sha256']),'dimension':v.shape[1]})