"""Shared atomic publication for the formal DNA06 codon-local interface."""
from __future__ import annotations
import hashlib, json, os, shutil, tempfile
from pathlib import Path
import numpy as np


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()


def publish_forward(output: Path, *, model_key: str, model_identity: dict,
                    source_paths: dict[str,Path], primary, repeat,
                    validation: dict, code_paths: dict[str,Path]):
    if output.exists(): raise ValueError('output path exists')
    keys,embeddings=primary; keys2,repeat_embeddings=repeat
    exact=keys==keys2 and np.array_equal(embeddings,repeat_embeddings)
    if not exact: raise RuntimeError('strict repeat mismatch')
    output.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix=f'.{output.name}.staging.',dir=output.parent))
    try:
        np.savez_compressed(stage/'EMBEDDINGS.npz',row_ids=np.asarray(keys),embeddings=embeddings)
        np.savez_compressed(stage/'REPEAT_EMBEDDINGS.npz',row_ids=np.asarray(keys2),embeddings=repeat_embeddings)
        val={'status':'PASS','interface_version':'0.1.2','representation':'codon-local','model_key':model_key,'rows':len(keys),'dimension':int(embeddings.shape[1]),'repeat_bitwise_identical':True,**validation}
        (stage/'FORWARD_VALIDATION.json').write_text(json.dumps(val,indent=2,sort_keys=True)+'\n',encoding='utf-8')
        names=('EMBEDDINGS.npz','REPEAT_EMBEDDINGS.npz','FORWARD_VALIDATION.json')
        files={n:{'sha256':sha256(stage/n),'bytes':(stage/n).stat().st_size} for n in names}
        manifest={'status':'PASS','interface_version':'0.1.2','representation':'codon-local','model_key':model_key,'model_identity':model_identity,
                  'source_sha256':{k:sha256(v) for k,v in source_paths.items()},
                  'code_sha256':{k:sha256(v) for k,v in code_paths.items()},'files':files}
        (stage/'DNA06_FORWARD_MANIFEST.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n',encoding='utf-8')
        (stage/'RUN_COMPLETE').write_text('PASS\n',encoding='utf-8')
        os.replace(stage,output)
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True); raise
