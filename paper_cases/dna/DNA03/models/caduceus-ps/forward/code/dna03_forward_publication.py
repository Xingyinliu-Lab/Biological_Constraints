"""Shared atomic publication for DNA03 v0.2.1 forward adapters."""
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
    keys,v1=primary; keys2,v2=repeat
    exact=keys==keys2 and np.array_equal(v1,v2)
    if not exact: raise RuntimeError('strict repeat mismatch')
    stage=Path(tempfile.mkdtemp(prefix=f'.{output.name}.staging.',dir=output.parent))
    try:
        np.savez_compressed(stage/'EMBEDDINGS.npz',view_ids=np.asarray(keys),embeddings=v1)
        np.savez_compressed(stage/'REPEAT_EMBEDDINGS.npz',view_ids=np.asarray(keys2),embeddings=v2)
        val={'status':'PASS','interface_version':'0.2.3','model_key':model_key,'views':len(keys),'dimension':int(v1.shape[1]),'repeat_bitwise_identical':True,**validation}
        (stage/'FORWARD_VALIDATION.json').write_text(json.dumps(val,indent=2,sort_keys=True)+'\n',encoding='utf-8')
        names=('EMBEDDINGS.npz','REPEAT_EMBEDDINGS.npz','FORWARD_VALIDATION.json')
        files={n:{'sha256':sha256(stage/n),'bytes':(stage/n).stat().st_size} for n in names}
        manifest={'status':'PASS','interface_version':'0.2.3','model_key':model_key,'model_identity':model_identity,
                  'source_sha256':{k:sha256(v) for k,v in source_paths.items()},
                  'code_sha256':{k:sha256(v) for k,v in code_paths.items()},'files':files}
        (stage/'DNA03_FORWARD_MANIFEST.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n',encoding='utf-8')
        (stage/'RUN_COMPLETE').write_text('PASS\n',encoding='utf-8')
        os.replace(stage,output)
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True); raise
