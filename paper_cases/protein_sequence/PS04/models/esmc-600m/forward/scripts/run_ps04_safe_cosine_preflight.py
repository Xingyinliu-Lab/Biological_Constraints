#!/usr/bin/env python3
"""Truth-blind repeatability preflight for the frozen PS04 cosine operator."""
import csv, gzip, hashlib, json, math, os
import numpy as np
from pathlib import Path
from safetensors import safe_open
from bc_safe_cosine_total_float32_v1 import safe_cosine_total_float32, OPERATOR_ID

ROOT = Path(os.environ['PS04_PREFLIGHT_ROOT'])
REL = ROOT / 'ps04_assay_refset_release_v1.0.1_L_4K_20261007/FORMAL_RELATION_LEDGER.tsv.gz'
ENDPOINTS = ROOT / 'ps04_assay_refset_release_v1.0.1_L_4K_20261007/FORMAL_ENDPOINTS.faa'
MODELS = {
    'ESM2': ROOT / 'ps04_model_forward_v0.1.0_20261001/esm2',
    'ProtT5': ROOT / 'ps04_model_forward_v0.1.0_20261001/prott5',
    'ESM-C': ROOT / 'ps04_model_forward_v0.1.0_20261001/esmc',
    'CARP-640': ROOT / 'ps04_carp640_forward_v0.1.0_20261001/carp',
}


def main():
    # This sample verifies numerical repeatability only. It is fixed before
    # model scoring and does not use roles, credits, or model outcomes.
    sample_size = 10000
    with gzip.open(REL, 'rt', newline='') as handle:
        rows = []
        for row in csv.DictReader(handle, delimiter='\t'):
            rows.append(row)
            if len(rows) == sample_size:
                break
    endpoint_ids = []
    with ENDPOINTS.open() as handle:
        for line in handle:
            if line.startswith('>'):
                endpoint_ids.append(line[1:].strip().split()[0])
    ids = sorted(set(endpoint_ids))
    sample_digest = hashlib.sha256("\n".join("\t".join(row[k] for k in (
        'anchor_id','positive_id','negative_id','negative_hsf_id','negative_topology_id')) for row in rows).encode()).hexdigest()
    out = {'artifact_id': 'PS04_SAFE_COSINE_PREFLIGHT_v0.1.1_20261007', 'operator_id': OPERATOR_ID,
           'sample_definition': 'first_10000_rows_in_frozen_formal_relation_ledger_order',
           'sample_relation_rows': len(rows), 'sample_relation_sha256': sample_digest,
           'leaf_sides_per_model': len(rows) * 2, 'formal_endpoint_count': len(ids), 'models': {}}
    for name, directory in MODELS.items():
        with (directory / 'embedding_index.tsv').open(newline='') as handle:
            index = {row['domain_id']: row['embedding_key'] for row in csv.DictReader(handle, delimiter='\t')}
        if set(ids) - set(index):
            raise ValueError(f'{name}: formal preflight endpoint missing from index')
        with safe_open(str(directory / 'embeddings.safetensors'), framework='np') as store:
            vectors = {identifier: store.get_tensor(index[identifier]) for identifier in ids}
        zero_norm_endpoint_count = sum(float(np.linalg.norm(vectors[identifier].astype(np.float32))) == 0.0 for identifier in ids)
        first, second = [], []
        for row in rows:
            anchor = vectors[row['anchor_id']]
            first.extend((safe_cosine_total_float32(anchor, vectors[row['positive_id']]), safe_cosine_total_float32(anchor, vectors[row['negative_id']])))
            second.extend((safe_cosine_total_float32(anchor, vectors[row['positive_id']]), safe_cosine_total_float32(anchor, vectors[row['negative_id']])))
        finite = all(math.isfinite(float(value)) for value in first + second)
        difference = max(abs(float(a) - float(b)) for a, b in zip(first, second))
        if not finite or difference != 0.0:
            raise ValueError((name, finite, difference))
        out['models'][name] = {'endpoint_count': len(ids), 'embedding_dim': int(vectors[ids[0]].shape[0]),
                               'finite': finite, 'max_repeat_abs_difference': difference,
                               'zero_norm_endpoint_count': int(zero_norm_endpoint_count)}
    out['delta_num'] = 0
    out['status'] = 'PASS'
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()

