#!/usr/bin/env python3
"""Read-back audit of all exported D0/D2 coefficients and result metadata.

This does not reconstruct matrices; it validates the actual deliverable files,
checks integer types strictly, and emits a compact index for all 269 blocks.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent


def strict_entries(entries,size):
    table={}
    for row in entries:
        assert len(row)==4 and all(type(x) is int for x in row), row
        v,i,j,c=row
        assert 0<=v<=941 and 0<=i<=j<size and c!=0, row
        assert (v,i,j) not in table, row
        table[v,i,j]=c
    return table


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--model',type=Path,required=True)
    args=p.parse_args()
    start=time.perf_counter()
    out=HERE/'results'
    with gzip.open(args.model,'rt') as f:model=json.load(f)
    d0=json.loads((out/'d0'/'summary.json').read_text())
    d2=json.loads((out/'d2_summary.json').read_text())
    assert d0['blocks_requested']==list(range(65))
    assert d0['passed']==65 and d0['mismatched']==0
    assert d2['fully_equal']==204 and d2['mismatching']==0
    with gzip.open(out/'d2_reconstructed.json.gz','rt') as f:
        d2_coeff={b['block']:b['entries'] for b in json.load(f)}
    assert set(d2_coeff)==set(range(65,269))
    index=[]
    for bid in range(269):
        if bid<65:
            details=d0['blocks'][bid]
            assert details['block']==bid
            with gzip.open(out/'d0'/f'block_{bid:03d}.rebuilt.json.gz','rt') as f:r=json.load(f)
            entries=r['entries'];size=r['size'];g=r['positive_scale_removed'];kept=r['kept_columns']
            original=details['full_columns'];raw=details['uncompressed_nonzero_keys']
            elapsed=details['wall_seconds']
            reported_digest=details['coefficients_sha256']
            assert g==details['gcd_computed'] and kept==details['kept_columns']
        else:
            details=d2['blocks'][bid-65]
            assert details['block']==bid
            entries=d2_coeff[bid];size=details['stored_columns']
            g=details['positive_scale_removed'];kept=details['kept_columns']
            original=details['original_columns'];raw=details['raw_nonzero_keys']
            elapsed=details['reconstruction_seconds']
            reported_digest=details['reconstructed_entries_sha256']
        actual=strict_entries(entries,size)
        b=model['blocks'][bid]
        expected=strict_entries(b['entries'],b['size'])
        assert actual==expected, bid
        assert size==b['size']==len(kept), bid
        assert g==b['label']['positive_scale_removed'] and type(g) is int and g>0, bid
        assert kept==sorted(set(kept)) and all(0<=k<original for k in kept),bid
        digest=hashlib.sha256(json.dumps(sorted(entries),separators=(',',':')).encode()).hexdigest()
        assert digest==reported_digest,bid
        index.append({'block':bid,'family':'D0' if bid<65 else 'D2',
                      'full_columns':original,'kept_size':size,
                      'gcd':g,'kept_columns':kept,'raw_nonzero_keys':raw,
                      'stored_nonzero_keys':len(actual),'sha256':digest,
                      'reconstruction_seconds':elapsed,'status':'PASS'})
    (out/'block_index.json').write_text(json.dumps(index,indent=2)+'\n')
    with (out/'block_index.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(index[0]));w.writeheader()
        for row in index:
            row=dict(row);row['kept_columns']=','.join(map(str,row['kept_columns']));w.writerow(row)
    result={'purpose':'Strict integer schema and read-back verification of exported coefficients',
            'blocks_verified':269,'all_equal':True,'total_nonzero_keys':sum(b['stored_nonzero_keys'] for b in index),
            'full_columns':sum(b['full_columns'] for b in index),
            'kept_columns':sum(b['kept_size'] for b in index),
            'total_wall_seconds':time.perf_counter()-start}
    (out/'delivery_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
