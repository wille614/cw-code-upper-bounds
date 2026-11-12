#!/usr/bin/env python3
"""Independent D0 driver: reconstruction -> exact Gram selection -> comparison.

The C++ engine has no access to expected entries, g, or stored block sizes.
An explicit coefficients file must exist before this driver examines the
corresponding model entries. No supplied numerical example is used.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

HERE = Path(__file__).resolve().parent


def exact_kept_columns(gram):
    """Exact rational Schur elimination in original column order.

    A zero PSD Schur diagonal means a dependent column; its entire residual
    row must be zero. A positive diagonal is the next independent column.
    """
    a = [[Fraction(x) for x in row] for row in gram]
    kept, pivots = [], []
    n = len(a)
    for k in range(n):
        p = a[k][k]
        if p < 0:
            raise ArithmeticError(f"negative Gram pivot {k}: {p}")
        if not p:
            if any(a[k][j] for j in range(k+1,n)):
                raise ArithmeticError(f"zero Gram pivot with nonzero row {k}")
            continue
        kept.append(k)
        pivots.append(str(p))
        for i in range(k+1,n):
            if not a[i][k]:
                continue
            for j in range(i,n):
                if a[j][k]:
                    a[i][j] -= a[i][k]*a[j][k]/p
                    a[j][i] = a[i][j]
    return kept,pivots


def prepare_labels(model_path, out):
    # Only labels are retained; no block entries are used in reconstruction.
    with gzip.open(model_path,'rt') as f:
        model = json.load(f)
    labels = [dict(block['label']) for block in model['blocks'][:65]]
    with out.open('w') as f:
        for bid,label in enumerate(labels):
            assert label['family']=='D0'
            shape = label['partition']+[0]*(4-len(label['partition']))
            f.write(' '.join(map(str,[bid,len(label['tableaux'])]+shape))+'\n')
            for tableau in label['tableaux']:
                f.write(' '.join(str(x) for row in tableau for x in row)+'\n')
    return labels


def sha_entries(entries):
    return hashlib.sha256(json.dumps(sorted(entries),separators=(',',':')).encode()).hexdigest()


def compress_block(bid,label,raw_prefix,variables):
    started = time.perf_counter()
    coeff={}
    with Path(str(raw_prefix)+'.raw.tsv').open() as f:
        for line in f:
            v,i,j,c=map(int,line.split())
            assert c and (v,i,j) not in coeff
            coeff[v,i,j]=c
    metrics=json.loads(Path(str(raw_prefix)+'.metrics.json').read_text())
    n=metrics['full_columns']
    sizes=[1 if sum(row.count(1) for row in tableau)==0 else 2
           for tableau in label['tableaux']]
    # Allowed contents are (12-t,t,t,6-t); t=0 corresponds to one word.
    if label['partition']==[18]:sizes.insert(0,0)
    assert len(sizes)==n
    g=math.gcd(*coeff.values())
    assert g>0
    gram=[[0]*n for _ in range(n)]
    arity=[0]+[(len(h)-1).bit_length() for h in variables]
    for (v,i,j),c in coeff.items():
        if sizes[i]==sizes[j] and sizes[i]>0 and arity[v]==sizes[i]:
            gram[i][j]+=c
            if j!=i:gram[j][i]+=c
    if sizes[0]==0:gram[0][0]=1
    kept,pivots=exact_kept_columns(gram)
    old_to_new={k:i for i,k in enumerate(kept)}
    entries=[]
    for (v,i,j),c in coeff.items():
        if i in old_to_new and j in old_to_new:
            assert c%g==0
            entries.append([v,old_to_new[i],old_to_new[j],c//g])
    entries.sort()
    with gzip.open(str(raw_prefix)+'.rebuilt.json.gz','wt',encoding='utf-8') as f:
        json.dump({'block':bid,'size':len(kept),'positive_scale_removed':g,
                   'kept_columns':kept,'entries':entries},f,separators=(',',':'))
    metrics.update({'block':bid,'family':'D0','full_columns':n,'kept_columns':kept,
                    'kept_size':len(kept),'gcd_computed':g,
                    'uncompressed_nonzero_keys':len(coeff),
                    'rebuilt_nonzero_keys':len(entries),
                    'coefficients_sha256':sha_entries(entries),
                    'gram':gram,'schur_pivots':pivots,
                    'compression_seconds':time.perf_counter()-started})
    return metrics,entries


def compare_block(metrics, entries, block):
    expected={}
    invalid=[]
    for e in block['entries']:
        if len(e)!=4 or not all(type(a) is int for a in e):
            invalid.append({'entry':e,'reason':'invalid integer quadruple'});continue
        v,i,j,c=e
        if not (0<=v<=941 and 0<=i<=j<block['size'] and c):
            invalid.append({'entry':e,'reason':'invalid index or explicit zero'})
        if (v,i,j) in expected:invalid.append({'entry':e,'reason':'duplicate key'})
        expected[v,i,j]=c
    actual={tuple(e[:3]):e[3] for e in entries}
    differences=[]
    for key in sorted(expected.keys()|actual.keys()):
        if expected.get(key,0)!=actual.get(key,0):
            differences.append({'key':list(key),'model':expected.get(key,0),'computed':actual.get(key,0)})
    metrics.update({'model_nonzero_keys':len(expected),'model_invalid_entries':invalid,
                    'coefficient_differences':differences,
                    'gcd_model':block['label']['positive_scale_removed'],
                    'size_model':block['size'],
                    'all_keys_values_equal':not differences and not invalid,
                    'gcd_equal':metrics['gcd_computed']==block['label']['positive_scale_removed'],
                    'size_equal':metrics['kept_size']==block['size']})
    metrics['pass']=all(metrics[x] for x in ('all_keys_values_equal','gcd_equal','size_equal'))
    return metrics


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--model',type=Path,required=True)
    ap.add_argument('--blocks',default='all',help='comma separated indices, or all')
    ap.add_argument('--workers',type=int,default=2)
    ap.add_argument('--reuse-raw',action='store_true')
    ap.add_argument('--no-compile',action='store_true')
    args=ap.parse_args()
    total_start=time.perf_counter()
    outdir=HERE/'results'/'d0';outdir.mkdir(parents=True,exist_ok=True)
    inputfile=outdir/'labels.txt'
    labels=prepare_labels(args.model,inputfile)
    variables=json.loads((HERE/'results'/'canonical_variables.json').read_text())
    if isinstance(variables,dict):variables=variables['variables']
    engine=HERE/'d0_engine'
    if not args.no_compile:
        subprocess.run(['g++','-O3','-std=c++17','-Wall','-Wextra',str(HERE/'d0_engine.cpp'),'-o',str(engine)],check=True)
    blocks=list(range(65)) if args.blocks=='all' else [int(s) for s in args.blocks.split(',')]
    # Small full systems first: stored sizes are not consulted for scheduling.
    blocks.sort(key=lambda b:(len(labels[b]['tableaux'])+(labels[b]['partition']==[18]),b))
    def compute(bid):
        prefix=outdir/f'block_{bid:03d}'
        if not (args.reuse_raw and Path(str(prefix)+'.raw.tsv').exists() and Path(str(prefix)+'.metrics.json').exists()):
            with Path(str(prefix)+'.engine.log').open('w') as log:
                subprocess.run([str(engine),str(inputfile),str(HERE/'results'/'orbit_lookup.tsv'),str(bid),str(prefix)],stdout=log,stderr=subprocess.STDOUT,check=True)
        metrics,entries=compress_block(bid,labels[bid],prefix,variables)
        print(f"RECONSTRUCTED D0 {bid}: full={metrics['full_columns']} kept={metrics['kept_size']} g={metrics['gcd_computed']} keys={len(entries)} seconds={metrics['wall_seconds']}",flush=True)
        return bid,metrics,entries
    computed={}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs={pool.submit(compute,bid):bid for bid in blocks}
        for future in as_completed(jobs):
            bid,metrics,entries=future.result()
            computed[bid]=(metrics,entries)
    # Model entries are first accessed only after every requested reconstruction.
    with gzip.open(args.model,'rt') as f:model=json.load(f)
    results=[]
    for bid in sorted(computed):
        metrics,entries=computed[bid]
        compare_block(metrics,entries,model['blocks'][bid])
        (outdir/f'block_{bid:03d}.result.json').write_text(json.dumps(metrics,indent=2)+'\n')
        results.append(metrics)
        print(f"{'PASS' if metrics['pass'] else 'FAIL'} D0 {bid}: differences={len(metrics['coefficient_differences'])} g_equal={metrics['gcd_equal']} size_equal={metrics['size_equal']}",flush=True)
    summary={'blocks_requested':sorted(blocks),'blocks_reconstructed':len(results),
             'passed':sum(r['pass'] for r in results),
             'mismatched':sum(not r['pass'] for r in results),
             'not_reconstructed':65-len(results),
             'full_columns':sum(r['full_columns'] for r in results),
             'kept_columns':sum(r['kept_size'] for r in results),
             'uncompressed_nonzero_keys':sum(r['uncompressed_nonzero_keys'] for r in results),
             'compared_nonzero_keys':sum(r['rebuilt_nonzero_keys'] for r in results),
             'sum_engine_wall_seconds':sum(r['wall_seconds'] for r in results),
             'sum_engine_cpu_seconds':sum(r['cpu_seconds'] for r in results),
             'sum_compression_seconds':sum(r['compression_seconds'] for r in results),
             'total_wall_seconds':time.perf_counter()-total_start,
             'worker_count':args.workers,
             'blocks':[{k:v for k,v in r.items() if k not in ('gram','schur_pivots')} for r in results]}
    (outdir/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='blocks'},indent=2),flush=True)
    if summary['mismatched']:raise SystemExit(1)

if __name__=='__main__':main()
