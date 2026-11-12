#!/usr/bin/env python3
"""Resolve the five raw-tableau shapes omitted from the stored D0 blocks.

Generate their tableaux independently, use the audited exact D0 engine, and
inspect every uncompressed coefficient. No model entries or artificial gcd
are involved. A zero Gram matrix proves every unordered-pair projection
column is zero.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

from check_labels import generate_shape

HERE=Path(__file__).resolve().parent
OMITTED=[(15,1,1,1),(9,3,3,3),(8,5,3,2),(7,7,3,1),(7,5,5,1)]


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--no-compile',action='store_true')
    args=ap.parse_args()
    start=time.perf_counter()
    dest=HERE/'results'/'omitted_shapes'
    dest.mkdir(parents=True,exist_ok=True)
    labels=dest/'labels.txt'
    candidates=[]
    with labels.open('w') as f:
        for i,shape in enumerate(OMITTED):
            bid=1000+i
            tableaux,counts=generate_shape(shape)
            candidates.append((bid,shape,tableaux,counts))
            f.write(' '.join(map(str,[bid,len(tableaux),*shape]))+'\n')
            for tab in tableaux:
                f.write(' '.join(map(str,sum(tab,[])))+'\n')
    engine=dest/'d0_engine'
    source=HERE/'d0_engine.cpp'
    if not args.no_compile:
        subprocess.run(['g++','-O3','-std=c++17','-Wall','-Wextra',str(source),'-o',str(engine)],check=True)
    variables=json.loads((HERE/'results'/'canonical_variables.json').read_text())
    arities=[0]+[(len(h)-1).bit_length() for h in variables]
    results=[]
    for bid,shape,tableaux,counts in candidates:
        run_start=time.perf_counter()
        prefix=dest/f'shape_{bid}'
        with Path(str(prefix)+'.engine.log').open('w') as log:
            subprocess.run([str(engine),str(labels),str(HERE/'results'/'orbit_lookup.tsv'),str(bid),str(prefix)],stdout=log,stderr=subprocess.STDOUT,check=True)
        metrics=json.loads(Path(str(prefix)+'.metrics.json').read_text())
        n=len(tableaux)
        row_sizes=[1 if sum(row.count(1) for row in tab)==0 else 2 for tab in tableaux]
        gram=[[0]*n for _ in range(n)]
        raw=[]
        with Path(str(prefix)+'.raw.tsv').open() as f:
            for line in f:
                v,i,j,c=map(int,line.split())
                assert c!=0
                raw.append([v,i,j,c])
                if row_sizes[i]==row_sizes[j] and arities[v]==row_sizes[i]:
                    gram[i][j]+=c
                    if i!=j:
                        gram[j][i]+=c
        result={
            'independent_test_id':bid,
            'partition':shape,
            'tableaux':tableaux,
            'tableaux_by_t':counts,
            'full_columns':n,
            'all_upper_triangle_pairs_evaluated':metrics['tableau_pairs']==n*(n+1)//2,
            'raw_nonzero_key_count':len(raw),
            'raw_all_zero':not raw,
            'gram':gram,
            'gram_all_zero':all(x==0 for row in gram for x in row),
            'elapsed_seconds':time.perf_counter()-run_start,
            'engine_metrics':metrics,
            'gcd_computed':False,
        }
        (dest/f'shape_{bid}.result.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
        results.append(result)
        print(json.dumps({k:v for k,v in result.items() if k not in ('tableaux','gram','engine_metrics')},ensure_ascii=False),flush=True)
    summary={
        'purpose':'Resolve exactly the five previously identified raw-tableau shapes omitted from D0.',
        'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'shape_count':len(results),
        'raw_all_zero_count':sum(r['raw_all_zero'] for r in results),
        'gram_all_zero_count':sum(r['gram_all_zero'] for r in results),
        'all_pairs_evaluated':all(r['all_upper_triangle_pairs_evaluated'] for r in results),
        'sum_engine_wall_seconds':sum(r['engine_metrics']['wall_seconds'] for r in results),
        'sum_engine_cpu_seconds':sum(r['engine_metrics']['cpu_seconds'] for r in results),
        'elapsed_seconds':time.perf_counter()-start,
        'conclusion':'For each omitted shape P^T P=0; hence P=0, and no nonzero projected D0 shape is missing among partitions with at most four rows and allowed contents.' if all(r['gram_all_zero'] for r in results) else 'At least one omitted shape has a nonzero projection; see individual results.',
        'model_entries_accessed':False,
        'results':results,
    }
    (dest/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='results'},ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':
    main()
