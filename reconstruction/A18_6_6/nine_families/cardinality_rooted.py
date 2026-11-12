#!/usr/bin/env python3
"""Independent exact reconstruction of cardinality-rooted (SPEC01 sections 2,13).

Only data/labels.json, independently built orbit data, and independently built
T_186 extension maps are consumed. No model entries or appendix examples are read.
The binary Young-pair formula below is derived from the literal row-symmetrized,
column-antisymmetrized definition and checked against its explicit expansion.
"""
from __future__ import annotations
from collections import defaultdict
from functools import lru_cache
from itertools import combinations, product
from math import comb, factorial, gcd
from pathlib import Path
import argparse
import json
import time

from common import (DATA, OrbitIndex, clean_entries, load_labels, restore_words,
                    gcd_values, write_rebuilt)


def _multinomial(counts):
    n = sum(counts)
    value = factorial(n)
    for count in counts:
        value //= factorial(count)
    return value


@lru_cache(None)
def binary_pair_polynomial(c, h, left_weight, right_weight):
    """Return ((histogram_00_01_10_11, integer coefficient), ...).

    A semistandard binary tableau of shape (c-h,h) has an all-1 bottom
    row. Top-row placements putting a 1 above a 1 cancel under the column
    transposition. Consequently its nonzero vector values are +/-1:
    the h paired columns must have distinct bits, the sign is (-1) to
    the number of top bits equal to 1, and the c-2h remaining top boxes
    have weight left_weight-h. A pair of such columns contributes
    2*(Y00*Y11-Y01*Y10); the unpaired boxes contribute the multinomial
    count for their joint histogram. No row factorial or word-order
    divisor is introduced.
    """
    assert 0 <= h <= c//2
    assert h <= left_weight <= c-h and h <= right_weight <= c-h
    length = c-2*h
    e, f = left_weight-h, right_weight-h
    answer = defaultdict(int)
    for both in range(max(0, e+f-length), min(e, f)+1):
        single = (length-e-f+both, f-both, e-both, both)
        multiplicity = _multinomial(single)
        for swapped in range(h+1):
            paired = (h-swapped, swapped, swapped, h-swapped)
            hist = tuple(a+b for a,b in zip(single,paired))
            answer[hist] += multiplicity*(1<<h)*((-1)**swapped)*comb(h,swapped)
    return tuple(sorted((hist, value) for hist,value in answer.items() if value))


def literal_binary_young(c, h, weight):
    """Literal distinct top-row permutations and column signed sums.

    This independent direct implementation is used only for algebra self-checks.
    Row-major positions: top 0..c-h-1, bottom c-h..c-1.
    """
    top = c-h
    result = defaultdict(int)
    for ones in combinations(range(top), weight-h):
        initial = sum(1<<i for i in ones)
        initial |= ((1<<h)-1)<<top
        for switches in range(1<<h):
            mask = initial
            for j in range(h):
                if switches>>j & 1:
                    b1, b2 = (mask>>j)&1, (mask>>(top+j))&1
                    if b1 != b2:
                        mask ^= (1<<j)|(1<<(top+j))
            result[mask] += (-1)**switches.bit_count()
    return {mask:value for mask,value in result.items() if value}


def explicit_pair_polynomial(c, left, right):
    result = defaultdict(int)
    full = (1<<c)-1
    for x,cx in left.items():
        for y,cy in right.items():
            h11 = (x&y).bit_count()
            h10 = (x&(full^y)).bit_count()
            h01 = ((full^x)&y).bit_count()
            hist = (c-h11-h10-h01,h01,h10,h11)
            result[hist] += cx*cy
    return tuple(sorted((hist,value) for hist,value in result.items() if value))


def self_check(max_c=6):
    start = time.perf_counter()
    pair_count = 0
    for c in range(max_c+1):
        for h in range(c//2+1):
            vectors = {w:literal_binary_young(c,h,w) for w in range(h,c-h+1)}
            for w,vector in vectors.items():
                assert all(value in (-1,1) for value in vector.values())
                assert all(mask.bit_count()==w for mask in vector)
            for left_weight in vectors:
                for right_weight in vectors:
                    expected=explicit_pair_polynomial(c,vectors[left_weight],vectors[right_weight])
                    actual=binary_pair_polynomial(c,h,left_weight,right_weight)
                    assert actual==expected, (c,h,left_weight,right_weight,actual,expected)
                    pair_count+=1
    return {'max_group_size':max_c,'ordered_tableau_pairs':pair_count,
            'seconds':time.perf_counter()-start,'result':'PASS'}


def _metadata_issue(issues, field, expected, actual):
    if expected != actual:
        issues.append({'field':field,'model':expected,'computed':actual})


def _tableau_weight(tableau, shape):
    assert len(tableau)==len(shape)
    assert all(len(row)==length for row,length in zip(tableau,shape))
    assert all(x in (0,1) for row in tableau for x in row)
    assert all(row==sorted(row) for row in tableau)
    for r in range(1,len(tableau)):
        assert all(tableau[r-1][j]<value for j,value in enumerate(tableau[r]))
    return sum(sum(row) for row in tableau)


def reconstruct_block(block, orbit, extension):
    label = block['label']
    source = label['source']
    issues = []
    for field, expected in [('family','D1'),('anchors',[63]),('k',3),('row_arity',1),
                            ('empty_column',False),
                            ('cells',[list(range(6,18)),list(range(6))])]:
        _metadata_issue(issues,'source.'+field,source.get(field),expected)
    _metadata_issue(issues,'bound',label['bound'],186)
    shapes=source['shapes']
    assert len(shapes)==2
    params=[]
    for cell,shape,n in zip(source['cells'],shapes,(12,6)):
        assert 1<=len(shape)<=2 and shape==sorted(shape,reverse=True)
        assert sum(shape)==n and len(cell)==n
        params.append((n,shape[1] if len(shape)==2 else 0))
    weights=[]
    arities=[]
    for column in source['columns']:
        assert len(column)==2
        current=tuple(_tableau_weight(tableau,shape) for tableau,shape in zip(column,shapes))
        assert sum(current)==6
        assert current[1]<=3 or current==(0,6)
        weights.append(current)
        arities.append(1 if current==(0,6) else 2)
    _metadata_issue(issues,'source.row_cardinalities',source['row_cardinalities'],arities)
    raw_source=defaultdict(int)
    dimension=len(weights)
    for i,left in enumerate(weights):
        for j in range(i,dimension):
            right=weights[j]
            external=binary_pair_polynomial(*params[0],left[0],right[0])
            internal=binary_pair_polynomial(*params[1],left[1],right[1])
            for h0,c0 in external:
                for h1,c1 in internal:
                    # Concatenation records the three bits (anchor,x,y).
                    words=restore_words(h0+h1)
                    v=orbit.v(words)
                    if v is not None:
                        raw_source[(v,i,j)] += c0*c1
    raw_source={key:value for key,value in raw_source.items() if value}
    gs=gcd_values(raw_source)
    assert gs>0
    _metadata_issue(issues,'source.positive_scale_removed',source['positive_scale_removed'],gs)
    normalized_source={key:value//gs for key,value in raw_source.items()}
    assert all(value%gs==0 for value in raw_source.values())
    raw_localized=defaultdict(int)
    for (v,i,j),coefficient in normalized_source.items():
        for target,value in extension[str(v)]:
            raw_localized[(target,i,j)] += coefficient*value
    raw_localized={key:value for key,value in raw_localized.items() if value}
    go=gcd_values(raw_localized)
    assert go>0
    _metadata_issue(issues,'positive_scale_removed',label['positive_scale_removed'],go)
    kept=sorted({i for v,i,j in raw_localized} | {j for v,i,j in raw_localized})
    _metadata_issue(issues,'kept_rows',label['kept_rows'],kept)
    _metadata_issue(issues,'size',block['size'],len(kept))
    remap={old:new for new,old in enumerate(kept)}
    normalized={(v,remap[i],remap[j]):coefficient//go
                for (v,i,j),coefficient in raw_localized.items()}
    assert all(coefficient%go==0 for coefficient in raw_localized.values())
    return {'block_id':block['block_id'],'size':len(kept),
            'entries':clean_entries(normalized),
            'computed_metadata':{'source_positive_scale_removed':gs,
                'positive_scale_removed':go,'kept_rows':kept,
                'source_dimension':dimension,'source_entry_count':len(normalized_source),
                'source_row_cardinalities':arities},
            'metadata_issues':issues}, {
                'block_id':block['block_id'],'size':dimension,
                'raw_entries':clean_entries(raw_source),
                'positive_scale_removed':gs,
                'normalized_entries':clean_entries(normalized_source)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-check-only',action='store_true')
    args=parser.parse_args()
    if args.self_check_only:
        check=self_check()
        (DATA/'cardinality_rooted_selfcheck.json').write_text(json.dumps(check,indent=2))
        print(json.dumps(check),flush=True)
        return
    start=time.perf_counter()
    orbit=OrbitIndex()
    extension=json.loads((DATA/'extension_map.json').read_text())
    if 'maps' in extension:
        extension=extension['maps']
    records=[]
    sources=[]
    for block in load_labels('cardinality-rooted'):
        t0=time.perf_counter()
        result,source=reconstruct_block(block,orbit,extension)
        result['seconds']=time.perf_counter()-t0
        records.append(result)
        sources.append(source)
        print(json.dumps({'block_id':result['block_id'],'entries':len(result['entries']),
                          'seconds':result['seconds'],'metadata_issues':result['metadata_issues']}),flush=True)
    elapsed=time.perf_counter()-start
    (DATA/'cardinality_rooted_sources.json').write_text(json.dumps(sources,separators=(',',':')))
    path=write_rebuilt('cardinality-rooted',records,elapsed,
                       method='Binary Young columns, exact joint histograms, T_186 substitution')
    print(json.dumps({'output':str(path),'blocks':len(records),'seconds':elapsed}),flush=True)


if __name__=='__main__':
    main()
