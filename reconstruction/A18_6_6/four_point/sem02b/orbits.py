#!/usr/bin/env python3
"""Independent exhaustive enumeration of compatible word-set orbits.

Only model.variables is read for comparison, after generating all candidates.
No block entries or supplied example coefficient is used.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

N, W, MAX_INTERSECTION = 18, 6, 3


def bit_patterns(k):
    return [tuple((p >> (k-1-j)) & 1 for j in range(k)) for p in range(1 << k)]


def permutation_maps(k):
    bits = bit_patterns(k)
    return [tuple(sum(bits[p][j] << (k-1-q) for q,j in enumerate(perm))
                  for p in range(1 << k))
            for perm in itertools.permutations(range(k))]


def permute_hist(h, mapping):
    out = [0] * len(h)
    for p, q in enumerate(mapping):
        out[q] = h[p]
    return tuple(out)


def canonical(h, maps):
    return min(permute_hist(h, m) for m in maps)


def extensions(h, k):
    """Append one word to an ordered k-word histogram, without model data.

    a[p] counts new-word ones among coordinates in old pattern p. All
    integral a with sum 6 and each old/new intersection <=3 are exhausted.
    """
    patterns = bit_patterns(k)
    n = len(h)
    suffix = [0] * (n+1)
    for p in range(n-1,-1,-1):
        suffix[p] = suffix[p+1] + h[p]
    a = [0]*n
    intersections = [0]*k
    def visit(p, remaining):
        if p == n:
            if remaining == 0:
                yield tuple(x for hp,ap in zip(h,a) for x in (hp-ap,ap))
            return
        lo = max(0,remaining-suffix[p+1])
        hi = min(h[p],remaining)
        for j, bit in enumerate(patterns[p]):
            if bit:
                hi = min(hi,MAX_INTERSECTION-intersections[j])
        for ap in range(lo,hi+1):
            a[p] = ap
            for j,bit in enumerate(patterns[p]):
                intersections[j] += bit*ap
            yield from visit(p+1,remaining-ap)
            for j,bit in enumerate(patterns[p]):
                intersections[j] -= bit*ap
    yield from visit(0,W)


def histogram_of_words(words):
    k = len(words)
    h = [0]*(1 << k)
    for i in range(N):
        p = 0
        for w in words:
            p = 2*p + ((w >> i)&1)
        h[p] += 1
    return tuple(h)


def words_of_histogram(h):
    k = (len(h)-1).bit_length()
    words = [0]*k
    start = 0
    for p,count in enumerate(h):
        for j in range(k):
            if p >> (k-1-j)&1:
                words[j] |= ((1<<count)-1)<<start
        start += count
    assert start == N
    return words


def ordered_four_lookup(variables):
    """All compatible ordered four-word histograms, allowing repetitions.

    Surjections from four positions onto each 1..4-word orbit representative
    cover all ordered quadruples, including every repeated-word case.
    """
    out = {}
    for v,h in enumerate(variables,1):
        k = (len(h)-1).bit_length()
        words = words_of_histogram(h)
        for assignment in itertools.product(range(k), repeat=4):
            if len(set(assignment)) != k:
                continue
            h4 = histogram_of_words([words[j] for j in assignment])
            previous = out.setdefault(h4,v)
            assert previous == v
    return out


def enumerate_orbits():
    all_ordered = [[(N-W,W)]]
    for k in range(1,4):
        next_layer = [out for h in all_ordered[-1] for out in extensions(h,k)]
        assert len(set(next_layer)) == len(next_layer)
        all_ordered.append(sorted(next_layer))
    layers = [sorted({canonical(h,permutation_maps(k)) for h in ordered})
              for k,ordered in enumerate(all_ordered,1)]
    return all_ordered,layers


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--model',default='upload/A18_6_6_model.json(1).gz')
    parser.add_argument('--out',default='sem02b/results')
    args=parser.parse_args()
    start=time.perf_counter()
    ordered,layers=enumerate_orbits()
    enumeration_seconds=time.perf_counter()-start
    variables=[h for layer in layers for h in layer]
    out=Path(args.out)
    out.mkdir(parents=True,exist_ok=True)
    (out/'canonical_variables.json').write_text(json.dumps(variables,separators=(',',':'))+'\n')
    lookup=ordered_four_lookup(variables)
    with (out/'orbit_lookup.tsv').open('w') as f:
        for h,v in sorted(lookup.items()):
            f.write('\t'.join(map(str,(*h,v)))+'\n')
    with (out/'orbit_lookup_packed.tsv').open('w') as f:
        for h,v in sorted(lookup.items()):
            assert max(h)<16
            packed=sum(count<<(4*p) for p,count in enumerate(h))
            f.write(f'{packed}\t{v}\n')
    # The independently constructed table is complete before model is opened.
    data=Path(args.model).read_bytes()
    raw=gzip.decompress(data)
    model=json.loads(raw)
    mv=model['variables']
    mismatches=[]
    for i in range(max(len(variables),len(mv))):
        got=list(variables[i]) if i<len(variables) else None
        expected=mv[i] if i<len(mv) else None
        if got!=expected:
            mismatches.append({'v':i+1,'computed':got,'model':expected})
    orbit_stats=[]
    for k,layer in enumerate(layers,1):
        maps=permutation_maps(k)
        totals=0
        stabilizers={}
        for h in layer:
            stab=sum(permute_hist(h,m)==h for m in maps)
            numer=math.factorial(N)
            denom=math.prod(math.factorial(x) for x in h)*stab
            assert numer%denom==0
            orbit_size=numer//denom
            totals += orbit_size
            stabilizers[str(stab)]=stabilizers.get(str(stab),0)+1
        orbit_stats.append({'k':k,'ordered_histograms':len(ordered[k-1]),
                            'orbits':len(layer),'sum_orbit_sizes':totals,
                            'word_stabilizer_counts':stabilizers})
    # First two orbit totals have separate elementary counts.
    expected_pairs=math.comb(N,W)*sum(math.comb(W,i)*math.comb(N-W,W-i)
                                     for i in range(MAX_INTERSECTION+1))//2
    assert orbit_stats[0]['sum_orbit_sizes']==math.comb(N,W)
    assert orbit_stats[1]['sum_orbit_sizes']==expected_pairs
    result={
        'method':'Exhaustive integer histogram refinement; canonization over all word permutations.',
        'model_gzip_sha256':hashlib.sha256(data).hexdigest(),
        'model_json_sha256':hashlib.sha256(raw).hexdigest(),
        'ordered_histogram_counts':[len(x) for x in ordered],
        'orbit_counts':[len(x) for x in layers],
        'variables_total':len(variables),
        'variables_exactly_equal':not mismatches,
        'variable_mismatches':mismatches,
        'orbit_sizes':orbit_stats,
        'ordered_four_lookup_count':len(lookup),
        'lookup_by_distinct_word_count':{str(k):sum(len(variables[v-1])==(1<<k) for v in lookup.values()) for k in range(1,5)},
        'enumeration_seconds':enumeration_seconds,
        'total_seconds':time.perf_counter()-start,
        'arithmetic':'Python arbitrary-precision integers; no floating arithmetic except elapsed wall time.',
        'model_blocks_entries_accessed':False,
        'ambiguities':[],
    }
    (out/'orbits.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
