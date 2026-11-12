"""Reconstruct the first five requested families, using definitions only."""
from collections import Counter
from itertools import combinations
from math import comb
from pathlib import Path
import argparse
import json
import time

from common import (N, W, V, BOUND, FULL, DATA, OrbitIndex, histogram,
                    restore_words, load_orbits, load_labels, clean_entries,
                    write_rebuilt)


def all_words():
    return [sum(1 << p for p in points) for points in combinations(range(N), W)]


def record(block, coeff, **metadata):
    return {'block_id': block['block_id'], 'size': 1,
            'entries': clean_entries({(v,0,0):c for v,c in coeff.items()}),
            'computed_metadata': metadata, 'metadata_issues': []}


def nonnegative():
    start = time.perf_counter()
    blocks = load_labels('nonnegative')
    assert [b['block_id'] for b in blocks] == list(range(269,1210))
    assert len(load_orbits()) == 941
    rows = [record(b, {b['block_id']-268:1}) for b in blocks]
    return write_rebuilt('nonnegative', rows, time.perf_counter()-start)


def cardinality_cap():
    start = time.perf_counter()
    blocks = load_labels('cardinality-cap')
    assert len(blocks) == 1
    bound = (18*124)//12
    assert bound == BOUND
    for b in blocks:
        assert b['label']['bound'] == bound
    rows = [record(b, {0:bound,1:-1}, bound=bound) for b in blocks]
    return write_rebuilt('cardinality-cap', rows, time.perf_counter()-start)


def extension_vector(words, index, universe):
    """Coefficients of L(X_words(186-M)) on a fixed support."""
    ans = Counter({index.v(words): BOUND})
    for u in universe:
        # Equality with an anchor word is allowed and collapses the support.
        if all(u == a or (u&a).bit_count() <= 3 for a in words):
            support = tuple(dict.fromkeys((*words,u)))
            v = index.by_hist[histogram(support)]
            ans[v] -= 1
    return {v:c for v,c in ans.items() if c}


def cardinality_scalar():
    start = time.perf_counter()
    index = OrbitIndex()
    universe = all_words()
    blocks = load_labels('cardinality-scalar')
    assert len(blocks) == 40
    expected = {tuple(r['h']) for r in index.records if r['k'] <= 3}
    assert set(tuple(b['label']['anchor']) for b in blocks) == expected
    rows, maps = [], {}
    for b in blocks:
        t = time.perf_counter()
        label = b['label']
        assert label['bound'] == BOUND
        words = restore_words(label['anchor'])
        v = index.v(words)
        coeff = extension_vector(words,index,universe)
        maps[str(v)] = sorted(coeff.items())
        row = record(b, coeff, anchor_variable=v)
        row['seconds'] = time.perf_counter()-t
        rows.append(row)
    (DATA/'extension_map.json').write_text(json.dumps(maps,separators=(',',':')))
    return write_rebuilt('cardinality-scalar',rows,time.perf_counter()-start,
                         universe_words=len(universe))


def matching_occupancy():
    start = time.perf_counter()
    index = OrbitIndex()
    blocks = load_labels('matching-occupancy')
    grouped = {}
    for b in blocks:
        l = b['label']
        key = (l['anchor'],l['lower'],l['upper'])
        grouped.setdefault(key,[]).append(b)
    computed = {}
    totals = [0,0,0,0]
    for (a,L,U), group in grouped.items():
        t = time.perf_counter()
        assert a == 63 and L & ~U == 0
        middle = U ^ L
        points = [i for i in range(N) if middle >> i & 1]
        m, k = len(points), W-L.bit_count()
        cap = m//3
        assert min(k,m-k) == 3 and cap == 3
        caps = [N-U.bit_count(),m,L.bit_count()]
        words = [L | sum(1<<p for p in sub) for sub in combinations(points,k)]
        # Multiplication by X_a makes every incompatible interval word vanish.
        words = [u for u in words if u==a or (u&a).bit_count()<=3]
        words.sort()
        neighbors = [[j for j in range(i+1,len(words))
                      if (words[i]&words[j]).bit_count() <= 3]
                     for i in range(len(words))]
        by_degree = [Counter() for _ in range(4)]
        by_degree[0][index.v((a,))] = 1
        totals[0] += 1
        for i,u in enumerate(words):
            by_degree[1][index.v((a,u))] += 1
            totals[1] += 1
            for j in neighbors[i]:
                v = words[j]
                by_degree[2][index.v((a,u,v))] += 1
                totals[2] += 1
                for l in neighbors[j]:
                    q = words[l]
                    if (u&q).bit_count() <= 3:
                        by_degree[3][index.v((a,u,v,q))] += 1
                        totals[3] += 1
        shared_seconds = time.perf_counter()-t
        for b in group:
            label = b['label']
            assert label['caps'] == caps and label['capacity'] == cap
            j = label['occupied']
            assert 0 <= j <= 2
            coeff = Counter()
            for degree in range(j,cap+1):
                multiplier = (-1)**(degree-j)*comb(degree,j)
                for v,c in by_degree[degree].items():
                    coeff[v] += multiplier*c
            row = record(b, coeff, capacity=cap,
                         interval_words=comb(m,k),
                         compatible_with_anchor=len(words))
            row['shared_group_seconds'] = shared_seconds
            computed[b['block_id']] = row
    return write_rebuilt('matching-occupancy',
                         [computed[b['block_id']] for b in blocks],
                         time.perf_counter()-start,
                         interval_groups=len(grouped),
                         enumerated_compatible_subsets_by_degree=totals)


def anticode_words(R, threshold):
    inside = [i for i in range(N) if R>>i&1]
    outside = [i for i in range(N) if not (R>>i&1)]
    for r in range(threshold,min(W,len(inside))+1):
        for a in combinations(inside,r):
            a_mask = sum(1<<i for i in a)
            for b in combinations(outside,W-r):
                yield a_mask | sum(1<<i for i in b)


def anticode_scalar():
    start = time.perf_counter()
    index = OrbitIndex()
    blocks = load_labels('anticode-vacancy-scalar')
    rows = []
    tested, retained, no_anchor = 0,0,0
    for b in blocks:
        t = time.perf_counter()
        label = b['label']
        D = restore_words(label['anchor'])
        assert len(D)==3
        v0 = index.v(D)
        R,r = label['subset'],label['r']
        threshold = 4+r
        assert r in (0,1,2) and R.bit_count()==4+2*r
        assert label['threshold'] == threshold and R & ~FULL == 0
        assert all((a&R).bit_count() < threshold for a in D)
        # Also verify the selected R is the initial segment within every cell.
        offset = 0
        for c in label['anchor']:
            segment = ((1<<c)-1)<<offset
            count = (R&segment).bit_count()
            assert R&segment == ((1<<count)-1)<<offset
            offset += c
        coeff = Counter({v0:1})
        for u in anticode_words(R,threshold):
            tested += 1
            if all(u==a or (u&a).bit_count()<=3 for a in D):
                retained += 1
                v = index.by_hist[histogram((*D,u))]
                coeff[v] -= 1
        no_anchor += 1
        row = record(b,coeff,anchor_variable=v0)
        row['seconds'] = time.perf_counter()-t
        rows.append(row)
    return write_rebuilt('anticode-vacancy-scalar',rows,time.perf_counter()-start,
                         anticode_words_tested=tested,
                         compatible_extensions=retained,
                         no_anchor_in_anticode_checked=no_anchor)


FUNCTIONS = {'nonnegative':nonnegative, 'cardinality-cap':cardinality_cap,
             'cardinality-scalar':cardinality_scalar,
             'matching-occupancy':matching_occupancy,
             'anticode-vacancy-scalar':anticode_scalar}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('families',nargs='*',choices=list(FUNCTIONS))
    args = parser.parse_args()
    for family in args.families or FUNCTIONS:
        path = FUNCTIONS[family]()
        result = json.loads(path.read_text())
        print(json.dumps({'family':family,'blocks':len(result['blocks']),
                          'seconds':result['seconds'],'output':str(path)}),flush=True)
