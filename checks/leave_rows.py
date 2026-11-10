from pathlib import Path

from itertools import combinations, permutations

from functools import lru_cache

from math import comb, factorial, gcd

import json, time, hashlib

def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0

def representative(hist):
    r = (len(hist) - 1).bit_length()
    words = [set() for _ in range(r)]
    p = 0
    for mask, count in enumerate(hist):
        for _ in range(count):
            for j in range(r):
                if mask & 1 << j:
                    words[j].add(p)
            p += 1
    return tuple((frozenset(x) for x in words))

def orbit_cardinality(hist):
    r = (len(hist) - 1).bit_length()
    stabilizer = 0
    for perm in permutations(range(r)):
        mapped = [0] * len(hist)
        for mask, v in enumerate(hist):
            new = sum(((mask >> j & 1) << perm[j] for j in range(r)))
            mapped[new] = v
        stabilizer += mapped == hist
    denominator = stabilizer
    for v in hist:
        denominator *= factorial(v)
    q, r = divmod(factorial(sum(hist)), denominator)
    assert not r
    return q

@lru_cache(None)
def graph_values(m, blocks):
    points = set(range(m))
    edges = {frozenset(x) for x in combinations(points, 2)}
    for block in blocks:
        for e in combinations(block, 2):
            edge = frozenset(e)
            assert edge in edges, 'not a packing'
            edges.remove(edge)
    adjacency = [set() for _ in range(m)]
    for u, v in edges:
        adjacency[u].add(v)
        adjacency[v].add(u)
    triangles = sum((b in adjacency[a] and c in adjacency[a] and (c in adjacency[b]) for a, b, c in combinations(range(m), 3)))
    squares = 0
    for a, b, c, d in combinations(range(m), 4):
        for order in ((a, b, c, d), (a, b, d, c), (a, c, b, d)):
            squares += all((order[(i + 1) % 4] in adjacency[order[i]] for i in range(4)))
    k = len(blocks[0]) if blocks else 4
    assert k == 4
    t = comb(m, 2) // 6 - len(blocks)
    F = triangles + 4 * t * (t - 2)
    return (F, len(blocks) * F, t * F, triangles + squares + 5 * t * t - 12 * t, 24 * (triangles + squares) + t * (t - 3) * (t - 5) * (5 * t - 26))

def residual_coefficient(words, m, ell, q):
    common = set.intersection(*(set(x) for x in words)) if words else set()
    if words:
        assert len(common) >= ell
        core = set(sorted(common)[:ell])
        union = set().union(*words) - core
        relabel = {x: i for i, x in enumerate(sorted(union))}
        residual = tuple((tuple(sorted((relabel[x] for x in u - core))) for u in words))
    else:
        residual = ()
    total = 0
    for size in range(len(words) + 1):
        for subset in combinations(residual, size):
            total += (-1) ** (len(words) - size) * graph_values(m, subset)[q]
    return total

def interval_count_dp(words, anchor, n, ell, b, j, h):
    dp = {(0, 0, 0, 0): 1}
    for x in range(n):
        can_core = all((x in u for u in words))
        can_excluded = all((x not in u for u in words))
        in_anchor = int(x in anchor)
        new = {}
        for (rc, zc, ac, az), value in dp.items():
            for key in ((rc, zc, ac, az), *([(rc + 1, zc, ac + in_anchor, az)] if can_core else []), *([(rc, zc + 1, ac, az + in_anchor)] if can_excluded else [])):
                if key[0] <= ell and key[1] <= b and (key[2] <= j) and (key[3] <= h):
                    new[key] = new.get(key, 0) + value
        dp = new
    return dp.get((ell, b, j, h), 0)

def check(model, detail):
    n=18; V=comb(n,6)
    histograms=model['variables']; reps=[representative(h) for h in histograms]
    sizes=[orbit_cardinality(h) for h in histograms]
    allpoints=set(range(n));results=[]
    for idx,block in enumerate(model['blocks']):
        label=block.get('label',{});family=label.get('family','')
        if not family.startswith('leave-motif'):continue
        ell,b,m=label['ell'],label['b'],label['m']; assert label['orientation']==0
        expected={}
        if family=='leave-motif-anchor':
            j,h=label['anchor_core'],label['anchor_excluded']
            for vid,(rep,o) in enumerate(zip(reps,sizes),1):
                value=0
                for i,anchor in enumerate(rep):
                    other=rep[:i]+rep[i+1:]
                    ways=interval_count_dp(other,anchor,n,ell,b,j,h)
                    if ways:value+=ways*residual_coefficient(other,m,ell,0)
                if value:expected[vid]=o*value
        else:
            q={'F':0,'N*F':1,'t*F':2,'S+5t2-12t':3,
               '24S+t(t-3)(t-5)(5t-26)':4}[label['polynomial']]
            value=graph_values(m,())[q]*V*comb(n,ell)*comb(n-ell,b)
            if value:expected[0]=value
            for vid,(rep,o) in enumerate(zip(reps,sizes),1):
                common=set.intersection(*(set(u) for u in rep));union=set().union(*rep)
                ways=choose(len(common),ell)*choose(n-len(union),b)
                if ways:
                    value=residual_coefficient(rep,m,ell,q)
                    if value:expected[vid]=o*ways*value
        divisor=0
        for value in expected.values():divisor=gcd(divisor,value)
        assert divisor>0
        expected={v:c//divisor for v,c in expected.items()}
        actual={v:c for v,i,j,c in block['entries']}
        assert expected==actual,(idx,expected,actual)
        assert divisor==label['removed_gcd'],idx
        results.append(dict(index=idx,label=label,nonzero=len(expected),status='PASS'))
        detail('leave inequality',idx,'PASS')
    assert len(results)==23
    return results

