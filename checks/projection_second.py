from fractions import Fraction as Q

from itertools import product, combinations

from math import comb

def inverse(m):
    n = len(m)
    a = [[Q(x) for x in row] + [Q(i == j) for j in range(n)] for i, row in enumerate(m)]
    for j in range(n):
        k = next((k for k in range(j, n) if a[k][j]))
        a[j], a[k] = (a[k], a[j])
        d = a[j][j]
        a[j] = [x / d for x in a[j]]
        for i in range(n):
            if i != j:
                c = a[i][j]
                a[i] = [x - c * y for x, y in zip(a[i], a[j])]
    return [r[n:] for r in a]

def dot(a, b):
    return sum((x * y for x, y in zip(a, b)))

def bilinear(a, m, b):
    return sum((a[i] * m[i][j] * b[j] for i in range(len(a)) for j in range(len(b))))

def feature(other, t, v):
    b, a, z, e, p = other
    h = a - t
    return [1, t, h, comb(t, 2), t * h, comb(h, 2)] + list(b) + [t * x for x in b] + [h * x for x in b] + [b[i] * b[j] for i, j in combinations(range(v), 2)]

def verify(case):
    edges, y, expected = CASES[case]
    v = max((max(i, j) for i, j, w in edges), default=-1) + 1
    u = 27 - v
    D = [[0] * v for _ in range(v)]
    for i, j, w in edges:
        D[i][j] = D[j][i] = w
    H = [[32 * (i == j) - 4 * D[i][j] for j in range(v)] for i in range(v)]
    inv = inverse(H)
    beta = [6 * sum(r) for r in inv]
    S = 32 - 36 * (sum(map(sum, inv)) + Q(u, 32))
    assert S > 0
    types = []
    cost = 5 - sum((w for _, _, w in edges))
    for b in product((0, 1), repeat=v):
        z = [1 - 2 * x for x in b]
        for a in range(u + 1):
            r = a + sum(b)
            if (r - 10) * (r - 11) // 2 > cost:
                continue
            e = 1 - dot(beta, z) - Q(3, 16) * (u - 2 * a)
            p = bilinear(z, inv, z) + Q(u, 32) + e * e / S
            if p <= 1:
                types.append((b, a, z, e, p))
    pairs = list(combinations(range(v), 2))
    rhs = [32, 13 * u, 5 * comb(u, 2)] + [13] * v + [5 * u] * v + [5 - D[i][j] for i, j in pairs]
    r0 = len(rhs)
    q = 6 + 3 * v + comb(v, 2)
    assert len(y) == r0 + len(types) * q
    products = []
    for index, (b, a, z, e, p) in enumerate(types):
        single = [1, a, comb(a, 2)] + list(b) + [a * x for x in b] + [b[i] * b[j] for i, j in pairs]
        target = [31, 12 * a, 13 * (u - a), 4 * comb(a, 2), 5 * a * (u - a), 5 * comb(u - a, 2)] + [13 - x for x in b] + [a * (5 - x) for x in b] + [5 * (u - a)] * v + [5 - D[i][j] - b[i] * b[j] for i, j in pairs]
        block = y[r0 + index * q:r0 + (index + 1) * q]
        products.append(dot(y[:r0], single) - dot(block, target))
    paircount = 0
    for i, (b, a, z, e, p) in enumerate(types):
        for j in range(i, len(types)):
            bb, aa, zz, ee, pp = types[j]
            for t in range(max(0, a + aa - u), min(a, aa) + 1):
                off = bilinear(z, inv, zz) + Q(u - 2 * a - 2 * aa + 4 * t, 32) + e * ee / S
                if off * off > (1 - p) * (1 - pp):
                    continue
                paircount += 1
                value = dot(y[r0 + i * q:r0 + (i + 1) * q], feature(types[j], t, v))
                if i != j:
                    value += dot(y[r0 + j * q:r0 + (j + 1) * q], feature(types[i], t, v))
                products.append(value)
    positive = dot(y[:r0], rhs)
    actual = (len(types), paircount, len(y), len(products), max(products), positive)
    assert actual == tuple(expected), (case, actual, expected)
    assert max(products) <= 0 and positive > 0
    detail(case, actual)

from collections import Counter, defaultdict

def fingerdetail(D):
    deg = [sum(row) for row in D]
    return tuple(sorted(((deg[i], tuple(sorted(((D[i][j], deg[j]) for j in range(len(D)) if D[i][j])))) for i in range(len(D)))))

def iso(D, E):
    n = len(D)
    if n != len(E):
        return False
    sig = lambda A, i: (sum(A[i]), tuple(sorted((w for w in A[i] if w))))
    sd = [sig(D, i) for i in range(n)]
    se = [sig(E, i) for i in range(n)]
    if sorted(sd) != sorted(se):
        return False
    cand = [[j for j in range(n) if sd[i] == se[j]] for i in range(n)]
    order = sorted(range(n), key=lambda i: (len(cand[i]), -sd[i][0]))
    mp = {}
    used = set()

    def dfs(k):
        if k == n:
            return True
        i = order[k]
        for j in cand[i]:
            if j in used or any((D[i][a] != E[j][b] for a, b in mp.items())):
                continue
            mp[i] = j
            used.add(j)
            if dfs(k + 1):
                return True
            del mp[i]
            used.remove(j)
        return False
    return dfs(0)

def partitions(m):
    if m == 0:
        yield ()
        return
    x = [0] * m

    def rec(k, mx):
        if k == m:
            yield tuple(x)
            return
        for v in range(mx + 2):
            if k % 2 and x[k - 1] == v:
                continue
            x[k] = v
            yield from rec(k + 1, max(mx, v))
    yield from rec(1, 0)

def graph_from_partition(x):
    n = 1 + max(x, default=-1)
    D = [[0] * n for _ in range(n)]
    for a, b in zip(x[::2], x[1::2]):
        assert a != b
        D[a][b] += 1
        D[b][a] += 1
    return D

def adjacency(edges):
    v = 1 + max((max(a, b) for a, b, w in edges), default=-1)
    D = [[0] * v for _ in range(v)]
    for a, b, w in edges:
        assert a != b and w > 0 and (D[a][b] == 0)
        D[a][b] = D[b][a] = w
    assert all((sum(row) > 0 for row in D))
    assert sum((w for a, b, w in edges)) <= 5
    return D

def check_all_graphs():
    bins = defaultdict(list)
    matrices = {i: adjacency(e) for i, e in GRAPHS.items()}
    assert set(GRAPHS) == set(range(102))
    for i, D in matrices.items():
        bins[sum(map(sum, D)) // 2, fingerdetail(D)].append(i)
    for bucket in bins.values():
        for i, j in combinations(bucket, 2):
            assert not iso(matrices[i], matrices[j]), ('duplicate', i, j)
    expected = [(1, 1), (1, 1), (7, 3), (87, 8), (1657, 23), (43833, 66)]
    seen = set()
    for t in range(6):
        count = 0
        hit = set()
        for partition in partitions(2 * t):
            count += 1
            D = graph_from_partition(partition)
            candidates = bins[t, fingerdetail(D)]
            matches = [i for i in candidates if iso(D, matrices[i])]
            assert len(matches) == 1, ('coverage failure', t, partition, matches)
            hit.add(matches[0])
            seen.add(matches[0])
        assert (count, len(hit)) == expected[t]
        detail('COVERAGE', t, count, len(hit))
    assert seen == set(GRAPHS)

def check_single(case):
    edges = GRAPHS[case]
    D = adjacency(edges)
    v = len(D)
    u = 27 - v
    H = [[32 * (i == j) - 4 * D[i][j] for j in range(v)] for i in range(v)]
    inv = inverse(H)
    assert all((sum((H[i][k] * inv[k][j] for k in range(v))) == (i == j) for i in range(v) for j in range(v)))
    beta = [6 * sum(row) for row in inv]
    S = 32 - 36 * (sum(map(sum, inv)) + Q(u, 32))
    if S < 0:
        assert case not in SINGLE_CERTIFICATES and case not in CASES
        return 'GRAM'
    cols = []
    for b in product((0, 1), repeat=v):
        z = [1 - 2 * x for x in b]
        base = bilinear(z, inv, z) + Q(u, 32)
        for a in range(u + 1):
            r = a + sum(b)
            if (r - 10) * (r - 11) // 2 > 5 - sum((w for _, _, w in edges)):
                continue
            e = 1 - dot(beta, z) - Q(3, 16) * (u - 2 * a)
            if S == 0:
                if e or base > 1:
                    continue
            elif base + e * e / S > 1:
                continue
            cols.append([1, a, comb(a, 2)] + list(b) + [a * x for x in b] + [b[i] * b[j] for i, j in combinations(range(v), 2)])
    assert len(cols) == EXPECTED_TYPES[case]
    if case in CASES:
        assert [list(edge) for edge in edges] == CASES[case][0]
        return 'PAIR'
    y, expected_rhs, expected_max = SINGLE_CERTIFICATES[case]
    rhs = [32, 13 * u, 5 * comb(u, 2)] + [13] * v + [5 * u] * v + [5 - D[i][j] for i, j in combinations(range(v), 2)]
    assert len(y) == len(rhs)
    values = [dot(y, col) for col in cols]
    assert max(values) == expected_max <= 0
    assert dot(y, rhs) == expected_rhs > 0
    return 'SINGLE'


def check(root, detail_output):
    global CASES, GRAPHS, SINGLE_CERTIFICATES, EXPECTED_TYPES, detail
    import json
    data=json.loads((root/'certificates/A32_16_13_second.json').read_text(encoding='utf-8'))
    CASES={int(k):v for k,v in data['CASES'].items()}
    GRAPHS={int(k):v for k,v in data['GRAPHS'].items()}
    SINGLE_CERTIFICATES={int(k):v for k,v in data['SINGLE_CERTIFICATES'].items()}
    EXPECTED_TYPES={int(k):v for k,v in data['EXPECTED_TYPES'].items()}
    detail=detail_output
    check_all_graphs()
    totals=Counter(check_single(k) for k in GRAPHS)
    assert totals=={'GRAM':72,'SINGLE':22,'PAIR':8}
    detail('Second certificate set:',dict(totals))
    for k in CASES:verify(k)
    return totals

