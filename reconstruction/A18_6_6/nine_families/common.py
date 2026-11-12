"""Independent exact combinatorial utilities for SEM02a.

No model coefficients or generator examples are imported by this module.
"""
from __future__ import annotations
from collections import defaultdict
from functools import lru_cache
from itertools import permutations, combinations
from math import factorial, comb, gcd
from pathlib import Path
import json
import time

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
N, W, V, BOUND = 18, 6, comb(18, 6), 186
FULL = (1 << N) - 1


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


@lru_cache(None)
def perm_maps(k):
    """Maps histogram positions under all permutations of word roles."""
    out = []
    for perm in permutations(range(k)):
        inv = [0] * (1 << k)
        for p in range(1 << k):
            q = 0
            for old in perm:
                q = (q << 1) | ((p >> (k - 1 - old)) & 1)
            inv[q] = p
        out.append(tuple(inv))
    return tuple(out)


def canonical(h):
    h = tuple(h)
    k = len(h).bit_length() - 1
    return min(tuple(h[p] for p in idx) for idx in perm_maps(k))


def restore_words(h):
    """Use consecutive coordinate segments, as in SPEC01 section 1."""
    k = len(h).bit_length() - 1
    words = [0] * k
    offset = 0
    for p, count in enumerate(h):
        segment = ((1 << count) - 1) << offset
        offset += count
        for j in range(k):
            if p & (1 << (k - 1 - j)):
                words[j] |= segment
    assert offset == N
    return tuple(words)


def histogram(words):
    words = tuple(dict.fromkeys(words))
    cells = [FULL]
    for word in words:
        nxt = []
        for cell in cells:
            inside = cell & word
            nxt.extend((cell ^ inside, inside))
        cells = nxt
    return tuple(cell.bit_count() for cell in cells)


def compatible(words):
    words = tuple(set(words))
    return all((a & b).bit_count() <= 3 for a, b in combinations(words, 2))


def orbit_size(h):
    h = tuple(h)
    k = len(h).bit_length() - 1
    stab = sum(all(h[p] == h[idx[p]] for p in range(len(h)))
               for idx in perm_maps(k))
    den = stab
    for a in h:
        den *= factorial(a)
    value, rem = divmod(factorial(N), den)
    assert rem == 0
    return value


def generate_variables():
    """Exhaust all histogram splits when adjoining a new compatible word.

    Every k-set can be obtained from a canonical (k-1)-set after relabeling
    coordinates. Choosing t[p] points in each old cell enumerates all such
    extensions; the bounds on its intersections enforce compatibility.
    """
    layers = [[(N-W, W)]]
    stats = []
    start = time.perf_counter()
    for new_k in range(2, 5):
        unique = set()
        candidates = 0
        t0 = time.perf_counter()
        for h in layers[-1]:
            old_k = new_k - 1
            active = [p for p, c in enumerate(h) if c]
            memberships = [[j for j in range(old_k)
                            if p & (1 << (old_k-1-j))] for p in active]
            suffix = [0] * (len(active)+1)
            for i in reversed(range(len(active))):
                suffix[i] = suffix[i+1] + h[active[i]]
            picked = [0] * len(h)
            intersections = [0] * old_k

            def visit(i, remaining):
                nonlocal candidates
                if i == len(active):
                    if remaining:
                        return
                    raw = tuple(x for p, c in enumerate(h)
                                for x in (c-picked[p], picked[p]))
                    candidates += 1
                    unique.add(canonical(raw))
                    return
                p = active[i]
                lo = max(0, remaining-suffix[i+1])
                hi = min(h[p], remaining)
                for j in memberships[i]:
                    hi = min(hi, 3-intersections[j])
                for a in range(lo, hi+1):
                    picked[p] = a
                    for j in memberships[i]:
                        intersections[j] += a
                    visit(i+1, remaining-a)
                    for j in memberships[i]:
                        intersections[j] -= a
                picked[p] = 0

            visit(0, W)
        layers.append(sorted(unique))
        stats.append({'k': new_k, 'split_candidates': candidates,
                      'orbits': len(unique), 'seconds': time.perf_counter()-t0})
    all_h = [h for layer in layers for h in layer]
    records = []
    for v, h in enumerate(all_h, 1):
        words = restore_words(h)
        assert all(w.bit_count() == W for w in words)
        assert compatible(words) and len(set(words)) == len(words)
        assert histogram(words) == h and canonical(h) == h
        k = len(words)
        stab = sum(tuple(h[p] for p in idx) == h for idx in perm_maps(k))
        records.append({'v': v, 'k': k, 'h': list(h), 'words': list(words),
                        'orbit_size': orbit_size(h), 'stabilizer': stab})
    return records, {'layer_counts': [len(x) for x in layers], 'layers': stats,
                     'generation_seconds': time.perf_counter()-start}


def load_orbits(path=None):
    return json.loads(Path(path or DATA/'orbits.json').read_text())


class OrbitIndex:
    def __init__(self, records=None):
        self.records = records if records is not None else load_orbits()
        self.by_hist = {}
        self.by_v = {r['v']: r for r in self.records}
        for r in self.records:
            h, k = r['h'], r['k']
            for idx in perm_maps(k):
                raw = tuple(h[p] for p in idx)
                prior = self.by_hist.setdefault(raw, r['v'])
                assert prior == r['v']

    def v(self, words):
        """Return 0 for the empty support, None for an incompatible support."""
        words = tuple(dict.fromkeys(words))
        if not words:
            return 0
        if not compatible(words):
            return None
        return self.by_hist[histogram(words)]


@lru_cache(None)
def J(m, k):
    if k < 0 or k > m:
        return 0
    if min(k, m-k) < 3:
        return 1
    return min(m*J(m-1, k-1)//k, m*J(m-1, k)//(m-k))


def load_labels(family=None):
    labels = json.loads((DATA/'labels.json').read_text())
    return [r for r in labels if family is None or r['label']['family'] == family]


def clean_entries(entries):
    """Convert dict keyed by (v,i,j) to canonical nonzero integer records."""
    return [[v, i, j, c] for (v, i, j), c in sorted(entries.items()) if c]


def gcd_values(entries):
    ans = 0
    for c in entries.values():
        ans = gcd(ans, abs(c))
    return ans


def write_rebuilt(family, records, elapsed, **extra):
    path = DATA / ('rebuilt_' + family + '.json')
    payload = {'family': family, 'seconds': elapsed, 'blocks': records, **extra}
    path.write_text(json.dumps(payload, ensure_ascii=False, separators=(',',':')))
    return path


if __name__ == '__main__':
    DATA.mkdir(exist_ok=True)
    records, stats = generate_variables()
    (DATA/'orbits.json').write_text(json.dumps(records, separators=(',',':')))
    (DATA/'orbit_generation.json').write_text(json.dumps(stats, indent=2))
    print(json.dumps(stats, indent=2), flush=True)
