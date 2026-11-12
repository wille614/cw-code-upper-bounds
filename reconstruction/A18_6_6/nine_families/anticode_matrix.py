"""Independent reconstruction of SPEC01 section 7, all 15 blocks.

This program imports labels and independently generated orbit representatives,
but never opens the model or its coefficients.  It reconstructs full matrices
before checking the retained-row and scale metadata.

For a fixed compatible support S, count ordered roles (a,b,c), and for the
negative term (a,b,c,q), whose union is exactly S.  Sum over all coordinate
sets R of size 4+2r with |a & R| = s.  Multiplication by

    |Orb(S)| / (18564 * C(6,s) * C(12,4+2r-s))

converts this count to the coefficient at one fixed (a,R): the coordinate
permutation group acts transitively on these pairs.  This is an independently
derived extension of the special r=s=0 double count in the specification.

Positive terms need |S| <= 3.  R's with equal selected counts in each
coordinate-pattern cell are bundled with their exact binomial multiplicity.
In a negative term, q is distinct from a,b,c since those three are outside
the anticode.  Once q is fixed, there are only 15,73,66 choices of R for
r=0,1,2 respectively.  Thus all four-word orbits are completely enumerated
without enumerating enormous ordered pairs of ambient words.
"""

from collections import defaultdict
from functools import lru_cache
from itertools import combinations, product
from math import comb, gcd
import json
import time

from common import FULL, V, load_labels, load_orbits, clean_entries, write_rebuilt

FAMILY = 'anticode-vacancy-matrix'


def candidate_rows(r, s):
    """All rows from the defining inequalities, independently of label rows."""
    t = 4 + 2*r
    caps = (12-t+s, t-s, 6-s, s)
    rows = []
    for alpha in product(*(range(cap+1) for cap in caps)):
        if sum(alpha) != 6:
            continue
        if alpha[1] + alpha[3] >= 4+r:
            continue
        overlap = alpha[2] + alpha[3]
        if overlap != 6 and overlap > 3:
            continue
        rows.append(alpha)
    return caps, rows


@lru_cache(None)
def onto_triples(indices):
    desired = set(indices)
    return tuple(triple for triple in product(indices, repeat=3)
                 if set(triple) == desired)


def weighted_R(histogram, total):
    """Yield one R per cell-count vector and its actual-coordinate multiplicity."""
    active = []
    offset = 0
    for cap in histogram:
        if cap:
            active.append((cap, offset))
        offset += cap
    assert offset == 18
    suffix = [0] * (len(active)+1)
    for i in range(len(active)-1, -1, -1):
        suffix[i] = suffix[i+1] + active[i][0]

    def visit(i, remaining, R, weight):
        if i == len(active):
            if remaining == 0:
                yield R, weight
            return
        cap, offset = active[i]
        lo = max(0, remaining - suffix[i+1])
        hi = min(cap, remaining)
        for t in range(lo, hi+1):
            segment = ((1 << t)-1) << offset
            yield from visit(i+1, remaining-t, R | segment,
                             weight * comb(cap, t))
    yield from visit(0, total, 0, 1)


@lru_cache(None)
def anticode_R(q, r):
    """All R such that q belongs to A_R, with no model-dependent choices."""
    inside = tuple(1 << x for x in range(18) if q & (1 << x))
    outside = tuple(1 << x for x in range(18) if not q & (1 << x))
    assert len(inside) == 6 and len(outside) == 12
    total = 4 + 2*r
    answer = []
    for in_count in range(4+r, min(6, total)+1):
        out_count = total-in_count
        for ci in combinations(inside, in_count):
            for co in combinations(outside, out_count):
                answer.append(sum(ci)+sum(co))
    assert len(answer) == (15, 73, 66)[r]
    assert len(set(answer)) == len(answer)
    return tuple(answer)


def row_class(a, b, R):
    in_a = (a & b).bit_count()
    in_R = (R & b).bit_count()
    both = (a & b & R).bit_count()
    return (6-in_a-in_R+both, in_R-both, in_a-both, both)


def count_roles(words, triples, R, r, weight, infos, local):
    """Add each ordered role choice once, restricted to stored upper triangle."""
    R_intersections = tuple((u & R).bit_count() for u in words)
    cache = {}
    for a, b, c in triples:
        s = R_intersections[a]
        if s >= 4+r:
            continue
        if R_intersections[b] >= 4+r or R_intersections[c] >= 4+r:
            continue
        info = infos[r, s]
        row_index = info['row_index']
        for x in (b, c):
            if (a, x) not in cache:
                alpha = row_class(words[a], words[x], R)
                # For these compatible words the defining conditions imply
                # that this row belongs to the full, independently built set.
                cache[a, x] = row_index[alpha]
        i, j = cache[a, b], cache[a, c]
        if i <= j:
            local[r, s][i, j] += weight


def main():
    start = time.perf_counter()
    labels = load_labels(FAMILY)
    orbits = load_orbits()
    infos = {}
    for r in range(3):
        for s in range(4+r):
            caps, rows = candidate_rows(r, s)
            infos[r, s] = {
                'caps': caps, 'rows': rows,
                'row_index': {alpha: i for i, alpha in enumerate(rows)},
                'denominator': V * comb(6, s) * comb(12, 4+2*r-s),
                'entries': {},
            }
    assert len(infos) == 15
    expected_pairs = set(infos)
    found_pairs = {(b['label']['r'], b['label']['overlap']) for b in labels}
    assert found_pairs == expected_pairs and len(labels) == len(infos)
    stats = {'orbit_supports': 0, 'positive_R_profiles': 0,
             'positive_actual_R': 0, 'negative_R_choices': 0,
             'integer_division_checks': 0}
    positive_seconds = 0.0
    negative_seconds = 0.0

    for orbit in orbits:
        v = orbit['v']
        k = orbit['k']
        words = tuple(orbit['words'])
        local = {pair: defaultdict(int) for pair in infos}

        # Positive term: union of three ordered word roles is exactly S.
        t0 = time.perf_counter()
        if k <= 3:
            triples = onto_triples(tuple(range(k)))
            for r in range(3):
                actual_R = 0
                for R, weight in weighted_R(orbit['h'], 4+2*r):
                    stats['positive_R_profiles'] += 1
                    actual_R += weight
                    count_roles(words, triples, R, r, weight, infos, local)
                assert actual_R == comb(18, 4+2*r)
                stats['positive_actual_R'] += actual_R
        positive_seconds += time.perf_counter()-t0

        # Negative term: q belongs to the anticode and therefore cannot be
        # any of a,b,c; their union is precisely S\{q}.
        t0 = time.perf_counter()
        if k >= 2:
            for q_index, q in enumerate(words):
                other = tuple(i for i in range(k) if i != q_index)
                triples = onto_triples(other)
                for r in range(3):
                    for R in anticode_R(q, r):
                        stats['negative_R_choices'] += 1
                        count_roles(words, triples, R, r, -1, infos, local)
        negative_seconds += time.perf_counter()-t0

        for pair, values in local.items():
            info = infos[pair]
            denominator = info['denominator']
            for (i, j), count in values.items():
                if count == 0:
                    continue
                numerator = count * orbit['orbit_size']
                coefficient, remainder = divmod(numerator, denominator)
                assert remainder == 0, (v, pair, i, j, numerator, denominator)
                stats['integer_division_checks'] += 1
                if coefficient:
                    info['entries'][v, i, j] = coefficient
        stats['orbit_supports'] += 1
        if v % 100 == 0 or v == len(orbits):
            print(json.dumps({'family': FAMILY, 'orbits_done': v,
                              'seconds': time.perf_counter()-start}), flush=True)

    records = []
    for block in labels:
        lab = block['label']
        pair = lab['r'], lab['overlap']
        info = infos[pair]
        raw_entries = info['entries']
        g = 0
        live = set()
        for (v, i, j), coefficient in raw_entries.items():
            g = gcd(g, abs(coefficient))
            live.update((i, j))
        assert g > 0
        kept = sorted(live)
        row_map = {old: new for new, old in enumerate(kept)}
        entries = {(v, row_map[i], row_map[j]): coefficient//g
                   for (v, i, j), coefficient in raw_entries.items()}
        rows = [list(info['rows'][i]) for i in kept]
        issues = []
        def check(name, expected, actual):
            if expected != actual:
                issues.append({'field': name, 'model': actual,
                               'reconstructed': expected})
        r, s = pair
        check('caps', list(info['caps']), lab['caps'])
        check('threshold', 4+r, lab['threshold'])
        check('harmonic', [0, 0, 0, 0], lab['harmonic'])
        check('rows', rows, lab['rows'])
        check('size', len(kept), block['size'])
        check('positive_scale_removed', g, lab['positive_scale_removed'])
        a, R = lab['anchor'], lab['subset']
        check('anchor_weight', 6, a.bit_count())
        check('subset_size', 4+2*r, R.bit_count())
        check('overlap', s, (a & R).bit_count())
        check('coordinate_mask_validity', 0, (a | R) & ~FULL)
        actual_caps = [(FULL ^ (a | R)).bit_count(),
                       (R & ~a).bit_count(), (a & ~R).bit_count(),
                       (a & R).bit_count()]
        check('caps_from_masks', list(info['caps']), actual_caps)
        records.append({
            'block_id': block['block_id'], 'size': len(kept),
            'entries': clean_entries(entries),
            'computed_metadata': {
                'positive_scale_removed': g,
                'candidate_rows': [list(row) for row in info['rows']],
                'kept_rows': kept, 'rows': rows,
                'pair_count_denominator': info['denominator'],
                'nonzero_coefficients': len(entries),
            },
            'metadata_issues': issues,
        })
    stats['positive_seconds'] = positive_seconds
    stats['negative_seconds'] = negative_seconds
    elapsed = time.perf_counter()-start
    path = write_rebuilt(FAMILY, records, elapsed, statistics=stats,
                         reconstruction_method='ordered-role orbit double count')
    print(json.dumps({'family': FAMILY, 'seconds': elapsed,
                      'output': str(path), 'statistics': stats}), flush=True)


if __name__ == '__main__':
    main()
