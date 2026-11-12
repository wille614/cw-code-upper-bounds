#!/usr/bin/env python3
"""Rebuild SPEC01 section 9, without opening the model coefficient file.

Inputs are the independently enumerated orbits and the coefficient-free label
snapshot. All coefficient operations are Python integers.
"""
from __future__ import annotations
import json
import time
from math import gcd
from functools import reduce
from common import V, DATA, J, choose, load_labels, load_orbits, write_rebuilt


def metadata_check(issues, field, observed, computed):
    if observed != computed:
        issues.append({'location': field, 'model_label': observed,
                       'independent': computed})


def rebuild():
    start = time.perf_counter()
    orbits = load_orbits()
    labels = load_labels('coordinate-spectrum')
    records = []
    counts = {}
    total_integer_checks = 0
    for item in labels:
        tick = time.perf_counter()
        label = item['label']
        ell, b = label['inside'], label['outside']
        cap = J(18-ell-b, 6-ell)
        roots, sign = label['roots'], label['sign']
        issues = []
        metadata_check(issues, 'capacity', label['capacity'], cap)
        if roots != sorted(set(roots)) or not 1 <= len(roots) <= min(4, cap):
            issues.append({'location': 'roots', 'model_label': roots,
                           'independent': 'must be sorted distinct, degree 1..min(4,J)'})
        if any(not 0 <= r <= cap for r in roots) or sign not in (-1, 1):
            issues.append({'location': 'roots/sign', 'model_label': [roots,sign],
                           'independent': 'roots in 0..J; sign in {-1,+1}'})
        def p(x):
            ans = sign
            for r in roots:
                ans *= x-r
            return ans
        values = [p(x) for x in range(cap+1)]
        total_integer_checks += len(values)
        if min(values) < 0:
            issues.append({'location': 'nonnegativity', 'model_label': [roots,sign],
                           'independent': values})
        degree = len(roots)
        diffs = [sum((-1)**(k-j)*choose(k,j)*p(j) for j in range(k+1))
                 for k in range(degree+1)]
        raw = {}
        constant = p(0)*choose(18,ell)*choose(18-ell,b)*V
        if constant:
            raw[0] = constant
        for orb in orbits:
            k, h = orb['k'], orb['h']
            if k <= degree:
                c = diffs[k]*orb['orbit_size']*choose(h[-1],ell)*choose(h[0],b)
                if c:
                    raw[orb['v']] = c
        scale = reduce(gcd, (abs(c) for c in raw.values()), 0)
        if not scale:
            raise ValueError(f"identically zero scalar block {item['block_id']}")
        metadata_check(issues,'positive_scale_removed',label['positive_scale_removed'],scale)
        metadata_check(issues,'size',item['size'],1)
        entries = [[v,0,0,c//scale] for v,c in sorted(raw.items())]
        records.append({'block_id':item['block_id'],'size':1,'entries':entries,
                        'computed_metadata':{'capacity':cap,'positive_scale_removed':scale,
                                             'finite_differences':diffs},
                        'metadata_issues':issues,'seconds':time.perf_counter()-tick})
        counts[(ell,b,cap)] = counts.get((ell,b,cap),0)+1
    elapsed = time.perf_counter()-start
    path = write_rebuilt('coordinate-spectrum',records,elapsed,
                         arithmetic='Python arbitrary precision integers',
                         positivity_integer_evaluations=total_integer_checks,
                         geometries=[{'inside':l,'outside':b,'capacity':c,'blocks':n}
                                     for (l,b,c),n in sorted(counts.items())])
    print(json.dumps({'family':'coordinate-spectrum','blocks':len(records),'seconds':elapsed,
                      'coefficient_keys':sum(len(r['entries']) for r in records),
                      'metadata_issues':sum(len(r['metadata_issues']) for r in records),
                      'output':str(path)},ensure_ascii=False),flush=True)

if __name__ == '__main__':
    rebuild()
