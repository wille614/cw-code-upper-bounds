#!/usr/bin/env python3
"""Independent exact reconstruction of SPEC01 sections 1, 2, 4 (D2).

No model coefficient is passed to reconstruction.  Young vectors are derived
from distinct top-row permutations followed by signed column swaps.  All
arithmetic affecting the result is Python integer or fractions.Fraction.
Timing uses a floating point clock only; it never enters a coefficient.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
from itertools import combinations, product
import json
from math import comb, factorial, gcd, prod
from pathlib import Path
import time


@lru_cache(maxsize=None)
def local_pair(c: int, h: int, a: int, b: int):
    """Joint-bit histogram polynomial from explicit two-row Young vectors.

    Each of the h two-cell columns contributes 2*(Y00*Y11-Y01*Y10).
    The c-2h unpaired cells carry a-h and b-h ones.  The multinomial
    coefficient counts their ordered placements.  q counts joint ones.
    """
    assert 0 <= h <= c // 2
    assert h <= a <= c-h and h <= b <= c-h
    n, A, B = c-2*h, a-h, b-h
    out = Counter()
    for r in range(max(0, A+B-n), min(A, B)+1):
        free = (n-A-B+r, B-r, A-r, r)
        mult = factorial(n) // prod(factorial(x) for x in free)
        for s in range(h+1):
            q = r+h-s
            out[q] += (1 << h) * (-1)**s * comb(h, s) * mult
    return tuple(((c-a-b+q, b-q, a-q, q), value)
                 for q, value in sorted(out.items()) if value)


def explicit_young(c: int, h: int, a: int):
    """Definition-level local u: distinct row fillings and signed swaps."""
    ans = Counter()
    top = c-h
    bottom = sum(1 << (top+j) for j in range(h))
    for ones in combinations(range(top), a-h):
        initial = bottom | sum(1 << j for j in ones)
        for switches in range(1 << h):
            x = initial
            for j in range(h):
                if switches & (1 << j):
                    k = top+j
                    if ((x >> j) & 1) != ((x >> k) & 1):
                        x ^= (1 << j) | (1 << k)
            ans[x] += (-1)**switches.bit_count()
    return {x: v for x, v in ans.items() if v}


def selfcheck_local():
    """Exhaustive for every c=0..9 used or potentially used by D2."""
    started = time.perf_counter()
    cases = vectors = pairs = 0
    for c in range(10):
        for h in range(c//2+1):
            us = {a: explicit_young(c, h, a) for a in range(h, c-h+1)}
            vectors += len(us)
            for a, ux in us.items():
                for b, uy in us.items():
                    actual = Counter()
                    for x, vx in ux.items():
                        for y, vy in uy.items():
                            q = (x & y).bit_count()
                            hist = (c-a-b+q, b-q, a-q, q)
                            actual[hist] += vx*vy
                            pairs += 1
                    actual = {key: value for key, value in actual.items() if value}
                    expected = dict(local_pair(c, h, a, b))
                    assert actual == expected, (c, h, a, b, actual, expected)
                    cases += 1
    return {"status": "PASS", "coordinate_group_sizes": list(range(10)),
            "explicit_young_vectors": vectors, "local_table_pairs": cases,
            "explicit_ordered_string_pairs": pairs,
            "elapsed_seconds": time.perf_counter()-started}


def independent_catalog():
    """Generate all nonempty (anchor,harmonic,rows) directly from constraints."""
    out = []
    for t in range(3, 7):
        c = (12-t, t, t, 6-t)
        for h in product(*(range(x//2+1) for x in c)):
            rows = []
            for alpha in product(*(range(s, x-s+1) for s, x in zip(h, c))):
                if sum(alpha) != 6:
                    continue
                if alpha[2]+alpha[3] not in (0, 1, 2, 3, 6):
                    continue
                if alpha[1]+alpha[3] not in (0, 1, 2, 3, 6):
                    continue
                rows.append(alpha)
            if rows:
                out.append((c, h, rows))
    assert len(out) == 204
    return out


def load_lookup(path: Path):
    lookup = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        row = tuple(map(int, line.split()))
        assert len(row) == 17
        key, value = row[:16], row[16]
        assert key not in lookup
        lookup[key] = value
    return lookup


def subset_cardinality(variable: int):
    if variable == 1:
        return 1
    if variable <= 5:
        return 2
    if variable <= 40:
        return 3
    return 4


def rational_gram_scan(gram):
    """Exact PSD elimination in original left-to-right column order."""
    n = len(gram)
    a = [[Fraction(value) for value in row] for row in gram]
    kept, pivots = [], []
    for j in range(n):
        pivot = a[j][j]
        assert pivot >= 0, ("non-PSD Gram", j, pivot)
        if not pivot:
            assert all(a[k][j] == 0 for k in range(j, n))
            continue
        kept.append(j)
        pivots.append([pivot.numerator, pivot.denominator])
        for k in range(j+1, n):
            for l in range(k, n):
                a[k][l] -= a[k][j]*a[j][l]/pivot
                a[l][k] = a[k][l]
    return kept, pivots


def reconstruct(c, h, rows, lookup):
    """Return the full compressed integer block, without reading model values."""
    started = time.perf_counter()
    n = len(rows)
    anchor_a = (0, 0, c[2], c[3])
    anchor_b = (0, c[1], 0, c[3])
    cardinalities = [2 if alpha in (anchor_a, anchor_b) else 3 for alpha in rows]
    raw = {}
    gram = [[0]*n for _ in range(n)]
    unmerged_terms = invalid_terms = 0
    max_abs_raw = 0
    g = 0
    for i, alpha in enumerate(rows):
        for j in range(i, n):
            beta = rows[j]
            parts = [local_pair(cp, hp, ap, bp)
                     for cp, hp, ap, bp in zip(c, h, alpha, beta)]
            co = Counter()
            for term_tuple in product(*parts):
                unmerged_terms += 1
                hist = tuple(x for item, _ in term_tuple for x in item)
                # All other pairs are already known compatible from rows.
                intersection = sum(hist[p] for p in (3, 7, 11, 15))
                if intersection not in (0, 1, 2, 3, 6):
                    invalid_terms += 1
                    continue
                v = lookup.get(hist)
                assert v is not None, ("missing compatible orbit", hist)
                co[v] += prod(value for _, value in term_tuple)
            co = {v: value for v, value in co.items() if value}
            raw[i, j] = co
            for value in co.values():
                g = gcd(g, abs(value))
                max_abs_raw = max(max_abs_raw, abs(value))
            if cardinalities[i] == cardinalities[j]:
                gram[i][j] = gram[j][i] = sum(
                    value for v, value in co.items()
                    if subset_cardinality(v) == cardinalities[i])

    assert g > 0
    # A second Gram computation works directly with the supports of P.
    # Weight classes are disjoint except x=a and x=b map to the same D.
    direct_gram = [[0]*n for _ in range(n)]
    for i, alpha in enumerate(rows):
        direct_gram[i][i] = prod((1 << hp)*comb(cp-2*hp, ap-hp)
                                for cp, hp, ap in zip(c, h, alpha))
        if cardinalities[i] == 2:
            for j in range(n):
                if cardinalities[j] == 2:
                    direct_gram[i][j] = 1
    assert gram == direct_gram, ("P Gram check", c, h)
    kept, pivots = rational_gram_scan(gram)
    entries = []
    for ii, i in enumerate(kept):
        for jj in range(ii, len(kept)):
            j = kept[jj]
            for v, value in raw[i, j].items():
                assert value % g == 0
                entries.append([v, ii, jj, value//g])
    entries.sort()
    digest = hashlib.sha256(json.dumps(entries, separators=(",", ":")).encode()).hexdigest()
    details = {"original_columns": n, "kept_columns": kept,
               "stored_columns": len(kept), "row_cardinalities": cardinalities,
               "positive_scale_removed": g,
               "raw_nonzero_keys": sum(len(x) for x in raw.values()),
               "reconstructed_nonzero_keys": len(entries),
               "joint_histogram_product_terms": unmerged_terms,
               "incompatible_product_terms_discarded": invalid_terms,
               "max_absolute_raw_coefficient": max_abs_raw,
               "gram_from_coefficients_equals_explicit_support_gram": True,
               "gram": gram, "exact_gram_pivots": pivots,
               "reconstructed_entries_sha256": digest,
               "reconstruction_seconds": time.perf_counter()-started}
    return entries, details


def compare(entries, block, details, c, h, rows):
    started = time.perf_counter()
    expected_label = {"family": "D2", "anchor": list(c), "harmonic": list(h),
                      "rows": [list(x) for x in rows],
                      "positive_scale_removed": details["positive_scale_removed"]}
    issues = []
    if block["label"] != expected_label:
        issues.append({"kind": "label", "model": block["label"], "computed": expected_label})
    if block["size"] != details["stored_columns"]:
        issues.append({"kind": "size", "model": block["size"],
                       "computed": details["stored_columns"]})
    model = {}
    for v, i, j, value in block["entries"]:
        key = (v, i, j)
        if key in model:
            issues.append({"kind": "duplicate_model_key", "key": list(key)})
        if value == 0:
            issues.append({"kind": "explicit_zero_model_key", "key": list(key)})
        if not (0 <= i <= j < block["size"] and 0 <= v <= 941):
            issues.append({"kind": "invalid_model_key", "key": list(key)})
        model[key] = value
    actual = {(v, i, j): value for v, i, j, value in entries}
    for key in sorted(set(model) | set(actual)):
        if model.get(key, 0) != actual.get(key, 0):
            issues.append({"kind": "coefficient", "key": list(key),
                           "model": model.get(key, 0), "computed": actual.get(key, 0)})
    details.update({"status": "PASS" if not issues else "FAIL",
                    "model_nonzero_keys": len(model),
                    "key_sets_equal": model.keys() == actual.keys(),
                    "all_key_values_equal": model == actual,
                    "issues": issues,
                    "comparison_seconds": time.perf_counter()-started})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path,
                        default=Path("upload/A18_6_6_model.json(1).gz"))
    parser.add_argument("--lookup", type=Path,
                        default=Path("sem02b/results/orbit_lookup.tsv"))
    parser.add_argument("--output", type=Path, default=Path("sem02b/results"))
    parser.add_argument("--blocks", help="comma-separated global zero-based indices")
    parser.add_argument("--selfcheck-only", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    overall = time.perf_counter()
    selfcheck = selfcheck_local()
    print(json.dumps({"local_young_selfcheck": selfcheck}), flush=True)
    if args.selfcheck_only:
        (args.output/"d2_local_selfcheck.json").write_text(
            json.dumps(selfcheck, indent=2)+"\n")
        return
    lookup = load_lookup(args.lookup)
    catalog = independent_catalog()
    # The JSON parser necessarily loads the bytes of the whole model.
    # Reconstruct does not receive any model object or model coefficient.
    with gzip.open(args.model, "rt") as fh:
        model = json.load(fh)
    selected = set(map(int, args.blocks.split(","))) if args.blocks else set(range(65, 269))
    assert selected <= set(range(65, 269))
    sequence = sorted(selected, key=lambda b: (len(catalog[b-65][2]), b))
    reports, reconstructed = [], []
    log_path = args.output/("d2_blocks_partial.jsonl" if args.blocks else "d2_blocks.jsonl")
    with log_path.open("w") as log:
        for b in sequence:
            c, h, rows = catalog[b-65]
            entries, details = reconstruct(c, h, rows, lookup)
            # This is the first point at which the model's entries for b are used.
            compare(entries, model["blocks"][b], details, c, h, rows)
            details.update({"block": b, "anchor": list(c), "harmonic": list(h)})
            reports.append(details)
            reconstructed.append({"block": b, "entries": entries})
            log.write(json.dumps(details, separators=(",", ":"))+"\n")
            log.flush()
            print(json.dumps({k: details[k] for k in
                              ("block", "status", "original_columns", "stored_columns",
                               "reconstructed_nonzero_keys", "positive_scale_removed",
                               "reconstruction_seconds")}), flush=True)
    reports.sort(key=lambda x: x["block"])
    result = {"family": "D2", "method": "explicit Young-vector pair factorization",
              "coefficient_arithmetic": "unbounded Python int",
              "rank_arithmetic": "fractions.Fraction",
              "local_selfcheck": selfcheck,
              "blocks_in_scope": 204, "fully_reconstructed": len(reports),
              "fully_equal": sum(r["status"] == "PASS" for r in reports),
              "mismatching": sum(r["status"] != "PASS" for r in reports),
              "not_reconstructed": 204-len(reports),
              "sample_only_blocks": [], "run_order": sequence,
              "original_columns": sum(r["original_columns"] for r in reports),
              "kept_columns": sum(r["stored_columns"] for r in reports),
              "reconstructed_nonzero_keys": sum(r["reconstructed_nonzero_keys"] for r in reports),
              "reconstruction_seconds": sum(r["reconstruction_seconds"] for r in reports),
              "comparison_seconds": sum(r["comparison_seconds"] for r in reports),
              "overall_seconds": time.perf_counter()-overall,
              "blocks": reports}
    suffix = "_partial" if args.blocks else ""
    (args.output/f"d2_summary{suffix}.json").write_text(json.dumps(result, indent=2)+"\n")
    with gzip.open(args.output/f"d2_reconstructed{suffix}.json.gz", "wt") as fh:
        json.dump(sorted(reconstructed, key=lambda x: x["block"]), fh, separators=(",", ":"))
    print(json.dumps({k: v for k, v in result.items() if k not in ("blocks", "run_order")}), flush=True)
    if result["mismatching"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
