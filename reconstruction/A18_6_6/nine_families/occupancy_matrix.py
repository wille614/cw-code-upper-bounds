#!/usr/bin/env python3
"""Rebuild all 52 SPEC01 section 10 blocks using exact integer arithmetic.

The C++ helper exhausts coordinate subsets for each orbit and computes the
Möbius coefficients for the basis polynomials 1,N,N². It never reads entries.
The small fixed-support counts use int64, with an a priori bound 199,147,520.
Orbit-size multiplication and all final coefficients use unbounded Python ints.
"""
from __future__ import annotations
import json
import subprocess
import time
from math import gcd
from functools import reduce
from common import ROOT, DATA, V, J, load_labels, load_orbits, write_rebuilt
from coordinate_spectrum import metadata_check


def rebuild():
    start = time.perf_counter()
    orbits = load_orbits()
    labels = load_labels('occupancy-matrix')
    orbit_sizes = {r['v']:r['orbit_size'] for r in orbits}
    orbit_sizes[0] = V
    inp = DATA/'occupancy_orbit_input.txt'
    inp.write_text(str(len(orbits))+'\n'+''.join(
        f"{r['v']} {r['k']} "+' '.join(map(str,r['words']))+'\n' for r in orbits))
    exe = DATA/'occupancy_basis'
    compile_start = time.perf_counter()
    subprocess.run(['g++','-std=c++17','-O3','-mpopcnt',str(ROOT/'occupancy_basis.cpp'),'-o',str(exe)],check=True)
    compile_seconds = time.perf_counter()-compile_start
    raw_path, timing_path = DATA/'occupancy_basis.tsv', DATA/'occupancy_basis_timings.tsv'
    count_start = time.perf_counter()
    subprocess.run([str(exe),str(inp),str(raw_path),str(timing_path)],check=True)
    count_seconds = time.perf_counter()-count_start
    # Lift fixed representative coefficients by the complete orbit size.
    bases = {}
    for line in raw_path.read_text().splitlines():
        core, comp, v, i, j, x, y, z = map(int,line.split())
        scale = orbit_sizes[v]
        bases.setdefault((core,bool(comp)),{})[(v,i,j)] = (x*scale,y*scale,z*scale)
    geometry_timings = []
    for line in timing_path.read_text().splitlines():
        s,t,nt,pairs = line.split()
        geometry_timings.append({'core':int(s),'seconds':float(t),
                                 'all_coordinate_subsets':int(nt),
                                 'orbit_subset_pairs':int(pairs)})
    records = []
    for item in labels:
        tick = time.perf_counter()
        lab = item['label']
        core, comp = lab['core'],lab['complement']
        weight = 12 if comp else 6
        cap = J(18-core,weight-core)
        classes = list(range(max(0,weight-(18-core)),min(core,weight)+1))
        rows = len(classes)+1
        issues=[]
        metadata_check(issues,'capacity',lab['capacity'],cap)
        metadata_check(issues,'classes',lab['classes'],classes)
        kind, k = lab['kind']
        if kind=='gap':
            if not 0 <= k < cap:
                issues.append({'location':'kind','model_label':lab['kind'],
                               'independent':'gap index must lie in 0..J-1'})
            polynomial = (k*(k+1),-2*k-1,1)
        elif kind=='capacity':
            metadata_check(issues,'kind[1]',k,cap)
            polynomial = (0,cap,-1)
        else:
            raise ValueError(f"unknown kind {kind!r}")
        metadata_check(issues,'polynomial',lab['polynomial'],list(polynomial))
        if min(sum(a*x**i for i,a in enumerate(polynomial)) for x in range(cap+1))<0:
            issues.append({'location':'nonnegativity','model_label':lab['kind'],
                           'independent':'negative on integer 0..J'})
        raw = {}
        for key,basis in bases[(core,comp)].items():
            value = sum(x*y for x,y in zip(polynomial,basis))
            if value:
                raw[key] = value
        scale = reduce(gcd,(abs(c) for c in raw.values()),0)
        if not scale:
            raise ValueError(f"identically zero matrix {item['block_id']}")
        kept = sorted({i for v,i,j in raw}|{j for v,i,j in raw})
        remap = {i:new for new,i in enumerate(kept)}
        entries = [[v,remap[i],remap[j],c//scale] for (v,i,j),c in sorted(raw.items())]
        metadata_check(issues,'positive_scale_removed',lab['positive_scale_removed'],scale)
        metadata_check(issues,'kept_rows',lab['kept_rows'],kept)
        metadata_check(issues,'size',item['size'],len(kept))
        records.append({'block_id':item['block_id'],'size':len(kept),'entries':entries,
                        'computed_metadata':{'capacity':cap,'classes':classes,'polynomial':list(polynomial),
                                             'positive_scale_removed':scale,'kept_rows':kept,
                                             'uncompressed_size':rows},
                        'metadata_issues':issues,'seconds':time.perf_counter()-tick})
    elapsed = time.perf_counter()-start
    path = write_rebuilt('occupancy-matrix',records,elapsed,
                         compile_seconds=compile_seconds,count_seconds=count_seconds,
                         geometry_timings=geometry_timings,
                         arithmetic='C++ int64 bounded fixed-support counts, then Python arbitrary precision integers')
    print(json.dumps({'family':'occupancy-matrix','blocks':len(records),'seconds':elapsed,
                      'coefficient_keys':sum(len(r['entries']) for r in records),
                      'metadata_issues':sum(len(r['metadata_issues']) for r in records),
                      'output':str(path)},ensure_ascii=False),flush=True)

if __name__=='__main__':
    rebuild()
