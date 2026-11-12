#!/usr/bin/env python3
"""Reproduce the independent SEM02b verification from bundled inputs."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
CODE=ROOT/'sem02b'
RESULTS=CODE/'results'
MODEL=ROOT/'inputs'/'A18_6_6_model.json.gz'
EXPECTED='157bd486c563c49bdb65cec5b00e245e8f816b10ab67d721841c34b3c5aabe77'


def run(script,*args):
    command=[sys.executable,str(CODE/script),*map(str,args)]
    print('RUN',script,*map(str,args),flush=True)
    subprocess.run(command,cwd=ROOT,check=True)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--workers',type=int,default=4)
    p.add_argument('--check-recorded',action='store_true',
                   help='Only re-read and compare the supplied output files; do not reconstruct.')
    args=p.parse_args()
    assert hashlib.sha256(MODEL.read_bytes()).hexdigest()==EXPECTED,'Wrong model file'
    if not args.check_recorded:
        run('orbits.py','--model',MODEL,'--out',RESULTS)
        r=json.loads((RESULTS/'orbits.json').read_text())
        assert r['variables_exactly_equal'] and r['orbit_counts']==[1,4,35,901]
        run('check_labels.py','--model',MODEL,'--out',RESULTS/'labels.json')
        r=json.loads((RESULTS/'labels.json').read_text())
        assert r['all_tableaux_equal'] and r['blocks_checked']==65
        run('theory_selftest.py')
        run('d0.py','--model',MODEL,'--blocks','all','--workers',args.workers)
        run('check_omitted_shapes.py')
        r=json.loads((RESULTS/'omitted_shapes'/'summary.json').read_text())
        assert r['gram_all_zero_count']==5 and r['all_pairs_evaluated']
        run('d2.py','--model',MODEL,'--lookup',RESULTS/'orbit_lookup.tsv','--output',RESULTS)
    run('audit_results.py','--model',MODEL)
    print('PASS: all 941 variables and all 269 stored D0/D2 blocks (full reconstruction mode).'
          if not args.check_recorded else 'PASS: recorded coefficients of all 269 blocks re-read and compared.',flush=True)

if __name__=='__main__':main()
