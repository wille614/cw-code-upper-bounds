#!/usr/bin/env python3
"""One entry point for the exact checks of the second table update."""
import argparse,json,sys,time
from pathlib import Path
from checks import basic

ROOT=Path(__file__).resolve().parent

def name(e):return 'A({n},{d},{w})'.format(**e)
def check(e,entries,detail):
    kind=e['kind'];n,d,w=(e[k] for k in ('n','d','w'))
    if kind=='code' or n==64:
        words=basic.read_code(ROOT,e['code'],detail)
        if n==64:basic.translate64(words,detail)
        elif n in (30,32):basic.ternary(ROOT,n,words,detail)
        return 'code size, length, weights, every pair'+('; 128 translates and Johnson' if n==64 else '')
    if n==58:
        basic.plane58(detail);return 'projective plane and configuration arithmetic; cited nonexistence theorem'
    if kind=='packing':
        assert basic.packing(n,w,w-d//2+1,detail)==e['bound'];return 'all divisibility hypotheses and leave-bound arithmetic'
    if kind=='jpf':
        assert w-d//2+1==3 and basic.jpf(n,w,detail)==e['bound'];return 'Theorem 13 hypotheses and nested-integer arithmetic'
    if kind=='johnson':
        parent=next((x for x in entries if [x[k] for k in ('n','d','w')]==e['parent']),None)
        if parent:check(parent,entries,detail)
        else:detail('Tabulated premise:',e['parent'],'<=',e['premise'],'; source: Brouwer table, 2026-11-10')
        assert e['parent']==[n-1,d,w-1] and n*e['premise']//w==e['bound']
        detail('floor(',n,'*',e['premise'],'/',w,') =',e['bound'])
        return 'Johnson arithmetic'+(' and parent evidence' if parent else ' from cited table premise')
    if kind=='lp':
        basic.linear(ROOT,e,detail);return 'exact rational dual and every required shell bound'
    if kind=='projection':
        from checks import projection_first,projection_second
        projection_first.check(ROOT,detail);projection_second.check(ROOT,detail)
        detail('Both independent certificate sets exclude all 102 graphs: 72 Gram, 22 single-column, 8 paired-column.')
        detail('Integer upper bound: 26; the written reduction excludes larger codes by a 28-word subset.')
        return 'two certificate sets; complete 102-graph coverage; all integer products'
    if kind=='sdp':
        from checks import semidefinite
        if n==18:semidefinite.matrix18(ROOT,detail)
        else:semidefinite.gram21(ROOT,detail)
        return 'exact rational dual, positive semidefiniteness and all residuals'+('; all 23 leave rows' if n==18 else '; full-to-principal coefficient equality')
    raise ValueError('Unknown evidence type')

def main():
    if not __debug__:raise SystemExit('Assertions must be enabled: omit -O and -OO.')
    if sys.version_info<(3,11):raise SystemExit('Python 3.11 or later is required.')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--all',action='store_true',help='include the two large semidefinite certificates')
    parser.add_argument('--list',action='store_true',help='list entries and evidence files')
    parser.add_argument('parameters',nargs='*',type=int,metavar='n d w')
    args=parser.parse_args();entries=json.loads((ROOT/'data/entries.json').read_text(encoding='utf-8'))
    if args.parameters:
        if len(args.parameters)!=3 or args.all or args.list:parser.error('supply exactly n d w, without another option')
        selected=[e for e in entries if [e[k] for k in ('n','d','w')]==args.parameters]
        if not selected:parser.error('this entry is not in the update')
    else:selected=entries
    if args.list:
        if args.all:parser.error('--list and --all cannot be combined')
        for e in selected:print(name(e),'|','; '.join(e['evidence']),flush=True)
        print('TOTAL:',len(selected),'entries; the two large certificates also require all shards named in their model indexes.')
        return
    if args.all:
        estimates=json.loads((ROOT/'data/measurements.json').read_text(encoding='utf-8'))
        print('Including the two large certificates. Approximate requirement:',estimates['estimated_all_seconds'],'s;',estimates['estimated_all_peak_bytes'],'bytes of memory. Actual resource use depends on the machine.',flush=True)
    start=time.perf_counter();passed=0;skipped=0
    detail=print if args.parameters else lambda *a,**k:None
    for e in selected:
        if e.get('large') and not(args.all or args.parameters):
            print(name(e),'|',e['old'],'->',e['bound'],'| large certificate | SKIPPED | 0.000 s',flush=True);skipped+=1;continue
        t=time.perf_counter()
        try:
            what=check(e,entries,detail)
        except Exception as exc:
            print(name(e),'|',e['old'],'->',e['bound'],'| FAILED |',str(exc),flush=True)
            raise
        passed+=1
        print(name(e),'|',e['old'],'->',e['bound'],'|',what,'| PASS |',f'{time.perf_counter()-t:.3f} s',flush=True)
    print('SUMMARY:',passed,'PASS,',skipped,'SKIPPED, 0 FAILED;',f'{time.perf_counter()-start:.3f} s',flush=True)
    if skipped:print('For the two large certificates: python verify.py --all (or python verify.py 21 8 9 / python verify.py 18 6 6).',flush=True)

if __name__=='__main__':main()
