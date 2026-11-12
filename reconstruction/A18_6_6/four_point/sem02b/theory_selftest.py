#!/usr/bin/env python3
"""Independent exact small-instance check of SPEC §2 against literal Young vectors.
No model or specification examples are loaded. Python integers only.
"""
from collections import Counter
from itertools import permutations, product
from math import comb, factorial
import json, time
from pathlib import Path


def partitions(n, largest=None):
    if n == 0:
        yield ()
        return
    for first in range(min(n, n if largest is None else largest), 0, -1):
        for tail in partitions(n-first, first):
            yield (first,) + tail


def semistandard(shape, m):
    cells = [(r,c) for r,l in enumerate(shape) for c in range(l)]
    values = {}
    def rec(k):
        if k == len(cells):
            yield tuple(tuple(values[r,c] for c in range(l)) for r,l in enumerate(shape))
            return
        r,c = cells[k]
        low=max(values.get((r,c-1),0), values.get((r-1,c),-1)+1)
        for a in range(low,m):
            values[r,c]=a
            yield from rec(k+1)
        values.pop((r,c),None)
    yield from rec(0)


def permutation_sign(p):
    return -1 if sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))%2 else 1


def literal_u(table):
    shape=tuple(map(len,table))
    rows=[]; start=0
    for l in shape:
        rows.append(tuple(range(start,start+l)));start+=l
    cols=[tuple(rows[r][c] for r in range(len(rows)) if len(rows[r])>c) for c in range(shape[0] if shape else 0)]
    rowwords=[sorted(set(permutations(row))) for row in table]
    colperms=[[(p,permutation_sign(p)) for p in permutations(range(len(col)))] for col in cols]
    answer=Counter()
    for roworder in product(*rowwords):
        rho=tuple(a for row in roworder for a in row)
        for permchoice in product(*colperms):
            x=list(rho);sgn=1
            for col,(p,s) in zip(cols,permchoice):
                sgn*=s
                for t,target in enumerate(col):x[target]=rho[col[p[t]]]
            answer[tuple(x)]+=sgn
    return {k:v for k,v in answer.items() if v}


def multiply(a,b):
    c=Counter()
    for ka,va in a.items():
        for kb,vb in b.items():c[tuple(x+y for x,y in zip(ka,kb))]+=va*vb
    return {k:v for k,v in c.items() if v}


def base(shape,m):
    zero=(0,)*(m*m);answer={zero:1}
    for length in range(1,m+1):
        power=(shape[length-1] if length<=len(shape) else 0)-(shape[length] if length<len(shape) else 0)
        det={}
        for p in permutations(range(length)):
            key=list(zero)
            for a,b in enumerate(p):key[a*m+b]+=1
            det[tuple(key)]=factorial(length)*permutation_sign(p)
        for _ in range(power):answer=multiply(answer,det)
    return answer


def divided_move(poly,j,i,e,m,side):
    if e==0:return poly
    out=Counter()
    sources=[j*m+b for b in range(m)] if side==0 else [a*m+j for a in range(m)]
    targets=[i*m+b for b in range(m)] if side==0 else [a*m+i for a in range(m)]
    for key,v in poly.items():
        ranges=[range(min(e,key[s])+1) for s in sources]
        for ks in product(*ranges):
            if sum(ks)!=e:continue
            moved=list(key);coef=v
            for s,t,k in zip(sources,targets,ks):
                coef*=comb(key[s],k);moved[s]-=k;moved[t]+=k
            out[tuple(moved)]+=coef
    return {k:v for k,v in out.items() if v}


def lower(poly,table,m,side):
    for j in range(len(table)-1,-1,-1):
        for i in range(m-1,j,-1):poly=divided_move(poly,j,i,table[j].count(i),m,side)
    return poly


def explicit_pair(a,b,m):
    answer=Counter()
    for x,cx in a.items():
        for y,cy in b.items():
            key=[0]*(m*m)
            for s,t in zip(x,y):key[s*m+t]+=1
            answer[tuple(key)]+=cx*cy
    return {k:v for k,v in answer.items() if v}


def main():
    started=time.perf_counter();stats=[];total=0
    for m,max_n in [(2,7),(3,5),(4,4)]:
        for n in range(max_n+1):
            cases=0;shapes=0
            for shape in partitions(n):
                if len(shape)>m:continue
                tabs=list(semistandard(shape,m));us=[literal_u(t) for t in tabs];h=base(shape,m)
                for i,t in enumerate(tabs):
                    left=lower(h,t,m,0)
                    for j in range(i,len(tabs)):
                        lhs=lower(left,tabs[j],m,1);rhs=explicit_pair(us[i],us[j],m)
                        assert lhs==rhs,(m,n,shape,t,tabs[j])
                        cases+=1
                shapes+=1
            total+=cases
            stats.append(dict(m=m,n=n,shapes=shapes,pairs=cases))
            print(json.dumps(stats[-1]),flush=True)
    large=[]
    for shape in [(3,2,1),(2,2,2),(3,1,1,1),(4,2),(4,2,1),(3,3,1),(3,2,1,1),(2,2,2,1)]:
        m=4;tabs=list(semistandard(shape,m))
        def richness(t):
            return (sum(len(set(row)) for row in t),sum(a>r for r,row in enumerate(t) for a in row),t)
        chosen=sorted(set([tabs[0],max(tabs,key=richness),tabs[-1]]))
        us=[literal_u(t) for t in chosen];h=base(shape,m);cases=0
        for i,t in enumerate(chosen):
            left=lower(h,t,m,0)
            for j in range(i,len(chosen)):
                assert lower(left,chosen[j],m,1)==explicit_pair(us[i],us[j],m),(shape,t,chosen[j])
                cases+=1
        total+=cases
        detail=dict(m=m,n=sum(shape),shape=shape,tableaux=chosen,pairs=cases,support_sizes=list(map(len,us)))
        large.append(detail);print(json.dumps(detail),flush=True)
    result=dict(status='PASS',arithmetic='Python arbitrary precision integers',cases=total,complete_small_details=stats,selected_larger_details=large,seconds=time.perf_counter()-started)
    Path(__file__).with_name('results').mkdir(exist_ok=True)
    Path(__file__).with_name('results').joinpath('theory_selftest.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
if __name__=='__main__':main()
