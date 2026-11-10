"""Independent exact interpretation of the first integer certificate set."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb
import json

def inverse(matrix):
    # Fraction-free Gauss-Jordan elimination on [H | I].
    n=len(matrix)
    a=[list(row)+[int(i==j) for j in range(n)] for i,row in enumerate(matrix)]
    previous=1
    for k in range(n):
        pivot=a[k][k];assert pivot>0
        old=[row[:] for row in a]
        for i in range(n):
            if i==k:continue
            for j in range(2*n):
                if j==k:continue
                q,r=divmod(pivot*old[i][j]-old[i][k]*old[k][j],previous)
                assert r==0
                a[i][j]=q
            a[i][k]=0
        previous=pivot
    inv=[[F(a[i][n+j],a[i][i]) for j in range(n)] for i in range(n)]
    assert all(sum(matrix[i][k]*inv[k][j] for k in range(n))==(i==j) for i in range(n) for j in range(n))
    return inv

@lru_cache(None)
def fingerprint(edge_tuple):
    if not edge_tuple:return ()
    weights={tuple(sorted((a,b))):w for a,b,w in edge_tuple};vertices={v for a,b,w in edge_tuple for v in (a,b)};parts=[]
    while vertices:
        component={min(vertices)};before=None
        while before!=component:
            before=set(component)
            for a,b in weights:
                if a in component or b in component:component.update((a,b))
        vertices-=component;colors=defaultdict(list)
        for v in component:colors[tuple(sorted(w for e,w in weights.items() if v in e))].append(v)
        best=None
        for orders in product(*(permutations(colors[c]) for c in sorted(colors))):
            order=[x for group in orders for x in group]
            candidate=tuple(weights.get(tuple(sorted((a,b))),0) for a,b in combinations(order,2))
            if best is None or candidate<best:best=candidate
        parts.append((len(component),best))
    return tuple(sorted(parts))

def coverage(catalogue,detail):
    counts=[]
    for t in range(6):
        found=set();partitions=0
        def visit(labels,largest):
            nonlocal partitions
            if len(labels)==2*t:
                counts=Counter(tuple(sorted(labels[i:i+2])) for i in range(0,len(labels),2))
                edges=tuple(sorted((a,b,w) for (a,b),w in counts.items()))
                found.add(fingerprint(edges));partitions+=1;return
            for value in range(largest+2):
                if len(labels)%2 and value==labels[-1]:continue
                visit(labels+[value],max(largest,value))
        visit([],-1)
        claimed=[fingerprint(tuple(tuple(e) for e in g['edges'])) for g in catalogue if g['deficit']==t]
        assert len(set(claimed))==len(claimed) and set(claimed)==found
        counts.append(len(found));detail('First set coverage:',t,partitions,len(found))
    assert counts==[1,1,3,8,23,66]

def check(root,detail):
    data=json.loads((root/'certificates/A32_16_13_first.json').read_text(encoding='utf-8'))
    graphs=data['graphs'];assert [g['id'] for g in graphs]==list(range(102))
    coverage(graphs,detail);totals=Counter()
    assert 32*55+12*11==1892>5*comb(28,2)==1890
    for graph in graphs:
        case=graph['id'];edges={(a,b):w for a,b,w in graph['edges']}
        v=1+max((max(a,b) for a,b in edges),default=-1);u=27-v
        assert sum(edges.values())==graph['deficit']<=5
        assert {x for e in edges for x in e}==set(range(v))
        H=[[32*int(i==j)-4*edges.get(tuple(sorted((i,j))),0) for j in range(v)] for i in range(v)]
        inv=inverse(H);beta=[6*sum(row) for row in inv];S=32-36*(sum(map(sum,inv))+F(u,32))
        if S<0:
            totals['GRAM']+=1;continue
        badpairs=list(combinations(range(v),2));types=[];columns=[];errors=[];diagonals=[];signs=[]
        rhs=[32,13*u,5*comb(u,2)]+[13]*v+[5*u]*v+[5-edges.get((i,j),0) for i,j in badpairs]
        for mask in range(1<<v):
            bits=[mask>>i&1 for i in range(v)];sign=[1-2*x for x in bits]
            norm=sum(sign[i]*inv[i][j]*sign[j] for i in range(v) for j in range(v))+F(u,32)
            for a in range(u+1):
                r=a+sum(bits)
                if (r-10)*(r-11)//2>5-sum(edges.values()):continue
                error=1-sum(x*y for x,y in zip(beta,sign))-F(3,16)*(u-2*a)
                if S==0:
                    if error or norm>1:continue
                    diagonal=norm
                else:
                    diagonal=norm+error*error/S
                    if diagonal>1:continue
                types.append((mask,a));signs.append(sign);errors.append(error);diagonals.append(diagonal)
                columns.append([1,a,comb(a,2)]+bits+[a*x for x in bits]+[bits[i]*bits[j] for i,j in badpairs])
        if str(case) in data['single']:
            c=data['single'][str(case)];y=c['coefficients'];assert len(y)==len(rhs)
            dots=[sum(x*z for x,z in zip(col,y)) for col in columns];right=sum(x*z for x,z in zip(rhs,y))
            assert max(dots)==c['maximum_column_dot']<=0 and right==c['positive_rhs']>0
            totals['SINGLE']+=1;continue
        c=data['pair'][str(case)];y=c['coefficients'];offset=len(rhs);local=6+3*v+len(badpairs);T=len(types)
        assert S>0 and len(y)==offset+T*local and T==c['types']
        products=[]
        for i,(mask,a) in enumerate(types):
            bits=[mask>>k&1 for k in range(v)]
            target=[31,12*a,13*(u-a),4*comb(a,2),5*a*(u-a),5*comb(u-a,2)]+[13-x for x in bits]+[a*(5-x) for x in bits]+[5*(u-a)]*v+[5-edges.get((k,l),0)-bits[k]*bits[l] for k,l in badpairs]
            products.append(sum(x*z for x,z in zip(columns[i],y[:offset]))-sum(x*z for x,z in zip(target,y[offset+i*local:offset+(i+1)*local])))
        count=0
        for i,(mask,a) in enumerate(types):
            for j in range(i,T):
                aa=types[j][1];part=sum(signs[i][r]*inv[r][s]*signs[j][s] for r in range(v) for s in range(v))
                for t in range(max(0,a+aa-u),min(a,aa)+1):
                    off=part+F(u-2*a-2*aa+4*t,32)+errors[i]*errors[j]/S
                    if off*off>(1-diagonals[i])*(1-diagonals[j]):continue
                    value=0
                    for first,second in ((i,j),(j,i)) if i!=j else ((i,j),):
                        bits=[types[second][0]>>k&1 for k in range(v)];h=types[second][1]-t
                        feature=[1,t,h,comb(t,2),t*h,comb(h,2)]+bits+[t*x for x in bits]+[h*x for x in bits]+[bits[k]*bits[l] for k,l in badpairs]
                        value+=sum(x*z for x,z in zip(feature,y[offset+first*local:offset+(first+1)*local]))
                    products.append(value);count+=1
        right=sum(x*z for x,z in zip(rhs,y[:offset]))
        assert max(products)==c['maximum_column_dot']<=0 and right==c['positive_rhs']>0 and count==c['pair_variables']
        detail('First pair certificate:',case,'types',T,'pairs',count,'max column',max(products),'positive rhs',right)
        totals['PAIR']+=1
    assert totals=={'GRAM':72,'SINGLE':22,'PAIR':8}
    detail('First certificate set:',dict(totals))
    return totals
