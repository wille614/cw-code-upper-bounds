"""Finite checks of explicit codes, constructions and arithmetic corollaries."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import comb
import json

def read_code(root, item, detail):
    n,d,w=item['n'],item['d'],item['w']
    lines=(root/item['file']).read_text(encoding='utf-8').splitlines()
    if lines[0]=='$BASE=16':
        assert all(s and set(s)<=set('0123456789ABCDEFabcdef') for s in lines[1:])
        assert all(len(s)==(n+3)//4 for s in lines[1:])
        words=[int(s,16) for s in lines[1:]]
    else:
        assert all(len(s)==n and set(s)<={'0','1'} for s in lines)
        words=[int(s,2) for s in lines]
    assert len(words)==item['count']==len(set(words))
    assert all(0<=x<1<<n and x.bit_count()==w for x in words)
    histogram=[0]*(n+1)
    for i,x in enumerate(words):
        for j in range(i):histogram[(x^words[j]).bit_count()]+=1
    observed={str(i):count for i,count in enumerate(histogram) if count}
    assert sum(histogram)==comb(len(words),2)
    assert min(map(int,observed))>=d
    assert observed==item['distance_distribution']
    detail('Words:',len(words),'length:',n,'weight:',w,'distinct: yes')
    detail('All unordered pairs:',sum(histogram),'distance distribution:',json.dumps(observed,sort_keys=True))
    return words

def translate64(words, detail):
    base=(((0,0),(1,0),(2,1),(1,2),(6,2),(4,5)),((0,0),(2,0),(0,1),(2,3),(1,5),(5,7)))
    differences=Counter(((a-c)%8,(b-d)%8) for block in base for a,b in block for c,d in block if (a,b)!=(c,d))
    group=set(product(range(8),repeat=2));hole=set(product((0,4),repeat=2))
    assert set(differences)==group-hole and set(differences.values())=={1}
    translates={sum(1<<(63-((x+a)%8+8*((y+b)%8))) for x,y in block) for block in base for a,b in group}
    assert len(translates)==128 and translates==set(words)
    bound=64*((63)//5)//6
    assert bound==128
    detail('60 ordered differences cover exactly the complement of {0,4} x {0,4}.')
    detail('128 distinct translates equal the supplied code; Johnson: floor(64/6 * floor(63/5)) =',bound)

def ternary_words(g):
    result=[]
    for a in product(range(3),repeat=5):
        word=[0]*(len(g)+4)
        for i,x in enumerate(a):
            for j,y in enumerate(g):word[i+j]=(word[i+j]+x*y)%3
        result.append(tuple(word))
    assert len(set(result))==243
    return result

def ternary(root, n, words, detail):
    data=json.loads((root/f'data/A{n}_ternary.json').read_text(encoding='utf-8'))
    outer=ternary_words(data['polynomial_ascending'])
    dist=Counter(sum(x!=0 for x in c) for c in outer)
    assert {str(k):v for k,v in dist.items()}==data['outer_weight_distribution']
    assert min(k for k in dist if k)>=(5 if n==30 else 6)
    image={sum(1<<(3*i+a) for i,a in enumerate(word)) for word in outer}
    if n==32:image={x>>1 for x in image if not x&1}
    final=set(words);removed={int(x,16) for x in data['removed_hex']};added={int(x,16) for x in data['added_hex']}
    assert removed==image-final and added==final-image
    assert final==(image-removed)|added
    if n==30:
        assert len(image)==243 and not removed and len(added)==130
        assert all(sum(bool(x&(7<<(3*i))) for i in range(10))<=5 for x in added)
        six_sets=set()
        for x in added:
            for subset in combinations([i for i in range(30) if x>>i&1],6):
                assert subset not in six_sets
                six_sets.add(subset)
        assert len(six_sets)==27300
        detail('All 130 added words occupy at most five triples; all 27300 six-subsets are distinct.')
    else:assert len(image)==162
    detail('Ternary outer code: length',len(outer[0]),'dimension 5; weight distribution',dict(sorted(dist.items())))
    detail('Unit-vector image after the stated shortening:',len(image),'removed:',len(removed),'added:',len(added),'final:',len(final))

def plane58(detail):
    points=[p for p in product(range(7),repeat=3) if any(p) and next(x for x in p if x)==1]
    lines=[frozenset(i for i,p in enumerate(points) if sum(a*b for a,b in zip(p,q))%7==0) for q in points]
    assert len(points)==len(lines)==57 and len(set(lines))==57
    assert all(len(x)==8 for x in lines) and all(len(a&b)==1 for a,b in combinations(lines,2))
    assert 57//7==8 and 58*8//8==58 and 8*8-8+2==58
    detail('PG(2,7): 57 points, 57 distinct lines, line size 8, every line pair meets once; append one zero coordinate.')
    detail('Point degree <= floor(57/7)=8; size 58 forces the excluded symmetric configuration 58_8.')
    detail('Literature premise: Davydov et al., Theorem 6.1, citing Gropp, Theorem 2.4; k=8 lies in 5..10.')

def packing(n,k,t,detail):
    assert 3<=t<=k<n
    for s in range(1,t):
        top,div=comb(n-s,t-s),comb(k-s,t-s)
        assert top%div==0
        detail('s =',s,':',top,'divisible by',div)
    total,q=comb(n,t),comb(k,t)
    assert total%q
    bound=total//q-2
    detail('Total',total,'block size',q,'remainder',total%q,'floor(total/block)-2 =',bound)
    return bound

def jpf(n,k,detail):
    assert n>k>=4 and (n-2)%(k-2)==0
    assert ((n-1)*(n-2))%((k-1)*(k-2))==0
    assert (n*(n-1)*(n-2))%(k*(k-1)*(k-2))!=0
    j=n*((n-1)*((n-2)//(k-2))//(k-1))//k
    detail('Theorem 13 hypotheses: PASS; nested Johnson =',j,'; subtract 3 =',j-3)
    return j-3

def binomial(n,k):return comb(n,k) if n>=0 and 0<=k<=n else 0
def coefficient(n,w,k,i):
    numerator=sum((-1)**(i-j)*binomial(w-j,i-j)*binomial(w-k,j)*binomial(n-w+j-k,j) for j in range(i+1))
    return Fraction(numerator,comb(w,i)*comb(n-w,i))

def linear(root, item, detail):
    c=json.loads((root/item['certificate']).read_text(encoding='utf-8'))
    n,d,w=(c[x] for x in ('n','d','w'));ii=list(range(d//2,w+1))
    assert [n,d,w]==[item[x] for x in ('n','d','w')] and c['distances']==ii
    if 'constraints' not in c:
        rows=[[-coefficient(n,w,k,i) for i in ii] for k in range(1,w+1)];rhs=[Fraction(1)]*w
    else:
        dependencies=json.loads((root/'data/shell_bounds.json').read_text(encoding='utf-8'))
        def norm(p):
            m,a,n,b,d=p;pairs=sorted(((m,min(a,m-a)),(n,min(b,n-b))))
            return [*pairs[0],*pairs[1],d]
        @lru_cache(None)
        def shell(key):
            x=dependencies['T'][key];m,a,v,b,d=x['parameters'];method=x['method'];bound=x['bound']
            assert norm(x['parameters'])==x['parameters']
            if method=='A_SINGLE':
                assert a==0 or b==0
                name=f'A{v}_{d}_{b}' if a==0 else f'A{m}_{d}_{a}'
                assert x['A_ref']==name and bound==dependencies['A'][name]['bound']
            elif method=='DISJOINT':assert 2*(a+b)==d and bound==min(t//s for t,s in ((m,a),(v,b)) if s)
            else:
                p,mul,den={'J11':([m-1,a-1,v,b,d],m,a),'J10':([m-1,a,v,b,d],m,m-a),'J21':([m,a,v-1,b-1,d],v,b),'J20':([m,a,v-1,b,d],v,v-b)}[method]
                assert norm(p)==dependencies['T'][x['parent']]['parameters']
                assert x['multiplier']==mul and x['denominator']==den and bound==mul*shell(x['parent'])//den
            detail(key,'<=',bound,method)
            return bound
        rows=[];rhs=[]
        for label,y in zip(c['constraints'],c['dual']):
            if label['type']=='Delsarte':rows.append([-coefficient(n,w,label['k'],i) for i in ii]);rhs.append(Fraction(1))
            elif label['type']=='shell':
                index=label['i'];bound=shell(label['proof'])
                assert norm([w,index,n-w,index,d])==dependencies['T'][label['proof']]['parameters']
                rows.append([Fraction(i==index) for i in ii]);rhs.append(Fraction(bound))
            else:raise ValueError('Unrecognized inequality')
    y=list(map(Fraction,c['dual']));assert len(y)==len(rows) and all(a>=0 for a in y)
    columns=[sum(y[k]*rows[k][j] for k in range(len(rows))) for j in range(len(ii))]
    assert all(x>=1 for x in columns)
    bound=1+sum(a*b for a,b in zip(rhs,y));assert bound==Fraction(c['objective'])
    assert bound.numerator//bound.denominator==item['bound']==c['integer_bound']
    detail('Nonnegative dual multipliers:',*[str(v) for v in y])
    detail('Each distance column has weighted coefficient >= 1:',*[str(v) for v in columns])
    detail('Exact rational upper bound:',str(bound),'; floor:',item['bound'])
    return bound
