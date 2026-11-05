"""Independent exact replay: standard-library rationals and the SDPA text input.

Does not import the generator, compression algorithm or certificate producer.
"""
from fractions import Fraction as Q
from pathlib import Path
import gzip
import hashlib
import json
import sys
import time

def require(obj,key):
    assert isinstance(obj,dict) and key in obj,key
    return obj[key]

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return h.hexdigest()

def verify(folder):
    start=time.perf_counter()
    with gzip.open(folder/'certificate.json.gz','rt',encoding='utf-8') as f:cert=json.load(f)
    assert require(cert,'schema')=='rational-ldlt-v1'
    dat=folder/'compressed.dat-s'
    assert digest(dat)==require(cert,'sdpa_sha256')
    assert digest(folder/'compressed_model.json')==require(cert,'model_sha256')
    dual=require(cert,'dual_blocks');matrices=[]
    for number,item in enumerate(dual):
        A=[[Q(v) for v in row] for row in require(item,'matrix')]
        L=[[Q(v) for v in row] for row in require(item,'L')]
        D=[Q(v) for v in require(item,'D')]
        n=len(A)
        assert len(L)==len(D)==n and all(len(row)==n for row in A)
        assert all(len(L[i])==i+1 and L[i][i]==1 for i in range(n))
        assert all(d>=0 for d in D)
        for i in range(n):
            for j in range(i+1):
                value=sum((L[i][k]*D[k]*L[j][k] for k in range(j+1)),Q(0))
                assert value==A[i][j]==A[j][i],('EXACT_KERNEL_DISAGREEMENT',number,i,j)
        matrices.append(A)
    with dat.open('r',encoding='utf-8') as f:
        line=f.readline()
        while line.startswith('*') or line.startswith('"'):line=f.readline()
        nvars=int(line);nblocks=int(f.readline());sizes=list(map(int,f.readline().split()));objective=list(map(Q,f.readline().split()))
        assert nblocks==len(sizes)==len(matrices)
        assert sizes==[len(a) for a in matrices] and len(objective)==nvars
        obj=require(cert,'objective_variable')
        assert objective==[Q(-int(i==obj-1)) for i in range(nvars)]
        residual=[Q(0)]*(nvars+1);residual[obj]=1
        seen=set();last_block=0;coefficient_count=0
        for line in f:
            if not line.strip():continue
            fields=line.split();assert len(fields)==5
            v,b,i,j=map(int,fields[:4]);c=Q(fields[4])
            assert 0<=v<=nvars and 1<=b<=nblocks and 1<=i<=j<=sizes[b-1]
            assert b>=last_block,'SDPA blocks must be in the declared monotone order'
            if b!=last_block:seen.clear();last_block=b
            key=(v,b,i,j);assert key not in seen;seen.add(key)
            coefficient_count+=1
            contribution=c*matrices[b-1][i-1][j-1]*(1 if i==j else 2)
            residual[v]+=(-contribution if v==0 else contribution)
    upper=require(cert,'box_bound')
    # Bound is the size of the ambient constant-weight space, checked independently.
    meta=json.loads((folder/'compressed_model.json').read_text(encoding='utf-8'))
    from math import comb
    assert upper==comb(require(meta,'n'),require(meta,'w'))
    assert nvars==len(require(meta,'variables'))
    correction=upper*sum((max(x,0) for x in residual[1:]),Q(0))
    bound=residual[0]+correction
    assert list(map(str,residual[1:]))==require(cert,'residuals')
    assert bound==Q(require(cert,'bound')) and correction==Q(require(cert,'residual_compensation'))
    assert residual[0]==Q(require(cert,'constant_term'))
    result={'status':'PASS','arithmetic':'Python fractions.Fraction; no floating-point acceptance',
            'bound':str(bound),'integer_bound':bound.numerator//bound.denominator,
            'blocks':nblocks,'coefficients_checked':coefficient_count,'wall_seconds':format(time.perf_counter()-start,'.6f'),
            'certificate_sha256':digest(folder/'certificate.json.gz')}
    (folder/'verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result),flush=True)
    return result

if __name__=='__main__':
    assert len(sys.argv)==2
    verify(Path(sys.argv[1]))
