"""Exact rational dual checks, using only the Python standard library."""
from fractions import Fraction as Q
from math import comb
import gzip,hashlib,json
from .positive import psd
from .leave_rows import check as check_leave

def digest(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def matrix18(root,detail):
    raw=gzip.decompress((root/'attachments/A18_6_6_model.json.gz').read_bytes())
    m=json.loads(raw)
    with gzip.open(root/'attachments/A18_6_6_certificate.json.gz','rt',encoding='utf-8') as f:c=json.load(f)
    assert c['schema']=='rational-matrix-certificate' and c['model_sha256']==hashlib.sha256(raw).hexdigest()
    assert (m['n'],m['d'],m['w'])==(18,6,6)
    assert c['objective_variable']==m['objective_variable']==1 and c['box_bound']==m['box_bound']==comb(18,6)==18564
    assert m['cardinality_bound']==186 and len(m['variables'])==941
    assert len(c['dual_blocks'])==len(m['blocks'])==16387
    check_leave(m,detail)
    coeff=[Q(0)]*(len(m['variables'])+1)
    for block,dual in zip(m['blocks'],c['dual_blocks']):
        matrix=[[Q(v) for v in row] for row in dual['matrix']]
        assert len(matrix)==block['size'];psd(matrix)
        for variable,i,j,value in block['entries']:
            assert 0<=variable<len(coeff) and 0<=i<=j<len(matrix) and type(value) is int
            coeff[variable]+=matrix[i][j]*value*(1 if i==j else 2)
    residual=[coeff[i]+int(i==1) for i in range(1,len(coeff))]
    assert residual==[Q(v) for v in c['residuals']]
    correction=18564*sum((max(v,Q(0)) for v in residual),Q(0));bound=coeff[0]+correction
    assert coeff[0]==Q(c['constant_term']) and correction==Q(c['residual_compensation']) and bound==Q(c['bound'])
    assert bound==Q(18397886022290943549115383775541,100000000000000000000000000000)
    assert bound.numerator//bound.denominator==183
    detail('All 23 leave rows agree coefficient by coefficient; all 16387 rational dual matrices are positive semidefinite.')
    detail('Constant:',str(coeff[0]),'positive-residual compensation:',str(correction))
    detail('Exact rational upper bound:',str(bound),'; floor: 183')
    return bound

def blocks(root,index):
    count=0;entries=0
    for part in index['shards']:
        path=(root/part['file']).resolve()
        assert path.is_relative_to(root.resolve()) and path.stat().st_size==part['bytes']
        assert digest(path)==part['sha256']
        assert part['first_block']==count
        before=count
        with gzip.open(path,'rt',encoding='utf-8') as f:
            for line in f:
                b=json.loads(line);assert set(b)=={'label','size','entries'}
                count+=1;entries+=len(b['entries'])
                yield b
        assert count-before==part['blocks']
    assert count==index['block_count'] and entries==index['coefficient_count']

def gram21(root,detail):
    fp=root/'data/A21_8_9_full.json';cp=root/'data/A21_8_9_principal.json'
    full=json.loads(fp.read_text(encoding='utf-8'));principal=json.loads(cp.read_text(encoding='utf-8'))
    with gzip.open(root/'attachments/A21_8_9_certificate.json.gz','rt',encoding='utf-8') as f:c=json.load(f)
    assert c['schema']=='rational-gram-certificate'
    assert c['full_model_sha256']==digest(fp) and c['model_sha256']==digest(cp)
    for key in ('n','d','w','normalization','level','objective_variable','box_bound','variables'):assert full[key]==principal[key]
    assert full['level']=='A4' and full['normalization']=='z(S)=binom(n,w)*average(1[S subset C])'
    assert (full['n'],full['d'],full['w'])==(c['n'],c['d'],c['w'])==(21,8,9)
    assert full['objective_variable']==c['objective_variable']==1
    box=comb(21,9);assert box==full['box_bound']==c['box_bound']==293930
    assert len(full['variables'])==10665 and full['block_count']==11372 and principal['block_count']==11368
    assert full['coefficient_count']==558841903 and principal['coefficient_count']==156029979
    den=int(c['factor_denominator']);assert den>0
    for factor in c['factors']:assert all(type(p) is int and p>=0 for p in factor['powers'])
    maximum=max(p for factor in c['factors'] for p in factor['powers']);common=den**2*2**(2*maximum)
    residual=[0]*(len(full['variables'])+1);residual[1]=common
    compressed=blocks(root,principal);b=next(compressed,None);index=0;nonnegative=set();checked=0
    for fi,original in enumerate(blocks(root,full)):
        if b is None:continue
        label=b['label']
        if label['family']=='nonnegative':
            if original['label']['family']!='nonnegative':continue
            assert b==original
            v=label['variable'];assert v not in nonnegative and b['size']==1 and b['entries']==[[v,0,0,1]]
            nonnegative.add(v)
        else:
            if label['source_block']!=fi:
                assert label['source_block']>fi;continue
            assert label['family']==original['label']['family'] in ('D0','D2')
            kept=label['kept_columns'];assert kept==sorted(set(kept)) and all(0<=x<original['size'] for x in kept)
            ix={old:i for i,old in enumerate(kept)};assert b['size']==len(kept)>0
            actual=iter(b['entries'])
            for v,i,j,a in original['entries']:
                if i in ix and j in ix:assert [v,ix[i],ix[j],a]==next(actual,None)
            assert next(actual,None) is None
        factor=c['factors'][index];n=b['size'];rr=factor['R'];pp=factor['powers']
        assert len(rr)==len(pp)==n and all(len(row)==n for row in rr)
        assert all(type(x) is int for row in rr for x in row)
        matrix=[[sum(rr[i][k]*rr[j][k] for k in range(n))*2**(2*maximum-pp[i]-pp[j]) for j in range(n)] for i in range(n)]
        for v,i,j,a in b['entries']:
            assert 0<=v<len(residual) and 0<=i<=j<n and type(a) is int
            residual[v]+=a*matrix[i][j]*(1 if i==j else 2);checked+=1
        index+=1;b=next(compressed,None)
    assert b is None and index==len(c['factors'])==11368 and checked==156029979
    assert nonnegative==set(range(1,len(residual)))
    residual=[Q(x,common) for x in residual]
    assert residual==list(map(Q,c['residuals'])) and all(x<=0 for x in residual[1:])
    correction=box*sum((max(x,0) for x in residual[1:]),Q(0));bound=residual[0]+correction
    assert correction==Q(c['residual_compensation'])==0 and bound==Q(c['bound'])
    assert bound==Q(57529569329575917506024543440753,200000000000000000000000000000)
    assert bound.numerator//bound.denominator==c['integer_bound']==287
    detail('Full model coefficients:',full['coefficient_count'],'; all retained principal coefficients:',checked)
    detail('11368 rational Gram factors: PASS; 10665 nonpositive residuals: PASS; compensation: 0')
    detail('Exact rational upper bound:',str(bound),'; floor: 287')
    return bound
