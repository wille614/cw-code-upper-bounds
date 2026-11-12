"""Strict comparisons after coefficient reconstructions have been frozen.

The builders never import this file. Model entries are read only here.
"""
from pathlib import Path
from collections import Counter
import gzip, hashlib, json, time, argparse
from common import DATA, ROOT

MODEL = ROOT/'upload/A18_6_6_model.json.gz'
EXPECTED_GZ = '157bd486c563c49bdb65cec5b00e245e8f816b10ab67d721841c34b3c5aabe77'
EXPECTED_JSON = 'ebf38af24dddea9230d2194faf394dc1666b75be5c40a574c6827ac2c8a0efc6'


def load_model():
    compressed = MODEL.read_bytes()
    assert hashlib.sha256(compressed).hexdigest() == EXPECTED_GZ
    raw = gzip.decompress(compressed)
    assert hashlib.sha256(raw).hexdigest() == EXPECTED_JSON
    return json.loads(raw)


def verify_variables(model):
    start = time.perf_counter()
    records = json.loads((DATA/'orbits.json').read_text())
    independent = [r['h'] for r in records]
    expected = model['variables']
    differences = [{'v':i+1,'model':a,'computed':b}
                   for i,(a,b) in enumerate(zip(expected,independent)) if a!=b]
    result = {'computed_count':len(independent),'model_count':len(expected),
              'layer_counts':dict(Counter(r['k'] for r in records)),
              'all_equal':independent==expected,'differences':differences,
              'metadata_checks':{k:model[k]==v for k,v in
                  [('n',18),('d',6),('w',6),('objective_variable',1),
                   ('box_bound',18564),('cardinality_bound',186)]},
              'seconds':time.perf_counter()-start}
    (DATA/'variable_comparison.json').write_text(json.dumps(result,indent=2))
    return result


def entry_map(entries, size):
    result,issues = {},[]
    for n,e in enumerate(entries):
        if not isinstance(e,list) or len(e)!=4 or any(type(x) is not int for x in e):
            issues.append({'position':n,'issue':'not four exact integers','entry':e})
            continue
        v,i,j,c = e
        key = (v,i,j)
        if not (0<=v<=941 and 0<=i<=j<size):
            issues.append({'position':n,'issue':'out of range','entry':e})
        if c==0:
            issues.append({'position':n,'issue':'explicit zero','entry':e})
        if key in result:
            issues.append({'position':n,'issue':'duplicate key','entry':e})
        result[key]=c
    return result,issues


def compare_family(model,path):
    start = time.perf_counter()
    frozen = path.read_bytes()
    obj = json.loads(frozen)
    family=obj['family']
    expected_ids=[i for i,b in enumerate(model['blocks']) if b['label']['family']==family]
    generated=obj['blocks']
    ids=[b['block_id'] for b in generated]
    assert len(ids)==len(set(ids))
    results=[]
    for b in generated:
        bid=b['block_id']
        actual=model['blocks'][bid]
        am,ai=entry_map(actual['entries'],actual['size'])
        bm,bi=entry_map(b['entries'],b['size'])
        differences=[{'key':list(k),'model':am.get(k,0),'computed':bm.get(k,0),
                      'model_key_present':k in am,'computed_key_present':k in bm}
                     for k in sorted(am.keys()|bm.keys()) if am.get(k,0)!=bm.get(k,0)]
        issues=[]
        if actual['label']['family']!=family:
            issues.append('wrong family')
        if actual['size']!=b['size']:
            issues.append({'model_size':actual['size'],'computed_size':b['size']})
        issues.extend(b.get('metadata_issues',[]))
        ok=not differences and not ai and not bi and not issues and am.keys()==bm.keys()
        canonical=json.dumps(sorted(b['entries']),separators=(',',':')).encode()
        results.append({'block_id':bid,'equal':ok,'size':b['size'],
                        'model_entries':len(am),'computed_entries':len(bm),
                        'key_sets_equal':am.keys()==bm.keys(),
                        'differences':differences,'model_structure_issues':ai,
                        'computed_structure_issues':bi,'metadata_issues':issues,
                        'computed_coefficients_sha256':hashlib.sha256(canonical).hexdigest()})
    result={'family':family,'expected_blocks':len(expected_ids),
            'reconstructed':len(ids),'equal':sum(r['equal'] for r in results),
            'mismatch':sum(not r['equal'] for r in results),
            'not_reconstructed':len(set(expected_ids)-set(ids)),
            'missing_ids':sorted(set(expected_ids)-set(ids)),
            'extra_ids':sorted(set(ids)-set(expected_ids)),
            'model_entries_compared':sum(r['model_entries'] for r in results),
            'reconstruction_seconds':obj['seconds'],
            'comparison_seconds':time.perf_counter()-start,
            'frozen_reconstruction_sha256':hashlib.sha256(frozen).hexdigest(),
            'blocks':results}
    (DATA/('comparison_'+family+'.json')).write_text(json.dumps(result,separators=(',',':')))
    return {k:v for k,v in result.items() if k!='blocks'}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('families',nargs='*')
    parser.add_argument('--variables',action='store_true')
    args=parser.parse_args()
    t=time.perf_counter(); model=load_model()
    print(json.dumps({'model_read_seconds':time.perf_counter()-t}),flush=True)
    if args.variables:
        print(json.dumps(verify_variables(model)),flush=True)
    for family in args.families:
        print(json.dumps(compare_family(model,DATA/('rebuilt_'+family+'.json'))),flush=True)
