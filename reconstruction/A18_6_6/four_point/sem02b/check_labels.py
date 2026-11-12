#!/usr/bin/env python3
"""Independently enumerate every D0 semistandard-tableau label.

Input coefficient entries are never accessed. The partition labels specify
which shapes to audit; for each shape, all SSYT of every allowed content are
generated without using the supplied tableau list. Actual and expected lists
are compared only after independent generation.
"""
from __future__ import annotations
import argparse
import gzip
import json
from pathlib import Path
import time

T_VALUES=(0,3,4,5,6)


def semistandard_tableaux(shape,content):
    shape=tuple(shape)
    content=list(content)
    if sum(shape)!=sum(content):
        return []
    rows=[[-1]*length for length in shape]
    cells=[(r,c) for r,length in enumerate(shape) for c in range(length)]
    answer=[]
    def visit(pos):
        if pos==len(cells):
            assert not any(content)
            answer.append([row.copy() for row in rows])
            return
        r,c=cells[pos]
        lo=0
        if c:
            lo=max(lo,rows[r][c-1])
        if r:
            lo=max(lo,rows[r-1][c]+1)
        for symbol in range(lo,len(content)):
            if content[symbol]:
                content[symbol]-=1
                rows[r][c]=symbol
                visit(pos+1)
                content[symbol]+=1
        rows[r][c]=-1
    visit(0)
    # Row-major recursion in increasing symbols has already produced this order.
    assert answer==sorted(answer,key=lambda tab:sum(tab,[]))
    return answer


def partitions(total,max_rows,max_part=None):
    if total==0:
        yield ()
        return
    if max_rows==0:
        return
    max_part=total if max_part is None else min(max_part,total)
    for first in range(max_part,0,-1):
        for rest in partitions(total-first,max_rows-1,first):
            yield (first,)+rest


def generate_shape(shape):
    by_t={t:semistandard_tableaux(shape,(12-t,t,t,6-t)) for t in T_VALUES}
    return sum((by_t[t] for t in T_VALUES),[]), {str(t):len(by_t[t]) for t in T_VALUES}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--model',default='upload/A18_6_6_model.json(1).gz')
    parser.add_argument('--out',default='sem02b/results/labels.json')
    args=parser.parse_args()
    start=time.perf_counter()
    # Load labels only into the working list. No block entries are accessed.
    with gzip.open(args.model,'rt') as f:
        model=json.load(f)
    labels=[(b,x['label']) for b,x in enumerate(model['blocks'])
            if x['label']['family']=='D0']
    del model
    blocks=[]
    generated={}
    for b,label in labels:
        shape=tuple(label['partition'])
        assert sum(shape)==18 and len(shape)<=4
        assert all(x>0 for x in shape)
        assert tuple(sorted(shape,reverse=True))==shape
        own,counts=generate_shape(shape)
        generated[shape]=(own,counts)
        actual=label['tableaux']
        differences=[]
        for i in range(max(len(own),len(actual))):
            ours=own[i] if i<len(own) else None
            theirs=actual[i] if i<len(actual) else None
            if ours!=theirs:
                differences.append({'tableau_index':i,'computed':ours,'model':theirs})
        blocks.append({'block':b,'partition':shape,'tableaux':len(own),
                       'columns_including_empty':len(own)+(shape==(18,)),
                       'tableaux_by_t':counts,'all_tableaux_equal':not differences,
                       'differences':differences})
    # An optional inventory of shapes before the unordered-pair projection.
    # A raw SSYT may project to zero; this inventory does not assert which
    # shapes should survive the projection.
    possible_shapes=[]
    all_partitions=list(partitions(18,4))
    for shape in all_partitions:
        if shape not in generated:
            generated[shape]=generate_shape(shape)
        own,counts=generated[shape]
        if own:
            possible_shapes.append(shape)
    stored_shapes={tuple(label['partition']) for _,label in labels}
    unstored=[shape for shape in possible_shapes if shape not in stored_shapes]
    result={
      'family':'D0',
      'blocks_checked':len(blocks),
      'block_numbers':[b for b,_ in labels],
      'all_tableaux_equal':all(b['all_tableaux_equal'] for b in blocks),
      'full_columns_total':sum(b['columns_including_empty'] for b in blocks),
      'tableaux_total':sum(b['tableaux'] for b in blocks),
      'partitions_of_18_with_at_most_4_rows':len(all_partitions),
      'partitions_with_any_allowed_raw_tableau':len(possible_shapes),
      'partitions_with_raw_tableaux_but_no_stored_block':unstored,
      'shape_inventory_caution':'Raw tableau existence alone does not imply a nonzero unordered-pair projection. No claim is made about completeness of nonzero projected shapes.',
      'elapsed_seconds':time.perf_counter()-start,
      'block_entries_accessed':False,
      'blocks':blocks,
    }
    dest=Path(args.out)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='blocks'},ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
