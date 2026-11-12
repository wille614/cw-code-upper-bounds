#!/usr/bin/env python3
"""Verify the supplied model and export coefficient-free labels.

The JSON parser decodes the container, but no numerical coefficient is accessed
or copied into the reconstruction inputs. Only compare.py uses entries.
"""
import hashlib
import json
from common import ROOT, DATA
from compare import load_model


def prepare():
    DATA.mkdir(exist_ok=True)
    model = load_model()
    labels = [dict(block_id=i, **{k:v for k,v in block.items() if k!='entries'})
              for i,block in enumerate(model['blocks'])]
    (DATA/'labels.json').write_text(json.dumps(labels,ensure_ascii=False,separators=(',',':')))
    hashes = {p.name:hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (ROOT/'upload').iterdir() if p.is_file()}
    (DATA/'input_hashes.json').write_text(json.dumps(hashes,indent=2))
    return model


if __name__=='__main__':
    model = prepare()
    print(json.dumps({'blocks':len(model['blocks']),'coefficient_free_labels_written':True}))
