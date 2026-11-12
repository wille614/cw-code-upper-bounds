#!/usr/bin/env python3
"""Reproduce all nine families offline, in the user's requested order.

Requires Python >=3.10 and g++ supporting C++17 and -mpopcnt (x86 POPCNT).
Run from any working directory: python /path/to/run_all.py
"""
import json
import subprocess
import sys
import time
from common import ROOT,DATA,generate_variables
from prepare_inputs import prepare
from compare import verify_variables,compare_family
from early_families import FUNCTIONS

ORDER = ['nonnegative','cardinality-cap','cardinality-scalar',
         'matching-occupancy','anticode-vacancy-scalar','coordinate-spectrum',
         'occupancy-matrix','cardinality-rooted','anticode-vacancy-matrix']
SCRIPTS = {'coordinate-spectrum':'coordinate_spectrum.py',
           'occupancy-matrix':'occupancy_matrix.py',
           'cardinality-rooted':'cardinality_rooted.py',
           'anticode-vacancy-matrix':'anticode_matrix.py'}


def main():
    start=time.perf_counter()
    model=prepare()
    records,stats=generate_variables()
    (DATA/'orbits.json').write_text(json.dumps(records,separators=(',',':')))
    (DATA/'orbit_generation.json').write_text(json.dumps(stats,indent=2))
    variables=verify_variables(model)
    assert variables['all_equal'] and all(variables['metadata_checks'].values())
    assert stats['layer_counts']==[1,4,35,901]
    summary={'variables':variables,'families':[]}
    for family in ORDER:
        print('RECONSTRUCT',family,flush=True)
        if family in FUNCTIONS:
            path=FUNCTIONS[family]()
        else:
            if family=='cardinality-rooted':
                subprocess.run([sys.executable,str(ROOT/SCRIPTS[family]),'--self-check-only'],check=True)
            subprocess.run([sys.executable,str(ROOT/SCRIPTS[family])],check=True)
            path=DATA/('rebuilt_'+family+'.json')
        # Reconstruction file has been written completely before entries are read.
        result=compare_family(model,path)
        summary['families'].append(result)
        print(json.dumps(result,ensure_ascii=False),flush=True)
        (DATA/'run_summary.json').write_text(json.dumps(summary,indent=2))
        if result['mismatch'] or result['not_reconstructed'] or result['extra_ids']:
            raise RuntimeError('Comparison failed. Inspect differences; do not fit the model.')
    summary['elapsed_seconds']=time.perf_counter()-start
    summary['equal_blocks']=sum(f['equal'] for f in summary['families'])
    summary['compared_nonzero_keys']=sum(f['model_entries_compared'] for f in summary['families'])
    assert summary['equal_blocks']==16095
    (DATA/'run_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps({k:v for k,v in summary.items() if k not in ('variables','families')}),flush=True)


if __name__=='__main__':
    main()
