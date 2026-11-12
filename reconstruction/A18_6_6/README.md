# Reconstructing the A(18,6,6) model

Two independently written programs reconstruct the integer coefficients of the
constraint matrices from their mathematical definitions and compare every key
and value with the stored model. Block numbers are zero-based.
The [label specification](label_specification.md) is translated from the
[definitive Chinese original](label_specification_zh.md).

| Reconstruction | Blocks | Result from a clean run | Wall time |
|---|---|---|---:|
| [Nine families](nine_families/SEM02a_report.md) | 269–16363 | `equal_blocks=16095`, `compared_nonzero_keys=693971` | 20.68 s |
| [Four-point matrices](four_point/SEM02b_report.md) | D0, 0–64, and D2, 65–268 | All 941 variables and all 269 blocks matched, `total_nonzero_keys=1380389` | 29.34 s |

The [public verifier](../../checks/leave_rows.py) recomputes the 613 nonzero
coefficients in the remaining 23 blocks (16364–16386) from their definitions.
Together the three checks cover all 16387 blocks and all 2074973 nonzero
coefficients in the [model](https://github.com/wille614/cw-code-upper-bounds/releases/download/v2.0.0/A18_6_6_model.json.gz).
The [proof](../../notes/A18_6_6.pdf) explains how the model and its dual certificate give the bound 183.

## Reproduce the comparisons

Download [SEM02a_programs.zip](https://github.com/wille614/cw-code-upper-bounds/releases/download/v2.0.2/SEM02a_programs.zip)
and [SEM02b_programs.zip](https://github.com/wille614/cw-code-upper-bounds/releases/download/v2.0.2/SEM02b_programs.zip). Each archive includes
its inputs and recorded outputs. Extract each into a separate directory. In
`SEM02a_audit/` and then in `SEM02b_programs/`, run

```sh
python3 run_all.py
```

Use Python 3.10 or later and a C++17-capable g++ on the executable search path.
The nine-family code uses x86 POPCNT, and the four-point code uses the compiler's
`__int128` type. No third-party Python packages are required.
The times above cover complete runs on x86-64 Windows with CPython 3.13.5 and
g++ 16.2.0. Both started with empty result directories and used the default
arguments (four workers for D0). The nine-family summary is in
`data/run_summary.json`. The four-point run prints its variable and block
counts and writes `sem02b/results/delivery_audit.json`.
