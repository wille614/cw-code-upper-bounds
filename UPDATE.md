# Proposed constant-weight-code table updates

Table values checked on 10 November 2026. In the copyable list, “upper” and “lower” specify the changed endpoint; “exact” replaces the whole interval. The four tables below give the complete current numerical cell.

```text
A(64,10,6) (exact): 126-128 -> 128
A(58,14,8) (exact): 57-58 -> 57
A(21,8,9) (upper): 294 -> 287
A(18,6,6) (upper): 186 -> 183
A(32,16,13) (upper): 27 -> 26
A(27,8,7) (upper): 500 -> 499
A(28,8,8) (upper): 1750 -> 1746
A(31,10,7) (upper): 109 -> 106
A(29,6,5) (upper): 364 -> 362
A(29,6,8) (upper): not listed -> 16965
A(29,6,9) (upper): not listed -> 39585
A(29,6,10) (upper): not listed -> 90465
A(29,6,11) (upper): not listed -> 162226
A(29,6,12) (upper): not listed -> 228441
A(29,6,13) (upper): not listed -> 289604
A(29,6,14) (upper): not listed -> 322266
A(24,6,12) (lower): 5617 -> 5625
A(25,6,11) (lower): 8720 -> 8733
A(25,6,12) (lower): 9895 -> 9910
A(26,6,11) (lower): 12222 -> 12225
A(26,6,13) (lower): 16225 -> 16254
A(28,8,11) (lower): 3791 -> 3832
A(28,8,14) (lower): 6359 -> 6392
A(30,10,10) (lower): 322 -> 373
A(32,12,11) (lower): 187 -> 215
A(47,8,7) (upper): 5095 -> 5094
A(33,10,7) (upper): 150 -> 136
A(34,10,8) (upper): 637 -> 578
A(38,6,5) (upper): 842 -> 840
A(53,6,5) (upper): 2341 -> 2339
A(52,10,7) (upper): 630 -> 628
A(62,10,7) (upper): 1079 -> 1077
A(30,6,6) (upper): 1820 -> 1810
A(39,6,6) (upper): 5473 -> 5460
A(54,6,6) (upper): 21069 -> 21051
A(53,10,8) (upper): 4173 -> 4160
A(63,10,8) (upper): 8497 -> 8481
A(34,10,7) (upper): 160 -> 155
A(37,10,7) (upper): 222 -> 195
A(47,10,7) (upper): 462 -> 449
A(33,10,8) (upper): 573 -> 511
A(35,10,8) (upper): 700 -> 678
A(38,10,8) (upper): 1053 -> 926
A(48,10,8) (upper): 2772 -> 2694
A(42,12,8) (upper): 157 -> 152
A(43,12,8) (upper): 192 -> 182
A(44,12,8) (upper): 234 -> 209
```

## Verification resources

Reference measurements: CPython 3.12.3, AMD EPYC 7282 16-Core Processor; one verification process, no third-party Python packages. Wall time includes interpreter startup; memory is maximum resident set size. These are measured reference values, not hardware-independent promises.

| Command | Wall seconds | Peak resident bytes |
|---|---|---|
| `python verify.py` | 38.18711095198523 | 22855680 |
| `python verify.py --all` | 633.0763699070085 | 13355847680 |

The per-entry seconds below are the measured checks inside the complete invocation. Memory figures are conservative references: the ordinary-command peak for its entries, and the complete-command peak for either large certificate. They are whole-process peaks, not separately measured per-entry allocations. Written theorems are supplied in the linked notes; the program checks their stated finite parts.

## Exact values

| Parameters | Current table cell | New endpoint/value | Method | Evidence | Command | Measured seconds; memory reference |
|---|---|---|---|---|---|---|
| A(64,10,6) | 126-128 | = 128 | Two relative difference blocks and Johnson counting | [notes/A64_10_6.md](notes/A64_10_6.md); [codes/a64.10.6.128](codes/a64.10.6.128) | `python verify.py 64 10 6` | 0.002 s; <= 22855680 bytes |
| A(58,14,8) | 57-58 | = 57 | Symmetric-configuration nonexistence and PG(2,7) | [notes/A58_14_8.md](notes/A58_14_8.md) | `python verify.py 58 14 8` | 0.003 s; <= 22855680 bytes |

## Requested entry

| Parameters | Current table cell | New endpoint/value | Method | Evidence | Command | Measured seconds; memory reference |
|---|---|---|---|---|---|---|
| A(21,8,9) | 280-294 | <= 287 | Full four-point bound A4; rational Gram dual | [notes/A21_8_9.md](notes/A21_8_9.md); [data/A21_8_9_full.json](data/A21_8_9_full.json); [data/A21_8_9_principal.json](data/A21_8_9_principal.json); [attachments/A21_8_9_certificate.json.gz](https://github.com/wille614/cw-code-upper-bounds/releases/download/v2.0.0/A21_8_9_certificate.json.gz) | `python verify.py 21 8 9` | 573.346 s; <= 13355847680 bytes |

## Main tables

| Parameters | Current table cell | New endpoint/value | Method | Evidence | Command | Measured seconds; memory reference |
|---|---|---|---|---|---|---|
| A(18,6,6) | 133-186 | <= 183 | A4 with local packing and valid counting inequalities | [notes/A18_6_6.md](notes/A18_6_6.md); [attachments/A18_6_6_certificate.json.gz](https://github.com/wille614/cw-code-upper-bounds/releases/download/v2.0.0/A18_6_6_certificate.json.gz); [attachments/A18_6_6_model.json.gz](https://github.com/wille614/cw-code-upper-bounds/releases/download/v2.0.0/A18_6_6_model.json.gz) | `python verify.py 18 6 6` | 21.169 s; <= 13355847680 bytes |
| A(32,16,13) | 25-27 | <= 26 | 102 defect graphs; two sets of integer separation certificates | [notes/A32_16_13.md](notes/A32_16_13.md); [certificates/A32_16_13_first.json](certificates/A32_16_13_first.json); [certificates/A32_16_13_second.json](certificates/A32_16_13_second.json) | `python verify.py 32 16 13` | 4.992 s; <= 22855680 bytes |
| A(27,8,7) | 303-500 | <= 499 | Divisible-leave lemma with t=4 | [notes/packing.md](notes/packing.md) | `python verify.py 27 8 7` | 0.000 s; <= 22855680 bytes |
| A(28,8,8) | 1057-1750 | <= 1746 | Johnson from 499 | [notes/packing.md](notes/packing.md) | `python verify.py 28 8 8` | 0.000 s; <= 22855680 bytes |
| A(31,10,7) | 62-109 | <= 106 | Delsarte inequalities plus a distance-shell bound | [notes/linear_programming.md](notes/linear_programming.md); [certificates/A31_10_7_lp.json](certificates/A31_10_7_lp.json); [data/shell_bounds.json](data/shell_bounds.json) | `python verify.py 31 10 7` | 0.001 s; <= 22855680 bytes |
| A(29,6,5) | 333-364 | <= 362 | Javadi–Poorhadi–Fallah, preprint Theorem 13 | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 29 6 5` | 0.000 s; <= 22855680 bytes |
| A(29,6,8) | 8148-... | <= 16965 | Johnson from the preceding length | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 29 6 8` | 0.000 s; <= 22855680 bytes |
| A(29,6,9) | 14781-... | <= 39585 | Johnson from the preceding length | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 29 6 9` | 0.000 s; <= 22855680 bytes |
| A(29,6,10) | 27183-... | <= 90465 | Johnson from the preceding length | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 29 6 10` | 0.000 s; <= 22855680 bytes |
| A(29,6,11) | 42042-... | <= 162226 | Delsarte inequalities plus a distance-shell bound | [notes/linear_programming.md](notes/linear_programming.md); [certificates/A29_6_11_lp.json](certificates/A29_6_11_lp.json); [data/shell_bounds.json](data/shell_bounds.json) | `python verify.py 29 6 11` | 0.018 s; <= 22855680 bytes |
| A(29,6,12) | 61823-... | <= 228441 | Delsarte inequalities | [notes/linear_programming.md](notes/linear_programming.md); [certificates/A29_6_12_lp.json](certificates/A29_6_12_lp.json); [data/shell_bounds.json](data/shell_bounds.json) | `python verify.py 29 6 12` | 0.010 s; <= 22855680 bytes |
| A(29,6,13) | 80837-... | <= 289604 | Delsarte inequalities | [notes/linear_programming.md](notes/linear_programming.md); [certificates/A29_6_13_lp.json](certificates/A29_6_13_lp.json); [data/shell_bounds.json](data/shell_bounds.json) | `python verify.py 29 6 13` | 0.008 s; <= 22855680 bytes |
| A(29,6,14) | 92366-... | <= 322266 | Delsarte inequalities | [notes/linear_programming.md](notes/linear_programming.md); [certificates/A29_6_14_lp.json](certificates/A29_6_14_lp.json); [data/shell_bounds.json](data/shell_bounds.json) | `python verify.py 29 6 14` | 0.012 s; <= 22855680 bytes |
| A(24,6,12) | 5617-15906 | >= 5625 | Explicit binary code; all pairwise distances checked | [notes/lower_bounds.md](notes/lower_bounds.md); [codes/a24.6.12.5625H](codes/a24.6.12.5625H) | `python verify.py 24 6 12` | 1.573 s; <= 22855680 bytes |
| A(25,6,11) | 8720-24630 | >= 8733 | Explicit binary code; all pairwise distances checked | [notes/lower_bounds.md](notes/lower_bounds.md); [codes/a25.6.11.8733H](codes/a25.6.11.8733H) | `python verify.py 25 6 11` | 3.764 s; <= 22855680 bytes |
| A(25,6,12) | 9895-30587 | >= 9910 | Explicit binary code; all pairwise distances checked | [notes/lower_bounds.md](notes/lower_bounds.md); [codes/a25.6.12.9910H](codes/a25.6.12.9910H) | `python verify.py 25 6 12` | 4.852 s; <= 22855680 bytes |
| A(26,6,11) | 12222-42017 | >= 12225 | Explicit binary code; all pairwise distances checked | [notes/lower_bounds.md](notes/lower_bounds.md); [codes/a26.6.11.12225H](codes/a26.6.11.12225H) | `python verify.py 26 6 11` | 7.329 s; <= 22855680 bytes |
| A(26,6,13) | 16225-61174 | >= 16254 | Explicit binary code; all pairwise distances checked | [notes/lower_bounds.md](notes/lower_bounds.md); [codes/a26.6.13.16254H](codes/a26.6.13.16254H) | `python verify.py 26 6 13` | 12.996 s; <= 22855680 bytes |
| A(28,8,11) | 3791-11896 | >= 3832 | Explicit binary code; all pairwise distances checked | [notes/lower_bounds.md](notes/lower_bounds.md); [codes/a28.8.11.3832H](codes/a28.8.11.3832H) | `python verify.py 28 8 11` | 0.744 s; <= 22855680 bytes |
| A(28,8,14) | 6359-22710 | >= 6392 | Explicit binary code; all pairwise distances checked | [notes/lower_bounds.md](notes/lower_bounds.md); [codes/a28.8.14.6392H](codes/a28.8.14.6392H) | `python verify.py 28 8 14` | 2.033 s; <= 22855680 bytes |
| A(30,10,10) | 322-1569 | >= 373 | Once-shortened ternary Golay code, unit-vector expansion and 130 added words | [notes/lower_bounds.md](notes/lower_bounds.md); [codes/a30.10.10.373](codes/a30.10.10.373); [data/A30_ternary.json](data/A30_ternary.json) | `python verify.py 30 10 10` | 0.039 s; <= 22855680 bytes |
| A(32,12,11) | 187-573 | >= 215 | Dual ternary Golay code, unit-vector expansion, binary shortening and 53 added words | [notes/lower_bounds.md](notes/lower_bounds.md); [codes/a32.12.11.215](codes/a32.12.11.215); [data/A32_ternary.json](data/A32_ternary.json) | `python verify.py 32 12 11` | 0.009 s; <= 22855680 bytes |

## Other tables

| Parameters | Current table cell | New endpoint/value | Method | Evidence | Command | Measured seconds; memory reference |
|---|---|---|---|---|---|---|
| A(47,8,7) | 2517-5095 | <= 5094 | Divisible-leave lemma with t=4 | [notes/packing.md](notes/packing.md) | `python verify.py 47 8 7` | 0.000 s; <= 22855680 bytes |
| A(33,10,7) | 72-150 | <= 136 | Delsarte inequalities plus a distance-shell bound | [notes/linear_programming.md](notes/linear_programming.md); [certificates/A33_10_7_lp.json](certificates/A33_10_7_lp.json); [data/shell_bounds.json](data/shell_bounds.json) | `python verify.py 33 10 7` | 0.001 s; <= 22855680 bytes |
| A(34,10,8) | 181-637 | <= 578 | Johnson from 136 | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 34 10 8` | 0.001 s; <= 22855680 bytes |
| A(38,6,5) | 714-842 | <= 840 | Javadi–Poorhadi–Fallah, preprint Theorem 13 | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 38 6 5` | 0.000 s; <= 22855680 bytes |
| A(53,6,5) | 2145-2341 | <= 2339 | Javadi–Poorhadi–Fallah, preprint Theorem 13 | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 53 6 5` | 0.000 s; <= 22855680 bytes |
| A(52,10,7) | 421-630 | <= 628 | Javadi–Poorhadi–Fallah, preprint Theorem 13 | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 52 10 7` | 0.000 s; <= 22855680 bytes |
| A(62,10,7) | 755-1079 | <= 1077 | Javadi–Poorhadi–Fallah, preprint Theorem 13 | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 62 10 7` | 0.000 s; <= 22855680 bytes |
| A(30,6,6) | 1316-1820 | <= 1810 | Johnson applied to the preceding theorem | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 30 6 6` | 0.000 s; <= 22855680 bytes |
| A(39,6,6) | 3822-5473 | <= 5460 | Johnson applied to the preceding theorem | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 39 6 6` | 0.000 s; <= 22855680 bytes |
| A(54,6,6) | 14469-21069 | <= 21051 | Johnson applied to the preceding theorem | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 54 6 6` | 0.000 s; <= 22855680 bytes |
| A(53,10,8) | 1823-4173 | <= 4160 | Johnson applied to the preceding theorem | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 53 10 8` | 0.000 s; <= 22855680 bytes |
| A(63,10,8) | 4228-8497 | <= 8481 | Johnson applied to the preceding theorem | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 63 10 8` | 0.000 s; <= 22855680 bytes |
| A(34,10,7) | 85-160 | <= 155 | Johnson from a tabulated or preceding upper bound | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 34 10 7` | 0.000 s; <= 22855680 bytes |
| A(37,10,7) | 134-222 | <= 195 | Johnson from a tabulated or preceding upper bound | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 37 10 7` | 0.000 s; <= 22855680 bytes |
| A(47,10,7) | 333-462 | <= 449 | Johnson from a tabulated or preceding upper bound | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 47 10 7` | 0.000 s; <= 22855680 bytes |
| A(33,10,8) | 165-573 | <= 511 | Johnson from a tabulated or preceding upper bound | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 33 10 8` | 0.000 s; <= 22855680 bytes |
| A(35,10,8) | 210-700 | <= 678 | Johnson from a tabulated or preceding upper bound | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 35 10 8` | 0.000 s; <= 22855680 bytes |
| A(38,10,8) | 285-1053 | <= 926 | Johnson from a tabulated or preceding upper bound | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 38 10 8` | 0.000 s; <= 22855680 bytes |
| A(48,10,8) | 905-2772 | <= 2694 | Johnson from a tabulated or preceding upper bound | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 48 10 8` | 0.000 s; <= 22855680 bytes |
| A(42,12,8) | 105-157 | <= 152 | Johnson from a tabulated or preceding upper bound | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 42 12 8` | 0.000 s; <= 22855680 bytes |
| A(43,12,8) | 111-192 | <= 182 | Johnson from a tabulated or preceding upper bound | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 43 12 8` | 0.000 s; <= 22855680 bytes |
| A(44,12,8) | 124-234 | <= 209 | Johnson from a tabulated or preceding upper bound | [notes/corollaries.md](notes/corollaries.md) | `python verify.py 44 12 8` | 0.000 s; <= 22855680 bytes |

For the complete model of A(21,8,9), the two evidence indexes enumerate every compressed attachment. Their checksums and those of the other large files are also in [attachments/SHA256SUMS](attachments/SHA256SUMS). A(21,8,10) is outside this update.
