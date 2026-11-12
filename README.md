# Constant-weight codes: second table update

This package contains 47 proposed updates to
[Brouwer's table](https://aeb.win.tue.nl/codes/Andw.html): two exact values,
explicit lower-bound codes, exact upper-bound certificates, and stated
consequences of published bounds. Start with [UPDATE.md](UPDATE.md),
whose opening list is ready to copy into a table update. The four groups
are exact values; A(21,8,9); Two-dimensional tables "Bounds on A(n,d,w)";
and Lists "Values of A(n,d,w)" for n >= 29. The self-contained
mathematical notes are in [notes/](notes/). The packing proof is available in
[Markdown](notes/packing.md) and [plain text](notes/packing.txt).

## Verification

Use CPython 3.11 or later, with assertions enabled. Only the Python
standard library is needed. From the directory containing this file:

```sh
python verify.py
python verify.py --all
python verify.py 64 10 6
python verify.py --list
```

The first command checks 45 entries and explicitly skips the two large
semidefinite certificates, for A(21,8,9) and A(18,6,6). The second checks
all 47 entries and prints reference time and memory before starting.
The third checks a selected entry in detail; substitute any parameters
listed by the fourth command. Lower-bound details include the complete
unordered distance distribution. Certificate details include exact
rational bounds and their floors. The 26 bound has two complete sets of
integer certificates, each with its own independent finite checker.

The full output of the first two commands appears in
[EXPECTED_OUTPUT.txt](EXPECTED_OUTPUT.txt). Elapsed times vary; mathematical
results, counts, entry order and status lines should agree. Neither a
solver nor a network connection is used by the verifier. The referenced
nonexistence theorems and existing table upper bounds are identified
as premises in the relevant notes; their finite arithmetic is checked.
Written proofs and exact model certificates have distinct roles and are
not described as proof-assistant formalizations.

## Large attachments

An assembled package contains `attachments/`. In a repository-only copy,
the compressed files are separate
[release attachments](https://github.com/wille614/cw-code-upper-bounds/releases/tag/v2.0.0); their exact names
and hashes are listed in [attachments/SHA256SUMS](attachments/SHA256SUMS).
Place them unchanged in `attachments/`. No decompression step is needed.
The two model indexes in `data/` identify every compressed coefficient
shard. Gzip headers, index names and certificate metadata use the public
filenames. The ordinary check does not need those large attachments.

All code files use the table's filename conventions. A suffix H denotes
hexadecimal words and a first line `$BASE=16`; other code files contain
binary words. The supporting data explain the two ternary outer codes,
their unit-vector maps and their supplementary words. There are no
optimization or code-construction utilities in this package; the finite
reconstructions inside the verifier check the stated certificates.

## Earlier results and attribution

The earlier three upper bounds and their paper remain available in
[version 1.2.0](https://github.com/wille614/cw-code-upper-bounds/releases/tag/v1.2.0),
with their data in
[version 1.0.0](https://github.com/wille614/cw-code-upper-bounds/releases/tag/v1.0.0).
The present update is supplied through code files, certificates and
short proof notes.

Copyright © 2026 Jiang Chang. Code is licensed under
[Apache License 2.0](LICENSE); certificates and mathematical data under
[CC BY 4.0](LICENSE-DATA). Sources and prior-code attributions are stated
in the individual notes.

The code and the text were prepared with AI assistance (GPT and Claude).
