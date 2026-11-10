# Nine explicit lower bounds

Each file contains one complete binary constant-weight code. Seven files
use the table's hexadecimal convention: the name ends in H, the first
line is `$BASE=16`, and subsequent lines are zero-padded hexadecimal
integers. The remaining two have one binary word per line. The verifier
checks the declared length, every weight, distinctness, and every
unordered word pair; it prints the complete distance distribution in
single-entry mode. The codes themselves, rather than a construction
heuristic, are the existence certificates.

## Set replacements in seven tabulated codes

In the same coordinate order as the indicated table file, write C_old
for the tabulated code and C for the supplied code. The exact identity
is C=(C_old minus D) union E, with D=C_old\C and E=C\C_old. The
following counts describe these concrete replacements; the coordinate
permutation is the identity throughout. Original file hashes and source
links are in `data/replacements.json`. The source attribution follows
Brouwer's public table [1]. No claim of authorship of the retained words
is implied. Set replacement and the elementary code operations are
discussed in the general constant-weight-code literature [2].

| Parameters | Source attribution | Original code | Retained | Removed | Added | Final |
|---|---|---:|---:|---:|---:|---:|
| (24,6,12) | Moshe Milshtein (Mi) | [5617](https://aeb.win.tue.nl/codes/cwc/d6/a24.6.12.5617H) | 5409 | 208 | 216 | 5625 |
| (25,6,11) | Antoine Olivier (Ol) | [8720](https://aeb.win.tue.nl/codes/cwc/d6/a25.6.11.8720H) | 8434 | 286 | 299 | 8733 |
| (25,6,12) | Pliny the Liberator (Pl) | [9895](https://aeb.win.tue.nl/codes/cwc/d6/a25.6.12.9895H) | 9677 | 218 | 233 | 9910 |
| (26,6,11) | Pliny the Liberator (Pl) | [12222](https://aeb.win.tue.nl/codes/cwc/d6/a26.6.11.12222H) | 12041 | 181 | 184 | 12225 |
| (26,6,13) | Pliny the Liberator (Pl) | [16225](https://aeb.win.tue.nl/codes/cwc/d6/a26.6.13.16225H) | 15464 | 761 | 790 | 16254 |
| (28,8,11) | Antoine Olivier (Ol) | [3791](https://aeb.win.tue.nl/codes/cwc/d8/a28.8.11.3791H) | 3393 | 398 | 439 | 3832 |
| (28,8,14) | Pliny the Liberator (Pl) | [6359](https://aeb.win.tue.nl/codes/cwc/d8/a28.8.14.6359H) | 6012 | 347 | 380 | 6392 |

## A(30,10,10) >= 373

Over F_3, g(x)=2+x²+2x³+x⁴+x⁵ is a generator polynomial of the
ternary Golay code. Indeed, g divides x¹¹−1, with quotient
1+x²+2x³+2x⁴+2x⁵+x⁶. The products b(x)g(x), deg(b)<=5, form
the cyclic [11,6,5] ternary Golay code. Restricting its last coordinate
to zero and deleting that coordinate gives the once-shortened ternary
Golay [10,5,5] code: its words are exactly b(x)g(x), deg(b)<=4.
Multiplication is injective because g is nonzero. The complete weight
distribution of these 243 words is

| Weight | 0 | 5 | 6 | 8 | 9 |
|---|---:|---:|---:|---:|---:|
| Number | 1 | 72 | 60 | 90 | 20 |

The unit-vector substitution maps a symbol a_i to the coordinate 3i+a_i
in a block of three, including the symbol zero. Distance doubles and
each resulting word has weight ten. This is the classical construction
of Kautz–Singleton [3]; see also Fu–Han Vinck–Shen [4]. It supplies 243
binary words. Add the 130 words listed in `data/A30_ternary.json` as
`added_hex`. The hexadecimal integers use bit position zero at the right,
so a support S is encoded by sum_(j in S) 2^j. Every added word occupies
at most five triples and therefore intersects any unit-vector word in
at most five positions. The 27300 six-subsets of the added supports are
all distinct, proving that no two added words intersect in six positions.
The two parts are disjoint, hence give 243+130=373 words at distance
at least ten. The verifier checks both the subset argument and the full
distance distribution, and compares this construction with the code file.

## A(32,12,11) >= 215

Put h(x)=(x−1)g(x)=1+2x+2x²+2x³+x⁴+x⁶ over F_3. The products
b(x)h(x), deg(b)<=4, form the [11,5,6] dual of the ternary Golay
code above. Each of the five shifted copies of h is orthogonal over
F_3 to each of the six shifted copies of g; their dimensions, five
and six, sum to eleven, proving equality with the dual code.
Its complete weight distribution is

| Weight | 0 | 6 | 9 |
|---|---:|---:|---:|
| Number | 1 | 132 | 110 |

Thus its minimum distance is six. Apply the same unit-vector map, giving
length 33, weight 11 and distance at least 12. Keep exactly the words
with first ternary symbol nonzero. Their binary coordinate numbered zero
(the position for symbol zero in the first triple) is then zero. Delete
that binary coordinate. There are 162 words, since each of the two
nonzero first symbols occurs 81 times. This is **shortening at a zero
binary coordinate after expansion**, not ordinary ternary shortening.
Weight and distance are preserved and the new length is 32. Original
binary coordinate j>=1 becomes coordinate j−1.

All 162 words of this core are present in the final file. Add the 53
words in `data/A32_ternary.json`, again encoded by hexadecimal support
integers with bit zero at the right. No core word is deleted. The full
pairwise check verifies the distances within and between both parts;
it also confirms their disjointness and their equality with the supplied
215-word set. This finite certificate establishes the lower bound without
requiring a general rule about all possible added words.

## Final files

| Parameters | Number of words | File |
|---|---:|---|
| (24,6,12) | 5625 | [`codes/a24.6.12.5625H`](../codes/a24.6.12.5625H) |
| (25,6,11) | 8733 | [`codes/a25.6.11.8733H`](../codes/a25.6.11.8733H) |
| (25,6,12) | 9910 | [`codes/a25.6.12.9910H`](../codes/a25.6.12.9910H) |
| (26,6,11) | 12225 | [`codes/a26.6.11.12225H`](../codes/a26.6.11.12225H) |
| (26,6,13) | 16254 | [`codes/a26.6.13.16254H`](../codes/a26.6.13.16254H) |
| (28,8,11) | 3832 | [`codes/a28.8.11.3832H`](../codes/a28.8.11.3832H) |
| (28,8,14) | 6392 | [`codes/a28.8.14.6392H`](../codes/a28.8.14.6392H) |
| (30,10,10) | 373 | [`codes/a30.10.10.373`](../codes/a30.10.10.373) |
| (32,12,11) | 215 | [`codes/a32.12.11.215`](../codes/a32.12.11.215) |

References:

1. A. E. Brouwer, *Bounds for binary constant weight codes*, accessed 10 November 2026. [Table](https://aeb.win.tue.nl/codes/Andw.html). The author labels identify the public source of each older code.
2. A. E. Brouwer, J. B. Shearer, N. J. A. Sloane and W. D. Smith, *A new table of constant weight codes*, IEEE Trans. Inform. Theory 36 (1990), 1334–1380. [DOI](https://doi.org/10.1109/18.59932).
3. W. H. Kautz and R. C. Singleton, *Nonrandom Binary Superimposed Codes*, IEEE Trans. Inform. Theory 10 (1964), 363–377, §V.B, pp.370–371 (unit-vector expansion and supplementary words). [DOI](https://doi.org/10.1109/TIT.1964.1053689).
4. F.-W. Fu, A. J. Han Vinck and S.-Y. Shen, *On the Constructions of Constant-Weight Codes*, IEEE Trans. Inform. Theory 44 (1998), 328–333, Proposition 3.3(1). [DOI](https://doi.org/10.1109/18.651060).
