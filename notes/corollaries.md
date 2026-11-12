# Consequences of published bounds

These are arithmetic consequences of prior theorems or of the explicitly
identified bounds in this package. The single-entry checks verify the
integer hypotheses and substitutions; citations identify any external
mathematical premise.

## Ten consequences of Javadi–Poorhadi–Fallah

Theorem 13 of the 2019 preprint [1] applies when n>k>=4,
k−2 divides n−2, (k−1)(k−2) divides (n−1)(n−2), and
k(k−1)(k−2) does not divide n(n−1)(n−2). It gives

    D(n,k,3) <= J(n,k,3)−3,
    J(n,k,3) = floor(n/k · floor((n−1)/(k−1) · floor((n−2)/(k−2)))).

Here D(n,k,3)=A(n,2(k−2),k). This upper-bound theorem does not assume
n sufficiently large. Its theorem number refers to the preprint, whose
journal version has the different title listed below. The five direct
substitutions and five Johnson consequences are:

| Direct entry | Nested Johnson arithmetic | New bound | Johnson consequence |
|---|---|---:|---|
| A(29,6,5) | floor(29/5 · floor(28/4 · 9)) = 365; 365−3 | 362 | A(30,6,6) <= floor(30·362/6) = 1810 |
| A(38,6,5) | floor(38/5 · floor(37/4 · 12)) = 843; 843−3 | 840 | A(39,6,6) <= floor(39·840/6) = 5460 |
| A(53,6,5) | floor(53/5 · floor(52/4 · 17)) = 2342; 2342−3 | 2339 | A(54,6,6) <= floor(54·2339/6) = 21051 |
| A(52,10,7) | floor(52/7 · floor(51/6 · 10)) = 631; 631−3 | 628 | A(53,10,8) <= floor(53·628/8) = 4160 |
| A(62,10,7) | floor(62/7 · floor(61/6 · 12)) = 1080; 1080−3 | 1077 | A(63,10,8) <= floor(63·1077/8) = 8481 |

## Ten Johnson consequences of the existing table

For any code C, the words containing a fixed coordinate, with that
coordinate removed, form an (n−1,d,w−1) code. Summing the resulting
upper bounds over all coordinates gives w|C|<=nA(n−1,d,w−1), hence

    A(n,d,w) <= floor(n A(n−1,d,w−1)/w).

This is Johnson's bound [2]. The ten substitutions, in table-update order,
are as follows. The first, second, third, fourth, eighth, ninth and tenth
rows use the public table [3]; the other rows propagate bounds from this
same list. The attribution for 124 is Serg Rodriguez, and the table's
observation giving 38 is attributed to János Tapolcai. These are the
table's public attributions, not citations to private correspondence.

| Entry | Premise | Complete arithmetic | New upper bound |
|---|---|---|---:|
| A(34,10,7) | A(33,10,6) <= 32 | floor(34·32/7) | 155 |
| A(37,10,7) | A(36,10,6) <= 37 | floor(37·37/7) | 195 |
| A(47,10,7) | A(46,10,6) <= 67 | floor(47·67/7) | 449 |
| A(33,10,8) | A(32,10,7) <= 124 | floor(33·124/8) | 511 |
| A(35,10,8) | A(34,10,7) <= 155 | floor(35·155/8) | 678 |
| A(38,10,8) | A(37,10,7) <= 195 | floor(38·195/8) | 926 |
| A(48,10,8) | A(47,10,7) <= 449 | floor(48·449/8) | 2694 |
| A(42,12,8) | A(41,12,7) <= 29 | floor(42·29/8) | 152 |
| A(43,12,8) | A(42,12,7) <= 34 | floor(43·34/8) | 182 |
| A(44,12,8) | A(43,12,7) <= 38 | floor(44·38/8) | 209 |

The premise A(41,12,7)<=29 is the table's consequence of the
nonexistence of two orthogonal Latin squares of order six. Its Johnson
consequence is A(42,12,7)<=34. The bound A(43,12,7)<=38 is the
observation attributed there to Tapolcai. No reproof of these external
table premises is claimed by the arithmetic checker.

## Four further substitutions

| Entry | Premise | Complete arithmetic | New upper bound |
|---|---|---|---:|
| A(29,6,8) | A(28,6,7)=4680 in the main table | floor(29·4680/8) | 16965 |
| A(29,6,9) | A(28,6,8)<=12285 in the main table | floor(29·12285/9) | 39585 |
| A(29,6,10) | A(28,6,9)<=31195 in the main table | floor(29·31195/10) | 90465 |
| A(34,10,8) | A(33,10,7)<=136, supplied LP certificate | floor(34·136/8) | 578 |

For the last row, the checker also verifies the rational parent
certificate. The seven consecutive upper bounds for n=29,d=6,w=8,...,14
are therefore 16965,39585,90465,162226,228441,289604,322266. The last
four are documented in [linear_programming.md](linear_programming.md).
The [packing bound A(27,8,7)<=499](packing.md) gives
A(28,8,8)<=floor(28·499/8)=1746 by the same Johnson bound.

References:

1. R. Javadi, E. Poorhadi and F. Fallah, *Packing cliques in 3-uniform hypergraphs*, J. Combin. Des. 28(8) (2020), 580–603. [DOI](https://doi.org/10.1002/jcd.21717). Theorem 13 is cited from *On the packing for triples*, [arXiv:1905.10807v1](https://arxiv.org/pdf/1905.10807v1) (2019).
2. S. M. Johnson, *A new upper bound for error-correcting codes*, IRE Trans. Inform. Theory 8 (1962), 203–207. [DOI](https://doi.org/10.1109/TIT.1962.1057714).
3. A. E. Brouwer, *Bounds for binary constant weight codes*, accessed 10 November 2026. [Table](https://aeb.win.tue.nl/codes/Andw.html).
