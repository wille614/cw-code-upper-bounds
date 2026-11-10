# Exact distance-distribution certificates

For a binary constant-weight code C of weight w<=n/2, let a_i be the
average number of words at distance 2i from a word of C. Thus a_0=1,
a_i>=0, a_i=0 for 0<i<d/2, and |C|=1+sum_(i>=d/2) a_i.
The Johnson-scheme inequalities of Delsarte [1], in the normalization
of Best et al. [2, equation (11), Theorem 21], are

    1 + sum_i q_k(i)a_i >= 0,  1<=k<=w,

where

    q_k(i) = sum_(j=0)^i (-1)^(i-j) binom(w-j,i-j)
                         binom(w-k,j) binom(n-w+j-k,j)
               / [binom(w,i) binom(n-w,i)].

Binomial coefficients with an out-of-range lower argument are zero.
All these numbers are rational and are recomputed by the checker.

Three certificates also use bounds on a distance shell. For a fixed
weight-w word, a word at distance 2i has i removed positions in its
support and i added positions in its complement. The shell is therefore
a doubly constant-weight code with part lengths w,n−w and weights i,i.
Writing T(m,a;v,b;d) for the maximum size of such a code, one has
a_i<=T(w,i;n−w,i;d), as in Agrell–Vardy–Zeger [5].
Complementation within either part preserves
distance, and the two parts may be interchanged. This explains the
normalized parameters in `data/shell_bounds.json`.

The required bounds use elementary disjointness and four shortening
inequalities, for example

    T(m,a;v,b;d) <= floor(m T(m−1,a−1;v,b;d)/a),
    T(m,a;v,b;d) <= floor(m T(m−1,a;v,b;d)/(m−a)),

with analogous inequalities in the second part. Count incidences of
ones or zeros and shorten at a fixed coordinate to prove each inequality.
If 2(a+b)=d, distinct words have disjoint supports, so T is at most
min(floor(m/a),floor(v/b)), omitting a part of weight zero. A zero-weight
part may be discarded. The three additional inequalities are:

| Parameters | Distance half i | Bound on a_i | Complete justification |
|---|---:|---:|---|
| (31,10,7) | 7 | 24 | T(7,0;24,7;10)=A(24,10,7)<=24, the existing table upper bound [3] |
| (33,10,7) | 5 | 93 | T(7,2;24,3;10)<=3; T(7,2;25,4;10)<=floor(25·3/4)=18; T(7,2;26,5;10)<=floor(26·18/5)=93 |
| (29,6,11) | 3 | 864 | T(9,2;16,1;6)<=4; T(10,3;16,1;6)<=floor(10·4/3)=13; T(11,3;16,1;6)<=floor(11·13/8)=17; T(11,3;17,2;6)<=floor(17·17/2)=144; T(11,3;18,3;6)<=floor(18·144/3)=864 |

Thus 106, 136 and 162226 use Delsarte inequalities **with additional
distance-shell bounds**. They are not claimed as values of the
unaugmented Delsarte program. Adding valid distance-distribution
inequalities is classical; see Best–Brouwer [4] and Best et al. [2,
Theorem 22]. The other three certificates use the ordinary Delsarte
inequalities alone.

## Dual certificate format and proof

Write every inequality in the form R_j a<=b_j: a Delsarte row is −q_k
with right side 1; a shell row is the corresponding coordinate vector
with its certified shell cap. A certificate supplies nonnegative rational
multipliers y_j such that sum_j y_j R_j,i>=1 for every permitted distance.
Then

    |C| = 1+sum_i a_i <= 1+sum_j y_j b_j.

Each JSON file gives `n,d,w`, the permitted half-distances, the dual
vector, its exact objective and its floor. An explicit `constraints` list
specifies the nonzero rows in an augmented certificate. Without that
list the dual vector is indexed by k=1,...,w. Unused zero-multiplier
shell rows are omitted. An optional primal vector is supplementary data;
the upper-bound proof depends only on the verified dual inequalities.

| Parameters | Exact rational bound | Floor |
|---|---:|---:|
| (31,10,7) | 2015/19 | 106 |
| (33,10,7) | 6688/49 | 136 |
| (29,6,11) | 8922459/55 | 162226 |
| (29,6,12) | 12564279/55 | 228441 |
| (29,6,13) | 10715355/37 | 289604 |
| (29,6,14) | 6123060/19 | 322266 |

`python verify.py n d w` checks the certificate with exact fractions,
every distance column, and the complete dependency chain of each active
shell cap. The cap A(24,10,7)<=24 is a stated table premise; all other
shell caps above terminate in the proved disjointness case. No optimizer
is needed. The integer floor is the advertised upper bound, without any
claim about exact optimality of a larger relaxation.

References:

1. P. Delsarte, *An algebraic approach to the association schemes of coding theory*, Philips Research Reports Supplements, No.10 (1973), §3.2 and §4.2. [Text](https://users.wpi.edu/~martin/RESEARCH/philips.pdf).
2. M. R. Best, A. E. Brouwer, F. J. MacWilliams, A. M. Odlyzko and N. J. A. Sloane, *Bounds for binary codes of length less than 25*, IEEE Trans. Inform. Theory 24 (1978), 81–93. [Text](https://neilsloane.com/doc/Me54.pdf).
3. A. E. Brouwer, *Bounds for binary constant weight codes*, accessed 10 November 2026. [Table](https://aeb.win.tue.nl/codes/Andw.html).
4. M. R. Best and A. E. Brouwer, *The triply shortened binary Hamming code is optimal*, Discrete Mathematics 17 (1977), 235–245. [DOI](https://doi.org/10.1016/0012-365X(77)90158-3).
5. E. Agrell, A. Vardy and K. Zeger, *Upper bounds for constant-weight codes*, IEEE Trans. Inform. Theory 46(7) (2000), 2373–2395. [DOI](https://doi.org/10.1109/18.887851); [Text](https://publications.lib.chalmers.se/records/fulltext/local_14986.pdf).
