# Divisible leaves and two packing bounds

A simple t-uniform hypergraph H is properly K_k^(t)-divisible if
binom(k−|S|,t−|S|) divides deg_H(S) for every nonempty S with |S|<t.
There is no divisibility condition on its total number E of edges.
Throughout k>=t>=2, Q=binom(k,t), b=binom(k−1,t−1), and
delta_v=deg_H(v)/b.

## Theorem 1

For every vertex v, delta_v <= E/Q. Consequently there is no such
hypergraph with 0<E<Q; at E=Q its nonisolated part is K_k^(t).

Proof. Induct on t. At t=2, if D=deg(v), every neighbor of v has degree
at least k−1. The sum of the degrees on v and its D neighbors is at least
kD. Thus 2E>=kD, which is the required inequality.

For t>=3 the link H_u is a simple (t−1)-uniform, properly
K_(k−1)^(t−1)-divisible hypergraph: the degree of T in the link equals
deg_H(T union {u}). Its edge count is b delta_u. Induction in H_u gives

    delta_uv := deg_H({u,v}) / binom(k−2,t−2) <= delta_u,

and the symmetric inequality with delta_v. Double counting gives
sum_(u!=v) delta_uv = (k−1)delta_v. Hence

    sum_u delta_u >= delta_v + sum_(u!=v) delta_uv = k delta_v,
    sum_u delta_u = tE/b = kE/Q.

This proves the inequality. Every positive delta_v is an integer.
Therefore 0<E<Q is impossible. At E=Q each positive delta_v equals 1
and their sum is k, so the support has exactly k vertices. Simplicity
and E=binom(k,t) then force every t-subset of the support to be an edge.

## Theorem 2

For t>=3, the possible edge counts satisfy E=0, E=Q, or E>=2Q.

Proof. Suppose 0<E<2Q. Theorem 1 makes every positive delta_v equal to 1.
Each nonempty link has binom(k−1,t−1) edges, so the equality case of
Theorem 1 at rank t−1>=2 makes it complete on k−1 vertices. Let N(v)
be that support and W={v} union N(v). For u,x in N(v), u!=x, there is
a (t−1)-subset of N(v) containing both, hence an edge containing v,u,x.
It follows that N(u) contains W\{u}. Both have size k−1, so they are
equal. Every t-subset of W is an edge, by taking the complete link at
any of its vertices. No edge meets both W and its complement, since
that would enlarge a neighbor set. Each nonempty component is therefore
K_k^(t). Since E<2Q, exactly one is possible, and E=Q.

The rank assumption matters: a four-cycle is properly triangle-divisible
and has four edges, between 3 and 6. Simplicity is also essential.

## Packing corollary

Let P be a t-(v,k,1) packing, t>=3. Suppose

    binom(k−s,t−s) divides binom(v−s,t−s),  1<=s<t,
    binom(k,t) does not divide binom(v,t).

Its leave H consists of the uncovered t-subsets and is simple. For every
proper nonempty S,

    deg_H(S) = binom(v−|S|,t−|S|)
               − binom(k−|S|,t−|S|) · #{B in P : S subset B}.

Thus H is properly divisible. Since E=binom(v,t)−Q|P| has nonzero
remainder modulo Q, E is neither 0 nor Q. Theorem 2 gives E>=2Q and

    |P| <= floor(binom(v,t)/binom(k,t)) − 2.

Weight k and distance 2(k−t+1) are equivalent to this packing condition.
For t=4,k=7 the substitutions are:

| v | binom(v,4) | Q | remainder | floor(binom(v,4)/Q)−2 |
|---|---:|---:|---:|---:|
| 27 | 17550 | 35 | 15 | 499 |
| 47 | 178365 | 35 | 5 | 5094 |

For v=27 the three proper divisibility checks are 20|2600, 10|300,
4|24; for v=47 they are 20|15180, 10|990, 4|44. Consequently
A(27,8,7)<=499 and A(47,8,7)<=5094. Johnson's incidence count also gives
A(28,8,8)<=floor(28·499/8)=1746. The single-entry verifier checks each
integer substitution and divisibility condition. The general theorems
are proved above, rather than inferred from finite tests.

Related observations appear in Brouwer [1], pp.89–91. For t=3,
Javadi–Poorhadi–Fallah [2, preprint Theorem 13] give a stronger result:
under the corresponding hypotheses their nested Johnson bound is reduced
by 3. The two t=4 substitutions above use Theorem 2, and the t=3
corollaries in [corollaries.md](corollaries.md) are attributed to [2].

References:

1. A. E. Brouwer, *Packing and covering of (k choose t)-sets*, in A. Schrijver (ed.), *Packing and Covering in Combinatorics*, Mathematical Centre Tracts 106 (1979), 89–97. [Volume](https://ir.cwi.nl/pub/10189/).
2. R. Javadi, E. Poorhadi and F. Fallah, *Packing cliques in 3-uniform hypergraphs*, J. Combin. Des. 28(8) (2020), 580–603. [DOI](https://doi.org/10.1002/jcd.21717). The theorem numbering is from *On the packing for triples*, [arXiv:1905.10807v1](https://arxiv.org/pdf/1905.10807v1), Theorem 13.
3. S. M. Johnson, *A new upper bound for error-correcting codes*, IRE Trans. Inform. Theory 8 (1962), 203–207. [DOI](https://doi.org/10.1109/TIT.1962.1057714).
