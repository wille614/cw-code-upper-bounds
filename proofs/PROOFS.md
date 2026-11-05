# Validity of the augmented four-point bounds

This document gives the mathematical justification for the constraints in the augmented four-point semidefinite programs for binary constant-weight codes with parameters $(18,6,7)$, $(19,6,7)$, and $(19,6,8)$. It explains why the moment vector of every such code is feasible and why the exact rational dual certificates therefore bound its cardinality.

## 1. Codes, the Johnson scheme, and moment normalization

Write $[n]=\{1,\ldots,n\}$ and identify a binary word of weight $w$ with its support in

$$
\mathcal F=\binom{[n]}w,\qquad V=|\mathcal F|=\binom nw.
$$

For $u,v\in\mathcal F$,

$$
d_H(u,v)=2\bigl(w-|u\cap v|\bigr).
$$

The relations indexed by $w-|u\cap v|$ form the Johnson scheme. A code is a subset $C\subseteq\mathcal F$ whose distinct words have Hamming distance at least $d$. Its largest possible cardinality is $A(n,d,w)$. Put $\delta=\lceil d/2\rceil$. Call a set of words **compatible** if it satisfies this distance condition.

For an arbitrary fixed code $C$, define

$$
X_u(C)=\mathbf 1_{\{u\in C\}},\qquad
X_S=\prod_{u\in S}X_u,\qquad X_\varnothing=1.
$$

Products are reduced using $X_u^2=X_u$ and $X_uX_v=0$ when distinct $u,v$ are incompatible. Thus $X_SX_T=X_{S\cup T}$, with an incompatible union interpreted as zero.

Let $G=S_n$ act by coordinate permutations and let $g$ be uniformly distributed in $G$. The scaled, invariant moments used here are

$$
z_C(S)=V\,\mathbb E_g[X_S(gC)],\qquad |S|\le4.
$$

They depend only on the $G$-orbit of the unordered set $S$. Incompatible moments are zero. In particular,

$$
z_C(\varnothing)=V,\qquad
z_C(\{u\})=|C|\quad(u\in\mathcal F),\qquad
0\le z_C(S)\le V.
$$

Write $z_1$ for the common singleton moment. The objective of the invariant program is $z_1$. Polak's unscaled moments are $x(S)=z(S)/V$, so $x(\varnothing)=1$ and $\sum_{u\in\mathcal F}x(\{u\})=z_1$.

**Constant convention.** The affine constant in the model is $1$. The notation $z(\varnothing)$ in this document means the fixed number $V$, not that affine constant. For example, a term $a\,z(\varnothing)$ contributes the constant coefficient $aV$.

For a polynomial whose Boolean reduction has degree at most four, define the linear moment functional

$$
\mathcal L_z\!\left(\sum_S a_SX_S\right)=\sum_S a_Sz(S).
$$

Consequently $\mathcal L_z(1)=V$ and

$$
\mathcal L_{z_C}(p)=V\,\mathbb E_g[p(gC)].
$$

Repeated words are always merged before applying this functional.

## 2. Polak's constraints and positivity under symmetry

### 2.1 The original four-point constraints

For an anchor $D\subseteq\mathcal F$, put

$$
\mathcal I_D=\{S\subseteq\mathcal F:D\subseteq S,\ |D|+2|S\setminus D|\le4\},
\qquad
M_D(z)_{S,T}=z(S\cup T).
$$

The following three model families are Polak's original constraints, with the normalization of Section 1:

| Model label | Constraint before block decomposition | Source |
| --- | --- | --- |
| `D0` | $M_\varnothing(z)\succeq0$ | [1], Sections 1.1 and 3.2 |
| `D2` | $M_D(z)\succeq0$ for compatible $D$ with $\lvert D\rvert=2$ | [1], Sections 1.1 and 3.1 |
| `nonnegative` | $z(S)\ge0$ for nonempty compatible sets with $\lvert S\rvert\le4$ | [1], final paragraphs of Section 1.1 |

Here and below $\succeq0$ means positive semidefinite. Polak's Section 1.1, equations (1)–(2) and Proposition 1.1, establish these constraints and their upper-bound property; they are used here without reproving them. The same section explains the sufficiency of even anchor sizes. Forbidden sets have moment zero throughout.

With $z(\varnothing)=V$ and orbit invariance, maximizing $z_1$ subject to these original constraints is Polak's four-point bound $A_4(n,d,w)$.

Polak's Section 1.2 justifies averaging over $S_n$ and replacing moments by orbit variables. Sections 2 and 2.1 describe fixed representation matrices and Young representatives; Sections 3.1 and 3.2 give the anchored and empty-anchor block reductions. These are the sources for the symmetry reduction of the underlying four-point program.

### 2.2 Averaging and projections for the added constraints

The added constraints use the following elementary principle. If $p(C)\ge0$ for every code, then $\mathcal L_{z_C}(p)\ge0$. If $h(C)\ge0$ and $f$ is a real vector of polynomials, then

$$
\mathcal L_{z_C}(hff^{\mathsf T})\succeq0,
$$

provided all entries have Boolean degree at most four. Indeed, for any real vector $a$ its quadratic form is

$$
V\,\mathbb E_g\!\left[h(gC)\bigl(a^{\mathsf T}f(gC)\bigr)^2\right]\ge0.
$$

For any fixed real matrix $P$, positive semidefiniteness implies

$$
K\succeq0\quad\Longrightarrow\quad P^{\mathsf T}KP\succeq0,
$$

because $a^{\mathsf T}P^{\mathsf T}KPa=(Pa)^{\mathsf T}K(Pa)\ge0$. This applies to Young projections, to selecting rows and columns, and to replacing word indicators by fixed sums of them. Averaging valid inequalities or positive semidefinite matrices also preserves positivity. Division by a positive common factor is harmless.

Thus selecting a subfamily of valid constraints, deleting zero rows, or retaining only some representation blocks cannot exclude the moment vector of a code. No completeness assertion about a selected subfamily is needed for the upper-bound argument. The general symmetry viewpoint is also treated in [5], Theorem 3.3 and Section 4.

## 3. Anticode vacancy scalars: `ACV-scalar`

Put $t=w-\delta+1$. For a coordinate set $R$ with $|R|=t+2r$, define

$$
Q(R,r)=\{u\in\mathcal F:|u\cap R|\ge t+r\}.
$$

Here $t\ge0$, $0\le r\le\delta-1$, and only feasible set sizes are considered. At the three parameters in this package, $\delta=3$ and $r\in\{0,1,2\}$. These are classical intersecting families; see [6], equation (1.9).

For $u,v\in Q(R,r)$,

$$
|u\cap v|\ge |u\cap R|+|v\cap R|-|R|\ge t,
\qquad d_H(u,v)\le2(\delta-1)<d.
$$

Hence every code contains at most one word of $Q=Q(R,r)$. The vacancy polynomial

$$
H_Q=1-\sum_{q\in Q}X_q
$$

therefore takes values in $\{0,1\}$ on every code.

For each retained compatible three-word anchor $D$, the scalar constraint is

$$
z(D)-\sum_{q\in Q}z(D\cup\{q\})\ge0.
$$

**Proof for every code.** The left side at $z_C$ is $\mathcal L_{z_C}(X_DH_Q)$. Both factors are nonnegative on every $gC$, so the inequality follows. The union has at most four words. Incompatible unions contribute zero, and repetitions use $X_q^2=X_q$.

If $D\cap Q\ne\varnothing$, then $X_DH_Q$ vanishes on every code. The retained scalar subfamily omits these cases, requires at least two compatible extension words, and removes duplicate rows. None of these omissions affects validity. The same proof would allow smaller anchors, but the scalar rows described here use the three-word anchors of this model family.

Coordinate permutations preserve the construction. Choosing representatives of the positions of $R$ relative to $D$ therefore gives valid inequalities in the orbit variables.

## 4. Anticode vacancy matrices: `ACV-matrix`

Fix $Q=Q(R,r)$ as above and a word $a\in\mathcal F$. Let

$$
\mathcal B_a=\{b\in\mathcal F:b=a\text{ or }d_H(a,b)\ge d\}.
$$

Define the matrix indexed by $b,c\in\mathcal B_a$ by

$$
K^{a,Q}(z)_{b,c}
=z(\{a,b,c\})-\sum_{q\in Q}z(\{a,b,c,q\}).
$$

The family `ACV-matrix` consists of its retained fixed projections

$$
P^{\mathsf T}K^{a,Q}(z)P\succeq0.
$$

**Proof for every code.** With $v=(X_b)_{b\in\mathcal B_a}$,

$$
K^{a,Q}(z_C)=V\,\mathbb E_g[X_aH_Qvv^{\mathsf T}]\succeq0.
$$

This follows from $X_aH_Q\ge0$. Each entry uses at most four distinct words. The index set requires compatibility with $a$; it does not require $b$ and $c$ to be compatible with one another. If they conflict, their entry is zero. The coordinate $b=a$ serves as the constant coordinate after multiplication by $X_a$.

For clarity, the projection in these models can be described explicitly. The stabilizer of $(a,R)$ permutes independently the four coordinate cells determined by membership in $a$ and $R$. At harmonic degree zero, words are grouped by their four cell weights. Let $\mathcal O_\beta$ be one such word orbit and put

$$
P_{b,\beta}=\mathbf1_{\{b\in\mathcal O_\beta\}}.
$$

An equivalent projected matrix has entries

$$
\sum_{b\in\mathcal O_\beta}\sum_{c\in\mathcal O_\gamma}
K^{a,Q}(z)_{b,c}.
$$

Degree-zero Young representatives give these columns up to fixed nonzero normalizations, which are further congruences. Only harmonic degree zero is used in this family here. For any fixed Young projection the positivity argument of Section 2.2 applies.

Rows indexed by words in $Q$ vanish on every code under the factor $H_Q$ and may be removed. If $a\in Q$, the entire matrix vanishes. Taking representatives under coordinate permutations, removing zero rows, and positive rescaling all preserve the claimed validity.

## 5. Matching-region occupancy projectors: `matching-occupancy`

Choose coordinate sets $L\subseteq U\subseteq[n]$ and set

$$
Q=\{u\in\mathcal F:L\subseteq u\subseteq U\},\qquad
m=|U\setminus L|,\qquad k=w-|L|.
$$

This family uses

$$
\min(k,m-k)=\delta,\qquad 3\delta\le m<4\delta.
$$

For words of $C\cap Q$, deletion of the fixed coordinates leaves $k$-subsets of an $m$-set with the same mutual distances. If $k=\delta$, two residual sets must be disjoint: their distance is at least $d$ and is even, so it is at least $2\delta$. If $m-k=\delta$, the same statement holds for the residual complements. Thus

$$
N=|C\cap Q|\le\left\lfloor\frac m\delta\right\rfloor=3.
$$

Define $e_0=1$ and

$$
e_\ell(X_Q)=
\sum_{\substack{E\subseteq Q,\ |E|=\ell\\E\text{ compatible}}}X_E.
$$

On a code, $e_\ell(X_Q)=\binom N\ell$. For $0\le j\le3$ put

$$
P_j=\sum_{\ell=j}^3(-1)^{\ell-j}\binom\ell j e_\ell(X_Q).
$$

If $N<j$ this is zero. Otherwise, the identity $\binom\ell j\binom N\ell=\binom Nj\binom{N-j}{\ell-j}$ gives

$$
P_j(C)=\binom Nj\sum_{r=0}^{N-j}(-1)^r\binom{N-j}r
=\mathbf1_{\{N=j\}}.
$$

The four polynomials are

$$
P_0=1-e_1+e_2-e_3,\quad
P_1=e_1-2e_2+3e_3,\quad
P_2=e_2-3e_3,\quad P_3=e_3.
$$

For a word anchor $a$, the retained constraints are

$$
\sum_{\ell=j}^3(-1)^{\ell-j}\binom\ell j
\sum_{\substack{E\subseteq Q,\ |E|=\ell\\E\text{ compatible}}}
z(\{a\}\cup E)\ge0,\qquad j=0,1,2.
$$

**Proof for every code.** This is $\mathcal L_{z_C}(X_aP_j)\ge0$, because $X_aP_j$ is a product of indicators. Every union contains at most four words. If $a\in Q$, it remains eligible in $E$; its repeated occurrence is merged by set union. Excluding it would change the polynomial.

The $j=3$ inequality is already a nonnegative sum of moments. It, zero rows, duplicates, and other rows with only nonnegative coefficients can be omitted. Region and anchor representatives are taken under coordinate permutations; averaging preserves every retained inequality.

## 6. Coordinate interval spectra: `coordinate-spectrum`

### 6.1 An elementary integer capacity

For $0\le k\le m$, define the shortening bound $J(m,d,k)$ recursively by

$$
J(m,d,k)=1\quad\text{if }\min(k,m-k)<\delta,
$$

and otherwise by

$$
J(m,d,k)=\min\left\{
\left\lfloor\frac{mJ(m-1,d,k-1)}k\right\rfloor,
\left\lfloor\frac{mJ(m-1,d,k)}{m-k}\right\rfloor
\right\}.
$$

**Proof of the capacity bound.** In the base case, the distance between two constant-weight words is at most $2\min(k,m-k)<d$, so a code has at most one word. Otherwise, for a length-$m$, weight-$k$ code $E$, let $E_i^1$ and $E_i^0$ be its words having a one or zero, respectively, at coordinate $i$. Deleting that fixed coordinate preserves distances. Induction on $m$ gives

$$
|E_i^1|\le J(m-1,d,k-1),\qquad
|E_i^0|\le J(m-1,d,k).
$$

Counting incidences yields

$$
k|E|=\sum_i|E_i^1|\le mJ(m-1,d,k-1),
$$

$$
(m-k)|E|=\sum_i|E_i^0|\le mJ(m-1,d,k).
$$

Taking integer parts proves the recursion. The cases $k=0$ and $k=m$ belong to the base case, so no division by zero occurs.

Now choose disjoint coordinate sets $L,Z$ of sizes $\ell,b$, respectively, and define

$$
Q(L,Z)=\{u\in\mathcal F:L\subseteq u,\ u\cap Z=\varnothing\}.
$$

The restricted code has residual length $m=n-\ell-b$ and weight $k=w-\ell$. Therefore its occupancy satisfies

$$
N_{L,Z}=\sum_{u\in Q(L,Z)}X_u\in\{0,1,\ldots,c\},
\qquad c=J(m,d,k).
$$

Only feasible intervals with $2\le c\le20$ are selected.

### 6.2 The selected nonnegative polynomials

For a selected degree $1\le q\le\min(4,c)$, choose a $q$-element subset $\mathcal R\subseteq\{0,\ldots,c\}$ and a sign $s\in\{1,-1\}$ such that

$$
p(x)=s\prod_{r\in\mathcal R}(x-r)\ge0
\quad\text{for every }x\in\{0,\ldots,c\}.
$$

These are a direct finite description of the relevant polynomial inequalities. To explain the selection, the cone of polynomials of degree at most $q$ nonnegative on $\{0,\ldots,c\}$ is pointed and has dimension $q+1$: evaluation is injective, and the constant polynomial $1$ satisfies all inequalities strictly. An extreme ray has at least $q$ independent active evaluation constraints. A nonzero polynomial of degree at most $q$ has at most $q$ distinct roots, so each extreme ray is a positive multiple of one of the displayed signed products. Conversely, its $q$ distinct roots give $q$ independent conditions, making the associated face one-dimensional. Lower-degree valid products are also allowable.

Expand a selected polynomial in the binomial basis:

$$
p(x)=\sum_{r=0}^q a_r\binom xr,\qquad
a_r=\Delta^r p(0),\qquad \Delta p(x)=p(x+1)-p(x).
$$

All $a_r$ are integers. The Boolean polynomial representing $\binom{N_{L,Z}}r$ is the sum of $X_T$ over compatible $r$-subsets of $Q(L,Z)$. Consequently these inequalities use at most four-word moments.

### 6.3 Exact orbit formula and its proof

Let $\Omega_r$ denote the orbits of compatible unordered $r$-sets of words. For a representative $T$, write

$$
o(T)=|\operatorname{orb}(T)|,\qquad
I(T)=\left|\bigcap_{u\in T}u\right|,\qquad
U(T)=\left|\bigcup_{u\in T}u\right|.
$$

There are $H_{\ell,b}=\binom n\ell\binom{n-\ell}b$ ordered disjoint pairs $(L,Z)$ of the specified sizes. The constraint with label `coordinate-spectrum` is

$$
a_0H_{\ell,b}V+
\sum_{r=1}^q a_r\sum_{T\in\Omega_r}
o(T)\binom{I(T)}\ell\binom{n-U(T)}b\,z(T)\ge0.
$$

An impossible binomial coefficient is interpreted as zero.

**Proof for every code.** For each fixed $(L,Z)$, the occupancy is an integer in $\{0,\ldots,c\}$, so $p(N_{L,Z})\ge0$. Apply $\mathcal L_{z_C}$ and sum over all such pairs. A fixed nonempty word set $T$ lies in the interval precisely when $L$ is chosen in its common intersection and $Z$ is chosen outside its union. There are exactly

$$
\binom{I(T)}\ell\binom{n-U(T)}b
$$

such choices. Grouping the word sets by orbit gives the displayed coefficients. The constant term is $a_0H_{\ell,b}V$ because $\mathcal L_z(1)=V$. This proves both the sign and the normalization of the row. Taking orbit representatives, deleting repeated rows, and positive rescaling preserve validity.

## 7. Coordinate occupancy matrices: `occupancy-matrix`

### 7.1 The localizing matrices, including complements

Choose an orientation $\varepsilon\in\{0,1\}$ and write

$$
u^{(0)}=u,\qquad u^{(1)}=[n]\setminus u,\qquad
w_0=w,\quad w_1=n-w.
$$

Complementation preserves Hamming distance. The indicator of the complemented word $u^{(1)}$ in the complemented code is still $X_u$; it is not $1-X_u$.

Fix a nonempty coordinate core $R$ of size $\ell\le w_\varepsilon$ and define

$$
N_R=\sum_{u:R\subseteq u^{(\varepsilon)}}X_u,\qquad
F_{R,i}=\sum_{u:|u^{(\varepsilon)}\cap R|=i}X_u,
\qquad v_R=(1,(F_{R,i})_i)^{\mathsf T}.
$$

The row indices $i$ range over

$$
\max(0,\ell+w_\varepsilon-n)\le i\le\min(\ell,w_\varepsilon).
$$

By Section 6.1,

$$
0\le N_R\le c,\qquad c=J(n-\ell,d,w_\varepsilon-\ell).
$$

The selected cores satisfy $2\le c\le20$. For each such orientation and core size the family uses

$$
p_j(x)=(x-j)(x-j-1)\quad(0\le j<c),\qquad
p_{\mathrm{end}}(x)=x(c-x).
$$

A precise scaled form of each matrix constraint is

$$
\mathcal O_{\varepsilon,\ell,p}(z)
=\sum_{\substack{R\subseteq[n]\\|R|=\ell}}
\mathcal L_z\bigl(p(N_R)v_Rv_R^{\mathsf T}\bigr)\succeq0.
$$

**Proof for every code.** If $N_R$ is an integer, then $N_R-j$ and $N_R-j-1$ cannot have opposite signs, so $p_j(N_R)\ge0$. The capacity bound also gives $p_{\mathrm{end}}(N_R)\ge0$. Each summand is therefore an average of nonnegative scalar multiples of outer products. Every entry has Boolean degree at most $2+1+1=4$.

The sum over all cores is invariant under coordinate permutations. Thus on code moments its normalization is exactly

$$
\mathcal O_{\varepsilon,\ell,p}(z_C)
=V\sum_{|R|=\ell}p(N_R(C))v_R(C)v_R(C)^{\mathsf T}
=V\binom n\ell\,\mathbb E_R[p(N_R(C))v_R(C)v_R(C)^{\mathsf T}].
$$

Both factors $V$ and $\binom n\ell$ are present. They are positive and do not affect positive semidefiniteness.

For these parameters the selected pairs $(\ell,c)$ are as follows:

| Parameters | Original orientation | Complement orientation | Number of matrices |
| --- | --- | --- | ---: |
| $(18,6,7)$ | $(3,15),(4,4)$ | $(6,19),(7,8),(8,3)$ | 54 |
| $(19,6,7)$ | $(3,20),(4,5)$ | $(7,19),(8,8),(9,3)$ | 60 |
| $(19,6,8)$ | $(4,15),(5,4)$ | $(7,9),(8,3)$ | 35 |

Each selected pair contributes $c$ adjacent-root polynomials and one endpoint polynomial.

### 7.2 Exact coefficients and repeated words

The entries above can be written directly in the existing orbit variables. Expand $N_R^p(v_R)_i(v_R)_j$ for $p=0,1,2$, taking the constant row to contribute no word slot. There are at most four ordered slots: $p$ region slots and zero, one, or two row slots. Repetitions are allowed.

Each tuple determines a unique underlying set $T$ of distinct words and a surjection from its slots onto $T$. Conversely, a set and such a surjection determine exactly one tuple. Incompatible $T$ contribute zero. This accounts for the diagonal terms in $N_R^2$ and for coincidences between row and region words.

For a representative $T$ and a fixed surjection, the core must lie in the intersection of all oriented region words; when there is no region slot, this eligible set is $[n]$. Partition it into four cells according to the two oriented row-word bits, using bit zero for a missing row. If the cell sizes are $c_{00},c_{01},c_{10},c_{11}$, the number of allowed cores is

$$
\sum_{(x_{00},x_{01},x_{10},x_{11})}
\prod_{a,b\in\{0,1\}}\binom{c_{ab}}{x_{ab}},
$$

where the sum imposes $\sum_{a,b}x_{ab}=\ell$ and, for present row slots, the required intersections $x_{10}+x_{11}=i$ and $x_{01}+x_{11}=j$. Multiply this count by $o(T)$, sum over the surjections, and use $z(T)$. The empty tuple contributes $\binom n\ell V$. Combining these three powers of $N_R$ with the coefficients of $p$ gives exactly $\mathcal O_{\varepsilon,\ell,p}(z)$.

For the complement orientation all coordinate bits are reversed before these counts are formed. The moment remains that of the original set $T$, with the same orbit size, so no additional variables are needed. The vector $v_R$ is itself a fixed linear image of the constant and individual-word indicators. Its use is therefore a congruence projection of a localizing matrix. Zero-row deletion and positive common-factor removal preserve validity.

## 8. The independently certified cardinality premise

Let $U_{\mathrm{base}}$ be the exact upper bound supplied by the baseline certificate for the unaugmented Polak four-point program. Every code is feasible for that program, so

$$
|C|\le U_{\mathrm{base}},\qquad
|C|\le B:=\lfloor U_{\mathrm{base}}\rfloor.
$$

The three baseline certificates in this package give, respectively,

$$
U_{\mathrm{base}}(18,6,7)=
\frac{13433029372173119186202113540456444959822875841641987}
{42949672960000000000000000000000000000000000000000},
\qquad B=312,
$$

$$
U_{\mathrm{base}}(19,6,7)=
\frac{1274550274376814049016394972588387247712471375324415947}
{2748779069440000000000000000000000000000000000000000},
\qquad B=463,
$$

$$
U_{\mathrm{base}}(19,6,8)=
\frac{841300017883100309254921891583162664400313}
{1250000000000000000000000000000000000000},
\qquad B=673.
$$

The floor is justified by the integer cardinality of a code. It is not an assertion that every feasible point of the original SDP has objective at most that floor. The baseline certificates establish the premise independently of the three cardinality families below.

Write

$$
M(X)=\sum_{u\in\mathcal F}X_u.
$$

For every code, including the empty code and codes smaller than $B$, the polynomial $B-M(X)$ is nonnegative. No assumption of equality $|C|=B$ is used.

## 9. Global cardinality cap: `cardinality-cap`

The constraint is

$$
B-z_1\ge0.
$$

**Proof for every code.** At $z_C$ this is $B-|C|\ge0$, established in Section 8. Equivalently,

$$
\frac1V\mathcal L_z(B-M(X))
=B-\frac1V\sum_{u\in\mathcal F}z(\{u\})=B-z_1.
$$

This also explains the affine constant $B$ in the scalar row. The expression is already invariant under coordinate permutations.

## 10. Anchored cardinality scalars: `cardinality-scalar`

For each retained nonempty compatible anchor $D$, with $k=|D|\le3$, the constraint is

$$
Bz(D)-\sum_{u\in\mathcal F}z(D\cup\{u\})\ge0.
$$

Equivalently, after separating repeated anchor words,

$$
(B-k)z(D)-\sum_{u\in\mathcal F\setminus D}z(D\cup\{u\})\ge0.
$$

**Proof for every code.** The first expression is $\mathcal L_{z_C}(X_D(B-M(X)))$. Its integrand is nonnegative by Section 8, and its degree is at most four. This applies even when an anchor word is absent from $C$, in which case the corresponding unaveraged product is zero.

For an exact orbit description, let $T$ be a compatible $k$-set and $W$ a compatible $(k+1)$-set orbit representative. Let

$$
m(T,W)=\#\{T'\subseteq W:|T'|=k,\ T'\in\operatorname{orb}(T)\},
\qquad
e(T,W)=\frac{o(W)m(T,W)}{o(T)}.
$$

This is an integer. Count incidences $(T',W')$ with $T'\subseteq W'$ and the specified orbit types: there are $o(W)m(T,W)$ incidences, and transitivity makes their number above each $T'$ equal. For fixed $T$, the containing sets $W'$ correspond exactly to compatible new words extending $T$. Hence the scalar row is

$$
(B-k)z(T)-\sum_{W\in\Omega_{k+1}}e(T,W)z(W)\ge0.
$$

In particular, for code-induced moments,

$$
|C|\,z_C(T)=kz_C(T)+\sum_{W\in\Omega_{k+1}}e(T,W)z_C(W).
$$

The last identity uses the fact that $|gC|=|C|$. On an arbitrary SDP vector, the right side denotes the linearized moment of $M(X)X_T$; it is not identified with the nonlinear product $z_1z(T)$. Orbit averaging and the stated row reductions preserve validity.

## 11. Rooted cardinality matrices: `cardinality-rooted`

Fix a word $a$ and use $\mathcal B_a$ from Section 4. Define

$$
K^{a,B}(z)_{b,c}
=Bz(\{a,b,c\})-\sum_{u\in\mathcal F}z(\{a,b,c,u\}),
\qquad b,c\in\mathcal B_a.
$$

The constraints are the retained Young projections

$$
P^{\mathsf T}K^{a,B}(z)P\succeq0.
$$

**Proof for every code.** If $v=(X_b)_{b\in\mathcal B_a}$, then

$$
K^{a,B}(z_C)
=V\,\mathbb E_g[X_a(B-M(X))vv^{\mathsf T}]\succeq0.
$$

The scalar multiplier is nonnegative on every $gC$. Equivalently, for this fixed code,

$$
K^{a,B}(z_C)=(B-|C|)\,\mathcal L_{z_C}(X_avv^{\mathsf T}).
$$

Every entry has at most four distinct words. Repeated words, including $b=a$ or $c=a$, are merged by union; a missing anchor gives a zero unaveraged matrix.

The stabilizer of $a$ is $S_w\times S_{n-w}$. The projections use tensor products of binary Young representatives on these two coordinate cells, with shapes $(w-h,h)$ and $(n-w-j,j)$, where $h,j$ are nonnegative integers with $2h\le w$ and $2j\le n-w$. The word-weight indices satisfy

$$
h\le s\le w-h,\qquad
j\le w-s\le n-w-j,\qquad
s=w\ \text{or}\ w-s\ge\delta,
$$

where $s=|a\cap b|$. These are fixed coefficient vectors, so Section 2.2 proves positivity of every retained block. The Young construction is the singleton-anchor case of [1], Section 3.1, applied here to the cardinality localizer.

For $T=\{a,b,c\}$ compatible and $k=|T|$, the orbit expression for each entry is exactly

$$
(B-k)z(T)-\sum_{W\in\Omega_{k+1}}e(T,W)z(W),
$$

using the incidence counts of Section 10. Since the congruence matrices are fixed, this linear substitution commutes with projection. If $T$ is incompatible, the entire entry is zero. Positive rescaling and deletion of zero rows preserve the constraint.

## 12. Feasibility, exact dual certificates, and the three bounds

The model labels are exhausted by the following list:

| Model label | Mathematical justification |
| --- | --- |
| `D0` | Polak's empty-anchor matrix; Section 2.1 |
| `D2` | Polak's two-anchor matrices; Section 2.1 |
| `nonnegative` | Polak's nonnegative orbit moments; Section 2.1 |
| `ACV-scalar` | Three-word-anchor vacancy scalars; Section 3 |
| `ACV-matrix` | One-word-anchor vacancy matrices; Section 4 |
| `matching-occupancy` | Capacity-three occupancy projectors; Section 5 |
| `coordinate-spectrum` | Nonnegative polynomials on integer interval occupancies; Section 6 |
| `occupancy-matrix` | Coordinate occupancy localizers in both orientations; Section 7 |
| `cardinality-cap` | Certified global cap; Sections 8–9 |
| `cardinality-scalar` | Nonempty-anchor cardinality localizers; Section 10 |
| `cardinality-rooted` | One-word-anchor cardinality matrix projections; Section 11 |

Let $C$ be any code with one of the three parameters. Polak's result makes $z_C$ feasible for the original four-point constraints. Sections 3–11 show that the same vector satisfies every added constraint. All moments used have at most four distinct words. Coordinate averaging, fixed congruences, and the selected positive scalings preserve these facts. The objective at this feasible point is $z_1=|C|$.

An exact rational dual certificate proves an upper bound on the objective of **every point satisfying all the constraints of the augmented SDP**, not only on code-induced moment vectors. This is weak duality with exact arithmetic: the dual matrices are positive semidefinite, and their inner products with feasible primal blocks are nonnegative. The recorded residual compensation bounds any remaining coefficient residual over the valid moment bounds. Consequently the final certified number, including that compensation, bounds the objective. No assertion that a numerical optimizer is exact, or that strong duality holds, is required.

The moment bounds $0\le z(S)\le V$ used for this purpose also hold throughout the underlying invariant four-point relaxation. To see the upper bounds, apply the $2\times2$ principal-minor inequality successively to the indices $\varnothing,\{u\}$ in $M_\varnothing$, then $\{u\},\{u,v\}$ in $M_\varnothing$, and then $D,D\cup\{w\}$ in $M_D$ for $|D|=2$. These give the bounds for singleton, pair, and triple moments, respectively. A quadruple entry between $D\cup\{w\}$ and $D\cup\{x\}$ is at most the geometric mean of the two bounded triple diagonal entries. Nonnegativity handles zero diagonal entries. These are consequences of the cited Polak constraints, so the compensation applies to all feasible points.

The three augmented certificates give the following exact rational numbers:

| Parameters | Certified rational upper bound | Exact decimal | Integer upper bound |
| --- | --- | --- | ---: |
| $(18,6,7)$ | $\dfrac{3868791767977833258518625959723}{12500000000000000000000000000}$ | $309.50334143822666068149007677784$ | 309 |
| $(19,6,7)$ | $\dfrac{46047385392564732804912812389921}{100000000000000000000000000000}$ | $460.47385392564732804912812389921$ | 460 |
| $(19,6,8)$ | $\dfrac{16683417296961520597884426836631}{25000000000000000000000000000}$ | $667.33669187846082391537707346524$ | 667 |

Applying each certificate to $z_C$ and taking the integer part proves

$$
\boxed{A(18,6,7)\le309,\qquad A(19,6,7)\le460,\qquad A(19,6,8)\le667.}
$$

## 13. Relation to established methods

These constraints apply established techniques. Multiplying valid inequalities by Boolean monomials and linearizing uses the Sherali–Adams approach [2]. The matrices obtained from $hff^{\mathsf T}$ are localizing matrices; see [3], Section 4, especially equation (4.4). Finite-support occupancy constraints and binomial moments belong to the discrete moment and linear programming framework represented by [4]. Group averaging and representation block decomposition are standard symmetry methods [5], specialized to constant-weight codes in [1]. The anticodes arise from the classical intersecting families displayed in [6]; only their elementary pairwise-intersection property is needed here. No novelty is claimed for these underlying techniques.

## References

1. Sven Polak, “Semidefinite programming bounds for constant weight codes,” *IEEE Transactions on Information Theory* **65** (2019), 28–38. [DOI](https://doi.org/10.1109/TIT.2018.2854800); [Article and bibliographic record](https://arxiv.org/abs/1703.05171); [full text](https://arxiv.org/pdf/1703.05171). Section and equation numbers cited above refer to this arXiv version.

2. Hanif D. Sherali and Warren P. Adams, “A Hierarchy of Relaxations between the Continuous and Convex Hull Representations for Zero-One Programming Problems,” *SIAM Journal on Discrete Mathematics* **3**(3) (1990), 411–430. [DOI](https://doi.org/10.1137/0403036).

3. Jean B. Lasserre, “Global Optimization with Polynomials and the Problem of Moments,” *SIAM Journal on Optimization* **11**(3) (2001), 796–817. [DOI](https://doi.org/10.1137/S1052623400366802).

4. András Prékopa, “Boole-Bonferroni Inequalities and Linear Programming,” *Operations Research* **36**(1) (1988), 145–162. [DOI](https://doi.org/10.1287/opre.36.1.145).

5. Karin Gatermann and Pablo A. Parrilo, “Symmetry groups, semidefinite programs, and sums of squares,” *Journal of Pure and Applied Algebra* **192**(1–3) (2004), 95–128. [Article and bibliographic record](https://arxiv.org/abs/math/0211450); [full text](https://arxiv.org/pdf/math/0211450).

6. Rudolf Ahlswede and Levon H. Khachatrian, “The Complete Intersection Theorem for Systems of Finite Sets,” *European Journal of Combinatorics* **18**(2) (1997), 125–136. [Full text on the author's university site](https://www.math.uni-bielefeld.de/ahlswede/homepage/public/122.pdf).
