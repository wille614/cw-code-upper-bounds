# Label and coefficient specification for the public A(18,6,6) model

The [Chinese original](label_specification_zh.md) is the definitive text and takes precedence if the two versions differ.

This specification applies to the attachment `A18_6_6_model.json.gz`. The SHA-256 of the compressed file is `157bd486c563c49bdb65cec5b00e245e8f816b10ab67d721841c34b3c5aabe77`, and the SHA-256 of the decompressed JSON bytes is `ebf38af24dddea9230d2194faf394dc1666b75be5c40a574c6827ac2c8a0efc6`. All block numbers below are **zero-based indices** into the `blocks` array. In the comparison, first compute from the labels and the variable dictionary, then open the block's `entries` for comparison. The appendix gives 11 complete worked examples from the model producer. Reviewers can independently reconstruct the remaining blocks from the same definitions.

| Block numbers (inclusive) | `family` | Count | Definition |
|---|---|---:|---|
|0–64|D0|65|§3|
|65–268|D2|204|§4|
|269–1209|nonnegative|941|§5|
|1210–14677|anticode-vacancy-scalar|13468|§6|
|14678–14692|anticode-vacancy-matrix|15|§7|
|14693–14875|matching-occupancy|183|§8|
|14876–16248|coordinate-spectrum|1373|§9|
|16249–16300|occupancy-matrix|52|§10|
|16301|cardinality-cap|1|§11|
|16302–16341|cardinality-scalar|40|§12|
|16342–16363|cardinality-rooted|22|§13|
|16364–16386|leave-motif, leave-motif-anchor, leave-motif-defect-two|6＋13＋4|A public verification from the definitions is already available in `checks/leave_rows.py`, distributed with the model.|

## 1 From codewords to integer coefficients

The coordinate set is $I=\{0,\ldots,17\}$, the word set is $F=\binom I6$, and $V=\#F=18564$. The integer mask of a coordinate set is $\sum_{i\in S}2^i$, with coordinate 0 in the least significant bit. Two distinct words can belong to the same code if and only if their intersection has size at most 3. $X_u$ is the indicator that word $u$ belongs to the code, with $X_S=\prod_{u\in S}X_u$ and $X_\varnothing=1$. To evaluate a polynomial, first use $X_u^m=X_u$ ($m\ge1$), then set monomials containing a forbidden pair of words to zero and collect terms with the same support. The localizing matrices in this specification involve at most four distinct words.

**Orbits and numbering.** Give $k$ distinct words a temporary ordering $(u_0,\ldots,u_{k-1})$. The pattern of coordinate $i$ is $p=\sum_{j=0}^{k-1}2^{k-1-j}\mathbf1_{i\in u_j}$. Let $h_p$ be the number of coordinates of pattern $p$, ordered by $p=0,\ldots,2^k-1$. The lexicographically smallest histogram over all $k!$ permutations of the words gives $H(S)$. This completely characterizes the orbit of the **unordered** word set under coordinate permutations. To realize a given histogram, proceed in increasing order of $p$ and assign to pattern $p$ the next $h_p$ consecutive coordinates, starting at coordinate 0. Recover the individual words from the bits in each pattern. For a given histogram $h$, write
$$
 O(h)=\frac{18!}{\bigl(\prod_p h_p!\bigr)\,\#\{\pi\in S_k:\pi h=h\}} .
$$
This is the number of unordered word sets in the orbit. The last factor in the denominator is the number of permutations of the **words** that preserve the histogram.

To generate the variable table independently, enumerate the nonnegative integer histograms for $k=1,2,3,4$ with total 18, marginal sum 6 for each bit, and at most 3 positions with both bits equal to 1 for each pair of bits. Take the canonical form under word permutations as above, remove duplicates, and sort first by $k$ and then lexicographically by histogram, in ascending order. In the file, `variables[v−1]` is canonical form number $v$. The four layers have sizes 1, 4, 35, 901 and numbers 1, 2–5, 6–40, 41–941, respectively. $(12,6)$ has number 1. Given a collection of words, remove duplicates before determining its orbit.

**Linearization.** For an actual code $C$, define $z(S)=V\,\mathbb E_{g\in S_{18}}X_S(gC)$. Let $Z_0=1$, and let $Z_v$ be the variable for a nonempty orbit. The linear operator $\mathcal L$ replaces $X_S$ by $Z_{v(S)}$, a forbidden term by 0, and the constant 1 by $V Z_0$. Thus $z(\varnothing)=18564$ and $Z_1=|C|$. The coefficient of a fixed $X_S$ goes directly into its orbit. If the polynomial has already been summed over all coordinate configurations of the same type, multiply by $O(h)$ only once when counting. Apply $\mathcal L$ entry by entry to a matrix. A nonnegative $p$ gives $\mathcal L(p)\ge0$, and a nonnegative $h$ gives $\mathcal L(hff^{\mathsf T})\succeq0$. The global cardinality block in §11 has an additional explicit factor $1/V$.

**File conventions.** The top-level `n,d,w` are 18, 6, 6, respectively, `objective_variable` is 1, `box_bound` is 18564, and `cardinality_bound` is 186. `size` is the order of the matrix actually stored. A quadruple $[v,i,j,c]$ in `entries` gives the integer coefficient $c$ of $Z_v$ in upper-triangular entry $(i,j)$, with $0\le i\le j<\text{size}$. Each nonzero key occurs once, and all others are zero. Copy the lower triangle by symmetry. For coefficient comparisons the key is $(v,i,j)$, irrespective of the order in which the array is traversed. Multiply off-diagonal entries by 2 in the dual trace inner product. When constructing $(i,j)$, still expand $f_i f_j$ in full. In particular, cross terms in a square carry their multiplicities.

Where present, `positive_scale_removed` is the greatest common divisor $g>0$ of the absolute values of all nonzero integer coefficients of the entire block at the specified stage, including the affine column. First collect equal keys and remove zeros, then divide all coefficients by $g$, preserving their signs. This factor applies to the entire block, not to each row separately. D0 and D2 take $g$ **before** selecting the principal submatrix, with no further reduction after selection. Sections §7, §10 and §13 separately specify the order of zero-row removal or of the two stages. Blocks without this field are stored using the explicit formulas in their respective sections and must not all be reduced to lowest integer ratios. The coefficients 840 and −840 in §8 are the actual stored values.

## 2 Young vectors, projection and compression

The explicit matrix $Q$ is given below, along with a coefficient formula that avoids listing all of its rows. Both descriptions use the same normalization. Number the cells of shape $\lambda$ row by row from top to bottom and from left to right within each row. For a given group of actual coordinates, assign them to cells in the order specified for that group. A semistandard tableau has weakly increasing rows and strictly increasing columns. For a tableau $\tau$, consider the fillings $\rho$ obtained by taking **all distinct rearrangements of the symbols in each row**, counting each filling once. Let $K_\lambda$ be the direct product of the groups of arbitrary cell permutations within each column. For a string $x$ of letters on these coordinates, define the integer
$$
 u_\tau(x)=\sum_{\rho}\sum_{\pi\in K_\lambda}\operatorname{sgn}(\pi)
                 \mathbf1\{x=\pi\rho\}.
$$
Multiply the signs of the column permutations. Take the tensor product of the vectors for distinct coordinate groups, that is, multiply their $u$ values. Decode a string into one or more words, then form an unordered word set $T$ as specified for each family. Sum the weights for the same $T$ to obtain the projection column $P_{T,\tau}$. There is no division by row-length factorials, normalization to unit column vectors, or division by the number of orderings (two) of an unordered pair.

For D0 and D2, take the original matrix to be
$$
 M_{T,T'}(z)=z(T\cup T'),
$$
A forbidden combination has value zero, and the empty set has value $V$. The original row and column sets are given in §3 and §4. Form $P$ from the complete column sequence in the label, and first compute all coefficients of $P^{\mathsf T}MP$ and their common divisor $g$. Scan $P$ from left to right, keeping a column whenever it lies outside the rational span of the columns already selected. The resulting sequence of zero-based indices is $K$. The file does not list `kept_columns` separately. Let
$$
 Q=P_{:,K}/\sqrt{g},\qquad B(z)=Q^{\mathsf T}M(z)Q .
$$
The integer implementation uses $P_{:,K}^{\mathsf T}M P_{:,K}/g$. If the original unscaled probability matrix $M(x)$ is used, with $x=z/V$, then $Q$ in the same expression must also be multiplied by $\sqrt V$.

Column selection can be implemented on a smaller scale without forming $P$. All sets $T$ in the support of a given column have the same cardinality $s_i$. To compute $P^{\mathsf T}P$, include only $T=T'$ in the original expansion of $(i,j)$. Thus the entry is 0 if $s_i\ne s_j$. If $s_i=s_j=s>0$, sum the unreduced coefficients of that entry **before column selection** for which the union has cardinality exactly $s$. Two sets with the same cardinality as their union are equal. The additional empty-set column has inner product 1 with itself and 0 with every other column. Scanning from left to right, use exact rational elimination or positive principal minors to determine which columns increase the rank of this Gram matrix. This gives the same $K$. Across the stored D0 blocks, 592 of the original 1142 columns are retained. The corresponding numbers for D2 are 1479 original and 1475 retained.

**Column-order check.** D0 block 1 has shape $(17,1)$. Let $j$ be the letter in the bottom cell. On a string of fixed content, the value of the explicit column is $\mathbf1_{x_{17}=j}-\mathbf1_{x_0=j}$. When $t>0$, also add the same difference for the string obtained by exchanging the two words (interchanging letters 1 and 2). Selecting linearly independent columns within each content retains the original indices $0,1,2,4,5,7,8,10$ among the 12 columns. In D2 only blocks 65, 104, 180 and 215 lose columns. In all four blocks $h$ is zero, and a completing word equal to either anchor word gives the same unit column for $D$. The later of the two columns has original index 1, 2, 2, 2, respectively. Delete that column in each block.

**Direct computation of the generating polynomial of a pair of tableaux.** For an alphabet of size $m$ and a matrix of indeterminates $Y=(Y_{ab})_{0\le a,b<m}$, the exponents of each monomial in $\sum_{x,y}u_\tau(x)u_\sigma(y)\prod_\ell Y_{x_\ell y_\ell}$ give the joint histogram of the letters at corresponding coordinates. The following expression computes this sum using integers. Pad $\lambda$ with zeros, numbering the row lengths from 1 in this expression,
$$
 H_\lambda(Y)=\prod_{l=1}^m\bigl(l!\det Y_{[0,l),[0,l)}\bigr)^{\lambda_l-\lambda_{l+1}} .
$$
Start from $H_\lambda$. First, for the left tableau $\tau$, take its rows $j$ from bottom to top (row numbers start at 0). Within each row, take the letters in the order $i=m-1,\ldots,j+1$ and let $e$ be the number of occurrences of letter $i$ in that row. Apply in succession
$$
 \frac1{e!}\left(\sum_b Y_{ib}\frac{\partial}{\partial Y_{jb}}\right)^e .
$$
Then, for the right tableau $\sigma$, apply the column operators $(\sum_a Y_{ai}\partial/\partial Y_{aj})^e/e!$ in the same order. An empty coordinate group contributes 1. When multiplying over several coordinate groups, use a distinct $Y$ for each group. Recover the words from the exponents and combine terms in the same orbit as in §1. The monomials already carry the permutation multiplicities, so do not multiply by $O(h)$ again. An implementation of this formula can also be checked against the explicit $u_\tau$ and vice versa. (For the four-point matrices and representative vectors, see §§3.1–3.2 and the tableau polynomials in the appendix.) [Polak's four-point bound](https://arxiv.org/abs/1703.05171)

## 3 D0 — empty anchor and four-letter alphabet

`family=D0`. `partition` is a partition of 18 into at most four rows, with positive row lengths in nonincreasing order. `tableaux` is the ordered list of tableaux of this shape, with letters $0,1,2,3$ decoded respectively as the bits $00,01,10,11$ of an ordered pair of words at each coordinate. Each tableau has content $(12-t,t,t,6-t)$, with $t\in\{0,3,4,5,6\}$. For a fixed shape, order first by $t$ as above, then within each content by the lexicographic order of the row-by-row reading. `positive_scale_removed` is defined in §2 and ranges from 64 to 860160 in this file.

The original rows are the empty set, singleton word sets and admissible two-word sets. If $\lambda=(18)$, insert a unit column for the empty set before all tableaux. For each remaining column, decode $x$ into two words $(b,c)$ and add $u_\tau(x)$ to original row $\{b,c\}$. When $b=c$, use the singleton row. The vector components can be written as $f_i=\sum_T P_{Ti}X_T$. The component for the empty-set column is 1. The stored block is $B(z)$ of §2. For example, in the single-row shape, before reduction the empty column paired with itself gives $18564 Z_0$, and its cross term with a tableau of content parameter $t$ is $18564\binom6t\binom{12}t z(\{b,c\})$, where $|b\setminus c|=t$. Use the same reduction and column-selection rules for the single-row shape and all other shapes.

**Block 64.** $\lambda=(6,6,6)$, and the three rows of the unique tableau consist of six 0s, six 1s and six 2s, respectively. The tableau-pair polynomial is $(6\det Y_{[0,3),[0,3)})^6$. Choose $\pi\in S_3$ six times, each time increasing the exponent at $(i,\pi(i))$ by 1 and multiplying the weight by $6\operatorname{sgn}(\pi)$. Decode each exponent array into four words using the four bits $(b,c,b',c')$, collapse repeated words, discard forbidden combinations and collect terms. The common divisor is 93312, and the sole column is retained. There are 37 terms, including $Z_2+2Z_6-12Z_7$. All 37 terms are given in the appendix.

## 4 D2 — two-word anchor and four groups with two-letter alphabets

`family=D2`. `anchor` stores the **four cell sizes** $c=(12-t,t,t,6-t)$, with $t=3,4,5,6$. Partition the coordinates into consecutive groups in order. The anchor words $(a,b)$ have bits $00,01,10,11$ in the respective groups. The four entries $h_p$ of `harmonic` satisfy $0\le h_p\le\lfloor c_p/2\rfloor$. `rows` is the complete list of weight vectors $\alpha=(\alpha_p)$, in lexicographic order. Require $\sum\alpha_p=6$ and $h_p\le\alpha_p\le c_p-h_p$. The distances from the new word $x$ to $a,b$ are $2(6-\alpha_2-\alpha_3)$ and $2(6-\alpha_1-\alpha_3)$, respectively, each either zero or at least 6. `positive_scale_removed` ranges from 1 to 256 in this file.

The original rows are $D=\{a,b\}$ and the admissible sets $D\cup\{x\}$. Group $p$ uses a two-letter tableau of shape $(c_p-h_p,h_p)$ with top row $0^{c_p-\alpha_p}1^{\alpha_p-h_p}$ and bottom row $1^{h_p}$. Omit empty rows. An empty group has only the empty string, of weight 1. Accumulate the product vector in row $T=D\cup\{x\}$ as in §2. If two distinct strings give $x=a$ or $x=b$, both go into the same row $D$. The component is $\sum_T P_{T,\alpha}X_T$, with $Q$ as in §2. In the generating polynomial, the joint bits in each group are in the order $00,01,10,11$ and represent the two new words. Concatenate the exponents from the four groups in order to form the histogram of the four bits $(a,b,x,y)$.

**Block 230.** $c=(6,6,6,0)$ and $h=(0,3,3,0)$, with the unique $\alpha=(0,3,3,0)$. The first group has six occurrences of $00$, the last group is empty, and the two middle groups each contribute $(2\det Y)^3$. For $p,q=0,\ldots,3$, the four exponents of the two groups are $(3-p,p,p,3-p)$ and $(3-q,q,q,3-q)$, respectively, with joint weight $64(-1)^{p+q}\binom3p\binom3q$. Collect the 16 terms by four-word orbit, divide by 64 and retain the sole column. The result is $Z_{27}+Z_{863}-6Z_{864}+6Z_{865}-2Z_{866}+9Z_{867}-18Z_{868}$.

## 5 nonnegative — nonempty orbits

The only label is `family=nonnegative`. From first to last, the blocks in this family represent $Z_1,\ldots,Z_{941}\ge0$, respectively. If $b$ is the block number in the whole file, then $v=b-268$. The position determines the object. Copying the label alone loses the index. Each block contains only $[v,0,0,1]$, with scale 1.

**Block 269.** This is the first block in the family, with histogram $(12,6)$ and $\mathcal L(X_{\{u\}})=Z_1$. The full list of coefficients is $[1,0,0,1]$.


## 6 anticode-vacancy-scalar — three-word anchor and anticode vacancy

`anchor` is the eight-entry canonical histogram of an admissible three-word set. Use the consecutive coordinate segments of §1 to recover the three-word set $D$ **in the same coordinate system**. `subset` is an 18-bit mask specifying $R\subset I$. Its numbering depends on the representative above. `r` takes the values 0, 1, 2, with $|R|=4+2r$ and `threshold` $=4+r$. Let
$$
 \mathcal A_R=\{u\in F: |u\cap R|\ge4+r\},\qquad H_R=1-\sum_{u\in\mathcal A_R}X_u .
$$
Any two words in this set share at least 4 points in $R$, so $H_R\ge0$. Store $\mathcal L(X_DH_R)$. For the coefficients, add 1 at $v(D)$, then subtract 1 at the corresponding index for every admissible $D\cup\{u\}$, with scale 1. A label corresponds to one fixed $R$, without summation over $R$. The model selects anticodes containing no anchor word. Within each pattern segment of the fixed anchor histogram, the initial coordinates represent the selected points. Only one copy of a constraint with the same coefficients is needed.

**Block 1215.** The histogram $(0,6,6,0,6,0,0,0)$ of $D$ gives three disjoint six-point segments. Here $r=1$, the mask of $R$ is 61505, and the threshold is 5. $\#\mathcal A_R=\binom65\binom{12}1+1=73$. Take $u$ in turn to be each of these words and retain those meeting each six-point segment in at most 3 points. Exactly 40 remain, and all yield orbit 42. The anchor orbit is 6. The whole block is $Z_6-40Z_{42}$.

## 7 anticode-vacancy-matrix — one-word anchor and anticode vacancy matrix

`r,threshold,subset` are defined as in §6. **In this family**, `anchor` is the mask of a six-point word $a$. We have `overlap` $=|a\cap R|=s$ and $0\le s\le\min(6,|R|,3+r)$, with nonnegative cell sizes. `caps` gives the four group sizes $(18-6-|R|+s,\,|R|-s,\,6-s,\,s)$ in the order $(\mathbf1_a,\mathbf1_R)=00,01,10,11$. All 15 blocks in this file have `harmonic` equal to $[0,0,0,0]$. `rows` lists the retained four-group weight vectors in lexicographic order. The value of `positive_scale_removed` is 1 throughout.

Start with $\sum\alpha_p=6$ and $0\le\alpha_p\le c_p$. Require $\alpha_1+\alpha_3<4+r$ and require $2(6-\alpha_2-\alpha_3)$ to be zero or at least 6. The set $F_\alpha$ consists of all words whose weights in the four groups are $\alpha$, and $f_\alpha=\sum_{b\in F_\alpha}X_b$ without division by the class size. The component for $F_\alpha=\{a\}$ is $X_a$ and is equivalent to a constant component after multiplication by the anchor indicator. The block is
$$
 \mathcal L(X_aH_R f f^{\mathsf T}).
$$
For $(\alpha,\beta)$, iterate over the ordered choices $b\in F_\alpha,c\in F_\beta$ and accumulate $z(\{a,b,c\})-\sum_{q\in\mathcal A_R}z(\{a,b,c,q\})$ as in §1. Compute the upper triangle for all candidate rows, take the common divisor of the entire block, remove every row that is zero for all variables together with its symmetric column, and renumber the rest in the original order. The label field `rows` already gives the retained sequence.

**Block 14678.** Here $a=258048,\ R=3840,\ r=s=0$, the threshold is 4, and the cell sizes are $(8,4,6,0)$. The 17 components are in the lexicographic order given by the preceding constraints (the appendix lists the actual label). For example, entry $(0,0)$ corresponds to $b=c=a$. The anticode words here consist of $R$ plus two points. The numbers with 0, 1 and 2 of those points in $a$ are 28, 48 and 15, respectively, so the entry is $Z_1-28Z_2-48Z_3-15Z_4$.

All entries can be recomputed by the following smaller double count. For each nonempty index $v$, fix a representative $S$. In $S$, make ordered choices of $a,b,c$ (repetitions allowed), adding a positive term only when $\{a,b,c\}=S$. Enumerate all four-point sets $R$ in $I\setminus a$ and locate $(i,j)$ from the cell weights of $b,c$. For the negative term, make ordered choices $a,b,c,q\in S$ with union support $\{a,b,c,q\}=S$. Here enumerate only the four-point sets $R$ in $q\setminus a$. In both counts, include only choices for which the rows of $b,c$ exist and $i\le j$. Subtract the second count from the first and multiply each entry by $O(H(S))/(18564\binom{12}4)=O(H(S))/9189180$. This recovers the coefficient for fixed $(a,R)$ from the sum over the entire orbit and all $(a,R)$. The resulting 10215 terms have common divisor 1 and are listed individually in the appendix.

## 8 matching-occupancy — an exact number of words in a matching interval

`anchor` is the mask of word $a$ and is fixed at 63. `lower,upper` are the masks of the lower endpoint $L$ and the upper endpoint $U$, respectively, with $L\subseteq U$. The three entries of `caps` are ordered as $(18-|U|,\,|U\setminus L|,\,|L|)$ and take the values $(6,9,3),(5,10,3),(4,11,3),(9,9,0)$ in the file. Let $\mathcal I=\{u\in F:L\subseteq u\subseteq U\}$, $N=\sum_{\mathcal I}X_u$ and $m=|U\setminus L|$. The remaining word weight $k=6-|L|$ satisfies $\min(k,m-k)=3$. By the distance constraint, the remaining three-point sets (taking complements within $U\setminus L$ if necessary) are pairwise disjoint, so `capacity` $=\lfloor m/3\rfloor=3$. The field `occupied` is $j=0,1,2$.

The object is $\mathcal L(X_a\,\mathbf1_{\{N=j\}})$. The exact polynomial of degree at most three is
$$
 \mathbf1_{\{N=j\}}=\sum_{k=j}^{3}(-1)^{k-j}\binom kj\binom Nk .
$$
In $\binom Nk$ each unordered set of $k$ **distinct** interval words is counted once. Then take the set union with $a$, discard forbidden combinations and assign the terms to their orbits. The formula still applies when $a\in\mathcal I$. There is no additional reduction.

**Block 14695.** Here $a=63,L=448,U=262080$, the three group sizes are $(6,9,3)$, and $j=2$. The 84 interval words consist of $L$ plus three of the nine middle points and are all compatible with $a$. There are $\binom93\binom63/2=840$ compatible pairs, which belong to orbit 11 together with the anchor. There are $9!/((3!)^3\,3!)=280$ partitions into three triples, which together with the anchor belong to orbit 88. Substitution in $\binom N2-3\binom N3$ gives $840Z_{11}-840Z_{88}$.

## 9 coordinate-spectrum — occupancy polynomials over all coordinate intervals of the same size

`inside` $=\ell$ and `outside` $=b$ are the numbers of points that must be included and avoided, ranging in the file over 0–3 and 0–9, respectively. Iterate over all ordered pairs $(R,Z)$ of disjoint sets with $|R|=\ell,|Z|=b$, counting each pair once. A local word must contain $R$ and avoid $Z$. Their count is $N_{R,Z}$. Let $m=18-\ell-b,k=6-\ell$. The value of `capacity` is the integer $J(m,k)$ defined as follows. Take 0 for a weight out of range, take 1 if $0\le k\le m$ and $\min(k,m-k)<3$, and otherwise take
$$
 J(m,k)=\min\left(\left\lfloor\frac{mJ(m-1,k-1)}{k}\right\rfloor,\ 
                         \left\lfloor\frac{mJ(m-1,k)}{m-k}\right\rfloor\right).
$$
This family uses $2\le J\le20$. The field `roots` lists 1 to $\min(4,J)$ distinct integers from $0,\ldots,J$ in increasing order. The field `sign` is $+1$ or $-1$ as specified in the label. The polynomial is $p(x)=\text{sign}\prod_{r\in\text{roots}}(x-r)$ and is nonnegative at the integers $0,\ldots,J$. The object is $\mathcal L(\sum_{R,Z}p(N_{R,Z}))$ and has no vector components.

Let $\Delta^k p(0)=\sum_{j=0}^k(-1)^{k-j}\binom kj p(j)$. Using $p(N)=\sum_k\Delta^kp(0)\binom Nk$ gives the unreduced constant coefficient $p(0)\binom{18}{\ell}\binom{18-\ell}b V$. For orbit $v$ with $k$-word representative histogram $h$, the corresponding coefficient is
$$
 \Delta^kp(0)\ O(h)\binom{h_{2^k-1}}\ell\binom{h_0}b .
$$
A binomial coefficient outside its range is 0. Finally, divide all coefficients by `positive_scale_removed` (observed minimum 55692 and maximum 137225088000 for this family). This sum ranges over all intervals and must not be restricted to one fixed interval.

**Block 14876.** Here $\ell=0,b=7,J=13,\text{roots}=[13],\text{sign}=-1$, so $p(N)=13-N$. There are 31824 intervals, the original constant is $13\cdot31824\cdot18564=7680149568$, and the singleton coefficient is $-18564\binom{12}7=-14702688$. The common divisor is 1336608, and the whole block is $5746Z_0-11Z_1$.

## 10 occupancy-matrix — occupancy matrix of a coordinate core

`core` is the core size $s$. The field `complement` is Boolean. When false, use $y(u)=u$ and oriented weight $w'=6$. When true, use $y(u)=I\setminus u$ and $w'=12$. Both cases use the original indicator variable $X_u$. Iterate over all $s$-point sets $T$. Let $N_T=\sum_{u:T\subseteq y(u)}X_u$. The value `capacity` $=J(18-s,w'-s)$ is given by §9. The triples $(s,\text{complement},J)$ used in the file are exactly $(2,\mathrm{false},20),(3,\mathrm{false},5),(7,\mathrm{true},13),(8,\mathrm{true},6),(9,\mathrm{true},3)$.

`classes` lists $\max(0,w'-(18-s)),\ldots,\min(s,w')$ in ascending order. **Component 0 is 1.** The remaining components are $f_{j,T}=\sum_{u:|y(u)\cap T|=j}X_u$ in the order of `classes`, including the class $j=0$ if present. Here $j$ is the intersection size. The stored index is its zero-based position in `classes` plus 1. The fields `kind` and `polynomial` together specify $p(x)=a+bx+cx^2$, with the three entries in order from constant to quadratic. The label `["gap",j]` corresponds to $(x-j)(x-j-1)$, with $0\le j<J$, and `["capacity",J]` corresponds to $x(J-x)$. The block is
$$
 \mathcal L\left(\sum_{|T|=s}p(N_T)\,f_T f_T^{\mathsf T}\right).
$$
First expand using the full set of components and take the common divisor `positive_scale_removed` (55692 to 81681600), then take all nonzero rows. The field `kept_rows` lists original indices in increasing order, with both rows and columns remapped to $0,\ldots,\text{size}-1$. No additional core weight vector, orbit averaging or division by a class size is needed.

The same orbit-by-orbit formula works for both orientations. For every subset $B$ of $S$, set $X_u$ to $\mathbf1_{u\in B}$ to obtain the numerical value $G_{ij}(B,T)$ in the expression above. The unreduced coefficient of variable $v(S)>0$ is
$$
 O(H(S))\sum_{|T|=s}\sum_{B\subseteq S}(-1)^{|S|-|B|}G_{ij}(B,T).
$$
This is the coefficient supported on exactly $S$ in the Boolean polynomial, counted over its orbit. For the empty set, store $V\sum_T G_{ij}(\varnothing,T)$ under key 0 instead.

**Block 16249.** Here $s=2$ in the original-word orientation, $J=20$ and $p(x)=x^2-x$. With $\text{classes}=[0,1,2]$, the four components are 1 and the counts of words whose intersections have sizes 0, 1 and 2, respectively. Use the formula above to sum over the 153 two-point sets and all subsets of the fixed representative, then divide all entries by 15315300. All four rows are retained. For example, only a pair of distinct words contributes to entry $(0,0)$, and an unordered pair of intersection size $i$ contributes $O(h)\,2\binom i2$. Thus the entry is $9Z_4+16Z_5$. The resulting 2238 terms are given in the appendix.

## 11 cardinality-cap — global code size

The sole value of `bound` is $B=186$. The known bound $A(17,6,6)\le124$ gives $B=\lfloor18\cdot124/12\rfloor$ by shortening on a zero coordinate, independently of the bound 183 to be proved. Let $M=\sum_{u\in F}X_u$. This block stores $\mathcal L(B-M)/V$ because $\mathcal L(M)=VZ_1$.

**Block 16301.** Divide both the constant term $186\cdot18564$ and $-18564Z_1$ by 18564 to get $186Z_0-Z_1$, with only two terms. The label does not need a separate scale field.

## 12 cardinality-scalar — the global bound multiplied by a nonempty anchor monomial

`bound` $=186$. The field `anchor` gives a histogram for each of the 40 admissible orbits on 1, 2 and 3 words, taken once each. Their lengths are 2, 4 and 8, respectively. Recover the fixed word set $D$ as in §1. Store $\mathcal L(X_D(B-M))$ without dividing by $V$ or any other factor. Add $B$ at $v(D)$, then iterate over all 18564 words $u$ and subtract 1 from the variable for each admissible $D\cup\{u\}$. All $|D|$ terms with $u\in D$ are subtracted from the anchor orbit itself.

**Block 16307.** Here $D$ consists of the three disjoint six-point segments in block 1215 of §6. The three terms with $u\in D$ change the coefficient of orbit 6 to $186-3=183$. The numbers of points chosen by the remaining compatible words in the three segments are permutations of $(3,3,0),(3,2,1),(2,2,2)$, corresponding to orbits 41, 42 and 43, respectively. Their counts are $3\binom63^2=1200$, $6\binom63\binom62\binom61=10800$ and $\binom62^3=3375$, respectively, and the whole block is $183Z_6-1200Z_{41}-10800Z_{42}-3375Z_{43}$.

## 13 cardinality-rooted — cardinality localization of the three-point matrix with a one-word anchor

The outer field is `bound` $=186$. The inner field `source` gives the three-point matrix for a fixed anchor, with the fields below.

| Field of `source` | Mathematical meaning and values in this attachment |
|---|---|
|`family,anchors,k,row_arity`|The values are `D1`, `[63]`, 3 and 1, respectively. The anchor set contains a single word $a$, the union for the source matrix contains at most three words, and each original row allows one completing word.|
|`cells`|Two ordered groups of coordinates, first 6 through 17 and then 0 through 5, giving $I\setminus a,a$. Within each group, align the coordinates with the cells of the Young diagram in label order.|
|`shapes`|One partition for each group, of sizes 12 and 6, respectively, each with at most two rows.|
|`columns`|Take the projection columns in the array order in the label. Each column contains two semistandard two-letter tableaux with shapes given by `shapes` and a total of 6 occurrences of 1. The decoded word $x$ must satisfy $x=a$ or be compatible with $a$.|
|`empty_column`|Fixed at `false`, with no additional unit column for the empty set.|
|`row_cardinalities`|Gives the cardinality of original row $\{a,x\}$ for each column. It is 1 when $x=a$ and 2 otherwise.|
|`positive_scale_removed`|The positive integer $g_s$ by which all linearized entries of the source block are divided before cardinality localization. The values used are 1, 8, 48, 64, 72, 80, 120, 128, 220, 240, 256, 264, 720, 1760.|

For the two tableaux in column $i$, use the product in §2 to obtain $u_i(x)$, and let $f_i=\sum_x u_i(x)X_x$. First obtain $g_s$ from all the integer coefficients of $\mathcal L(X_a f_i f_j)$ and divide every source coefficient by $g_s$. For each $Z_{v(S)}$ in the source block, apply the same linear substitution
$$
 T_B(Z_{v(S)})=B Z_{v(S)}-\sum_{u\in F} z(S\cup\{u\})
$$
with forbidden terms set to zero. Since $M$ is invariant under coordinate permutations, any representative of the orbit of $S$ may be used. The resulting matrix is $\mathcal L((B-M)X_af f^{\mathsf T})/g_s$. Divide the whole result by the outer `positive_scale_removed` (1 or 2 in this file), then delete all-zero rows and columns. The outer `kept_rows` lists the original zero-based source-column indices in increasing order. Restrict to these and renumber. This family only deletes zero rows. Do not additionally apply the rank-increase compression of §2.

**Block 16363.** The outer group has shape $(6,6)$ with six 0s in the top row and six 1s in the bottom row. The inner group has a single row of six 0s. Pair the outer coordinates as $(6,12),(7,13),\ldots,(11,17)$. Each of the 64 words chooses one point from each pair. Its sign is determined by the parity of the number of times a 1 is swapped into the top cell. All entries in the inner group are 0. Enumerate $64^2=4096$ ordered word pairs and accumulate $z(\{63,x,y\})$ with the product of their signs (63 here denotes a word by its mask). Divide the entire source block by 64 to obtain
$$
 Z_2+Z_6-6Z_7+15Z_8-20Z_{11}.
$$
For each of these five terms, apply $T_{186}$ by iterating over all 18564 words and sum by orbit. The outer common divisor is 1, and only row 0 is retained. There are 120 terms in the final result, including $184Z_2+182Z_6-1170Z_7$. The full table is in the appendix.

## Appendix. Complete upper triangles for 11 worked examples

In each family take the block of smallest order, breaking ties by the fewest nonzero terms and then by the earliest position in the file. The coefficients below are obtained by the computations in the respective sections. Each row gives one upper-triangular position, with continuation lines for positions with many terms. Each pair `v:c` denotes the coefficient $cZ_v$. All unlisted keys are zero.

For each block, the reported SHA-256 is calculated by sorting the quadruples $[v,i,j,c]$ lexicographically as integer quadruples and writing them in UTF-8 as a compact JSON array with no spaces or newlines. The comparison separately requires exact equality of the full lists of quadruples, including the nonzero counts and coordinates.

### Block 64 — D0

Order 1, 37 nonzero entries. Comparison result **PASS**. The coefficient SHA-256 is `9fa1da79fc901cba7779af2376f0e9f11c693dcc948572d6f9bf8ee1479208a9`.

Actual label

```json
{"family":"D0","partition":[6,6,6],"tableaux":[[[0,0,0,0,0,0],[1,1,1,1,1,1],[2,2,2,2,2,2]]],"positive_scale_removed":93312}
```

| Row $i$ | Column $j$ | Variable indices and coefficients |
|---:|---:|---|
|0|0|`2:1` `6:2` `7:-12` `8:30` `11:-40` `41:-40` `46:-240` `49:240` `60:180` `64:240` `72:-240` `77:-360`|
|0|0|`91:20` `94:180` `113:240` `123:240` `133:-720` `138:240` `165:360` `174:360` `182:-600` `225:-240` `236:-360` `257:360`|
|0|0|`263:-600` `290:720` `324:180` `326:-270` `379:40` `402:720` `458:-600` `472:360` `575:-120` `621:-360` `655:60` `763:120`|
|0|0|`863:-10`|

### Block 230 — D2

Order 1, 7 nonzero entries. Comparison result **PASS**. The coefficient SHA-256 is `216b89e3d31736771bca4a0b24eb7b4d89e30df4739f2e6e969e79c80a3b24f0`.

Actual label

```json
{"family":"D2","anchor":[6,6,6,0],"harmonic":[0,3,3,0],"rows":[[0,3,3,0]],"positive_scale_removed":64}
```

| Row $i$ | Column $j$ | Variable indices and coefficients |
|---:|---:|---|
|0|0|`27:1` `863:1` `864:-6` `865:6` `866:-2` `867:9` `868:-18`|

### Block 269 — nonnegative

Order 1, 1 nonzero entries. Comparison result **PASS**. The coefficient SHA-256 is `745540e3cf79e3a0170d63806f5c9290aa72dd0c76e7c1aff03c5775f9489463`.

Actual label

```json
{"family":"nonnegative"}
```

| Row $i$ | Column $j$ | Variable indices and coefficients |
|---:|---:|---|
|0|0|`1:1`|

### Block 1215 — anticode-vacancy-scalar

Order 1, 2 nonzero entries. Comparison result **PASS**. The coefficient SHA-256 is `8b3073237a43f01439b5bcc1e2e80af11921125466dae92be5fc89dff7edbdc2`.

Actual label

```json
{"family":"anticode-vacancy-scalar","anchor":[0,6,6,0,6,0,0,0],"r":1,"subset":61505,"threshold":5}
```

| Row $i$ | Column $j$ | Variable indices and coefficients |
|---:|---:|---|
|0|0|`6:1` `42:-40`|

### Block 14678 — anticode-vacancy-matrix

Order 17, 10215 nonzero entries. Comparison result **PASS**. The coefficient SHA-256 is `3a7487ca625893e76bc65faf6c97273c6498d3212c58b76da9c34f98f8b258d6`.

Actual label

```json
{"family":"anticode-vacancy-matrix","r":0,"overlap":0,"caps":[8,4,6,0],"harmonic":[0,0,0,0],"rows":[[0,0,6,0],[0,3,3,0],[1,2,3,0],[1,3,2,0],[2,1,3,0],[2,2,2,0],[2,3,1,0],[3,0,3,0],[3,1,2,0],[3,2,1,0],[3,3,0,0],[4,0,2,0],[4,1,1,0],[4,2,0,0],[5,0,1,0],[5,1,0,0],[6,0,0,0]],"anchor":258048,"subset":3840,"threshold":4,"positive_scale_removed":1}
```

| Row $i$ | Column $j$ | Variable indices and coefficients |
|---:|---:|---|
|0|0|`1:1` `2:-28` `3:-48` `4:-15`|
|0|1|`5:80` `27:-2240` `34:-1920` `38:-240`|
|0|2|`5:960` `21:-20160` `27:-6720` `28:-20160` `29:-20160` `34:-2880` `35:-2880` `36:-8640`|
|0|3|`4:480` `21:-10080` `28:-13440` `35:-2880`|
|0|4|`5:2240` `15:-33600` `21:-26880` `22:-40320` `24:-40320` `27:-2240` `28:-20160` `29:-13440` `31:-20160` `32:-6720`|
|0|5|`4:2520` `16:-37800` `21:-30240` `23:-60480` `24:-30240` `28:-20160` `30:-15120` `31:-20160`|
|0|6|`3:672` `15:-10080` `22:-20160` `28:-6720`|
|0|7|`5:1120` `11:-11200` `15:-33600` `18:-16800` `21:-6720` `22:-10080` `24:-20160` `26:-3360`|
|0|8|`4:3360` `12:-33600` `16:-50400` `17:-67200` `19:-33600` `21:-10080` `23:-60480` `24:-20160` `25:-26880` `26:-3360`|
|0|9|`3:2016` `12:-20160` `15:-30240` `17:-50400` `18:-10080` `22:-30240` `23:-20160` `24:-10080`|
|0|10|`2:224` `11:-2240` `15:-6720` `21:-3360`|
|0|11|`4:1050` `8:-6300` `12:-33600` `14:-8400` `16:-12600` `17:-16800` `19:-16800` `20:-1050`|
|0|12|`3:1680` `9:-10080` `12:-26880` `13:-33600` `14:-6720` `15:-10080` `17:-50400` `18:-6720` `19:-8400`|
|0|13|`2:420` `8:-2520` `11:-6720` `12:-10080` `15:-10080` `16:-6300`|
|0|14|`3:336` `7:-1008` `9:-10080` `10:-1008` `12:-6720` `13:-8400` `14:-3360`|
|0|15|`2:224` `7:-672` `8:-3360` `9:-4032` `11:-2240` `12:-10080`|
|0|16|`2:28` `6:-28` `7:-672` `8:-840` `9:-1008`|
|1|1|`5:80` `27:-2240` `34:-1920` `39:2160` `40:80` `865:-6720` `866:-2240` `868:-60480` `928:-17280`|
|1|2|`34:1920` `36:17280` `37:17280` `38:1920` `39:17280` `764:-40320` `765:-40320` `775:-362880` `776:-362880` `779:-362880` `864:-13440` `865:-13440`|
|1|2|`867:-120960` `868:-241920` `889:-120960` `890:-120960` `900:-40320` `902:-40320` `905:-241920` `906:-241920` `907:-362880` `927:-17280` `928:-51840` `933:-17280`|
|1|2|`934:-34560` `938:-69120`|
|1|3|`35:5760` `36:17280` `38:1920` `777:-120960` `778:-40320` `779:-362880` `899:-40320` `901:-13440` `905:-241920` `933:-17280`|
|1|4|`27:2240` `29:20160` `32:20160` `33:2240` `34:6720` `36:60480` `37:60480` `575:-33600` `576:-100800` `605:-302400` `606:-907200` `632:-302400`|
|1|4|`633:-907200` `636:-33600` `752:-120960` `760:-362880` `763:-26880` `764:-120960` `772:-120960` `774:-241920` `775:-725760` `777:-241920` `778:-26880` `779:-725760`|
|1|4|`796:-241920` `799:-725760` `813:-241920` `817:-725760` `822:-40320` `823:-120960` `825:-362880` `863:-2240` `864:-40320` `865:-26880` `866:-2240` `867:-60480`|
|1|4|`868:-60480` `873:-40320` `874:-20160` `875:-40320` `882:-60480` `889:-120960` `897:-120960` `899:-80640` `900:-120960` `901:-20160` `902:-80640` `904:-40320`|
|1|4|`905:-241920` `906:-241920` `907:-120960` `916:-120960` `919:-60480`|
|1|5|`28:20160` `31:60480` `32:20160` `35:20160` `36:60480` `624:-302400` `625:-302400` `634:-907200` `635:-1209600` `768:-120960` `771:-120960` `774:-241920`|
|1|5|`775:-725760` `776:-241920` `777:-241920` `779:-725760` `807:-725760` `813:-241920` `816:-1088640` `817:-362880` `823:-241920` `825:-362880` `873:-40320` `875:-40320`|
|1|5|`886:-60480` `889:-241920` `890:-120960` `896:-120960` `897:-120960` `899:-40320` `904:-40320` `905:-241920` `916:-120960`|
|1|6|`28:20160` `29:20160` `34:6720` `632:-302400` `633:-302400` `636:-100800` `796:-241920` `799:-362880` `822:-80640` `874:-20160` `882:-60480` `901:-6720`|
|1|7|`27:4480` `29:40320` `32:40320` `33:4480` `379:-44800` `392:-403200` `415:-403200` `434:-44800` `568:-201600` `575:-67200` `576:-67200` `593:-403200`|
|1|7|`605:-604800` `620:-403200` `632:-604800` `636:-134400` `648:-201600` `742:-40320` `752:-120960` `754:-80640` `763:-13440` `764:-40320` `765:-13440` `774:-120960`|
|1|7|`777:-120960` `778:-26880` `796:-241920` `798:-80640` `813:-241920` `815:-40320` `822:-40320` `823:-120960`|
|1|8|`21:13440` `24:40320` `26:13440` `28:40320` `31:120960` `32:40320` `402:-134400` `404:-403200` `425:-403200` `428:-1209600` `432:-134400` `435:-403200`|
|1|8|`571:-67200` `579:-201600` `602:-403200` `617:-1209600` `621:-201600` `622:-134400` `624:-1209600` `625:-201600` `628:-201600` `630:-604800` `634:-1814400` `635:-604800`|
|1|8|`646:-604800` `650:-604800` `735:-40320` `736:-26880` `747:-40320` `758:-120960` `763:-40320` `764:-201600` `765:-53760` `768:-362880` `769:-40320` `770:-161280`|
|1|8|`771:-120960` `772:-120960` `774:-120960` `775:-362880` `776:-120960` `794:-120960` `807:-725760` `811:-241920` `813:-241920` `815:-40320` `816:-362880` `817:-362880`|
|1|9|`22:40320` `24:40320` `28:40320` `29:40320` `425:-403200` `428:-403200` `432:-403200` `435:-403200` `602:-403200` `605:-604800` `606:-604800` `617:-604800`|
|1|9|`620:-201600` `628:-403200` `632:-604800` `633:-604800` `646:-604800` `648:-201600` `747:-40320` `752:-241920` `754:-80640` `758:-120960` `760:-362880` `769:-40320`|
|1|9|`794:-120960` `796:-241920` `798:-80640` `799:-362880`|
|1|10|`21:13440` `27:4480` `415:-134400` `434:-44800` `593:-201600` `636:-67200` `742:-40320` `778:-13440`|
|1|11|`21:16800` `24:50400` `26:16800` `225:-100800` `235:-302400` `246:-100800` `377:-67200` `391:-403200` `402:-268800` `403:-134400` `414:-201600` `425:-806400`|
|1|11|`427:-201600` `432:-268800` `565:-50400` `571:-67200` `572:-33600` `586:-50400` `602:-403200` `604:-100800` `621:-100800` `622:-134400` `623:-16800` `624:-302400`|
|1|11|`625:-100800` `628:-201600` `630:-201600`|
|1|12|`15:16800` `18:16800` `22:50400` `24:50400` `236:-100800` `240:-403200` `248:-302400` `395:-134400` `400:-604800` `402:-268800` `403:-67200` `404:-268800`|
|1|12|`423:-604800` `425:-806400` `427:-201600` `428:-806400` `561:-16800` `571:-134400` `572:-33600` `574:-100800` `575:-100800` `576:-168000` `579:-201600` `599:-151200`|
|1|12|`602:-403200` `604:-100800` `605:-302400` `606:-302400` `617:-604800` `620:-201600`|
|1|13|`15:16800` `21:16800` `235:-100800` `246:-100800` `391:-201600` `392:-268800` `414:-201600` `415:-268800` `565:-50400` `568:-201600` `586:-50400` `593:-201600`|
|1|14|`15:13440` `18:13440` `113:-40320` `118:-40320` `222:-80640` `229:-120960` `236:-201600` `237:-40320` `240:-201600` `376:-13440` `384:-40320` `395:-134400`|
|1|14|`396:-26880` `400:-201600` `402:-134400` `403:-67200` `404:-134400`|
|1|15|`11:4480` `15:13440` `113:-13440` `118:-40320` `222:-40320` `225:-67200` `229:-120960` `235:-201600` `376:-13440` `377:-67200` `379:-44800` `384:-40320`|
|1|15|`391:-201600` `392:-134400`|
|1|16|`11:2240` `41:-2240` `111:-13440` `113:-26880` `220:-6720` `222:-40320` `225:-33600`|
|2|2|`5:960` `21:-20160` `28:-20160` `29:40320` `32:60480` `33:6720` `34:24960` `35:-2880` `36:241920` `37:250560` `38:10560` `39:95040`|
|2|2|`40:960` `621:-100800` `622:-403200` `623:-100800` `652:-907200` `653:-3628800` `654:-907200` `658:-907200` `659:-3628800` `662:-100800` `735:-20160` `736:-80640`|
|2|2|`737:-20160` `738:-181440` `739:-725760` `743:-181440` `763:-80640` `764:-322560` `765:-80640` `774:-967680` `775:-2903040` `776:-725760` `777:-725760` `778:-80640`|
|2|2|`779:-2903040` `788:-362880` `804:-1451520` `805:-362880` `813:-967680` `815:-241920` `828:-725760` `829:-2903040` `831:-1451520` `833:-5806080` `834:-1451520` `835:-362880`|
|2|2|`836:-120960` `837:-725760` `838:-2177280` `839:-483840` `840:-2903040` `841:-2903040` `863:-6720` `864:-134400` `865:-94080` `866:-13440` `867:-362880` `868:-604800`|
|2|2|`869:-60480` `870:-241920` `873:-161280` `875:-161280` `876:-40320` `877:-120960` `879:-483840` `880:-967680` `881:-120960` `882:-120960` `883:-362880` `889:-483840`|
|2|2|`890:-120960` `893:-60480` `899:-241920` `900:-403200` `901:-40320` `902:-161280` `905:-967680` `906:-1451520` `907:-725760` `908:-241920` `910:-967680` `913:-241920`|
|2|2|`914:-241920` `915:-483840` `918:-20160` `919:-120960` `920:-120960` `921:-967680` `922:-967680` `923:-967680` `924:-483840` `927:-17280` `928:-51840` `929:-8640`|
|2|2|`930:-34560` `931:-138240` `932:-8640` `933:-34560` `934:-69120` `935:-69120` `936:-17280` `937:-138240`|
|2|3|`28:40320` `31:120960` `32:40320` `35:46080` `36:138240` `38:5760` `656:-604800` `657:-604800` `660:-1814400` `661:-2419200` `741:-120960` `742:-120960`|
|2|3|`744:-362880` `774:-241920` `775:-725760` `776:-241920` `777:-241920` `779:-725760` `789:-241920` `814:-241920` `830:-1451520` `831:-725760` `832:-2177280` `833:-1451520`|
|2|3|`834:-241920` `837:-725760` `840:-1451520` `872:-40320` `873:-40320` `874:-40320` `878:-241920` `879:-241920` `880:-241920` `882:-120960` `889:-241920` `890:-120960`|
|2|3|`894:-120960` `899:-40320` `905:-241920` `910:-483840` `913:-362880` `914:-120960` `921:-483840` `930:-34560` `931:-69120` `933:-17280`|
|2|4|`27:40320` `29:362880` `32:362880` `33:40320` `34:53760` `36:483840` `37:483840` `38:13440` `39:120960` `402:-403200` `403:-403200` `456:-3628800`|
|2|4|`457:-3628800` `480:-3628800` `481:-3628800` `489:-403200` `571:-201600` `572:-201600` `575:-201600` `576:-201600` `580:-1814400` `581:-1814400` `587:-1814400` `605:-2419200`|
|2|4|`606:-1814400` `610:-1814400` `615:-1814400` `620:-604800` `621:-403200` `622:-403200` `624:-604800` `630:-604800` `632:-1814400` `633:-1814400` `636:-201600` `652:-3628800`|
|2|4|`653:-3628800` `658:-3628800` `659:-3628800` `662:-403200` `677:-3628800` `678:-3628800` `679:-3628800` `681:-3628800` `690:-3628800` `692:-3628800` `695:-1814400` `697:-1814400`|
|2|4|`705:-604800` `706:-1814400` `707:-1814400` `709:-1814400` `710:-1814400` `713:-3628800` `715:-3628800` `717:-604800` `735:-80640` `736:-80640` `738:-725760` `739:-725760`|
|2|4|`743:-725760` `749:-725760` `750:-725760` `752:-604800` `754:-241920` `760:-362880` `763:-120960` `764:-322560` `765:-80640` `768:-241920` `770:-241920` `772:-120960`|
|2|4|`774:-1572480` `775:-2177280` `776:-725760` `777:-1209600` `778:-120960` `779:-1814400` `780:-1451520` `782:-2177280` `783:-1451520` `786:-1451520` `787:-1451520` `788:-725760`|
|2|4|`790:-362880` `791:-725760` `793:-1088640` `795:-725760` `796:-725760` `799:-1451520` `800:-1088640` `802:-725760` `803:-1451520` `804:-725760` `813:-1209600` `817:-1088640`|
|2|4|`818:-725760` `821:-362880` `822:-120960` `823:-483840` `825:-362880` `828:-1451520` `829:-1451520` `831:-2903040` `832:-725760` `833:-4354560` `835:-725760` `836:-362880`|
|2|4|`837:-1814400` `838:-2903040` `839:-241920` `840:-2177280` `841:-1451520` `847:-725760` `849:-1451520` `851:-725760` `853:-1451520` `854:-1451520` `856:-725760` `857:-362880`|
|2|4|`858:-1451520` `859:-362880` `860:-725760` `864:-13440` `865:-13440` `867:-120960` `868:-241920` `869:-120960` `870:-120960` `873:-80640` `875:-80640` `877:-362880`|
|2|4|`878:-241920` `879:-483840` `880:-725760` `881:-120960` `882:-241920` `883:-362880` `889:-241920` `890:-120960` `891:-241920` `892:-241920` `895:-120960` `897:-120960`|
|2|4|`898:-120960` `899:-40320` `900:-40320` `902:-40320` `904:-40320` `905:-483840` `906:-241920` `907:-362880` `908:-241920` `910:-483840` `911:-483840` `912:-483840`|
|2|4|`913:-120960` `914:-120960` `915:-241920` `916:-120960` `917:-483840`|
|2|5|`21:60480` `24:181440` `26:60480` `28:262080` `31:786240` `32:262080` `35:141120` `36:423360` `38:20160` `472:-604800` `473:-2419200` `474:-604800`|
|2|5|`484:-1814400` `485:-7257600` `486:-1814400` `487:-604800` `491:-2419200` `584:-604800` `585:-1209600` `586:-302400` `590:-907200` `591:-3628800` `592:-302400` `621:-302400`|
|2|5|`624:-2116800` `625:-604800` `627:-1209600` `631:-302400` `634:-3628800` `635:-2116800` `652:-907200` `655:-604800` `656:-4838400` `657:-1209600` `660:-7257600` `661:-4233600`|
|2|5|`672:-1814400` `689:-7257600` `690:-5443200` `691:-1814400` `693:-907200` `694:-907200` `695:-907200` `696:-1209600` `697:-604800` `703:-2419200` `704:-604800` `708:-3628800`|
|2|5|`711:-7257600` `712:-3628800` `713:-3628800` `714:-3024000` `715:-2419200` `716:-907200` `717:-907200` `735:-60480` `738:-181440` `740:-241920` `741:-967680` `742:-241920`|
|2|5|`744:-1451520` `757:-181440` `763:-120960` `764:-483840` `765:-120960` `766:-120960` `767:-483840` `768:-604800` `769:-120960` `771:-241920` `772:-181440` `774:-967680`|
|2|5|`775:-2177280` `776:-604800` `777:-604800` `778:-120960` `779:-1814400` `784:-725760` `785:-2903040` `786:-1451520` `788:-362880` `789:-1451520` `790:-362880` `791:-2177280`|
|2|5|`792:-725760` `793:-725760` `794:-483840` `795:-241920` `797:-1451520` `798:-241920` `801:-2903040` `802:-1451520` `803:-1451520` `807:-1451520` `808:-725760` `810:-181440`|
|2|5|`813:-1209600` `814:-483840` `815:-120960` `816:-1088640` `817:-1088640` `819:-544320` `821:-362880` `823:-120960` `824:-120960` `825:-181440` `830:-2903040` `831:-1451520`|
|2|5|`832:-2177280` `833:-4354560` `834:-483840` `837:-362880` `840:-725760` `843:-2903040` `847:-725760` `848:-1451520` `849:-1451520` `852:-1451520` `853:-1451520` `854:-1451520`|
|2|5|`856:-725760` `857:-362880` `859:-362880` `869:-60480` `871:-20160` `872:-241920` `873:-80640` `874:-80640` `875:-120960` `876:-20160` `878:-483840` `879:-241920`|
|2|5|`880:-483840` `882:-60480` `885:-60480` `889:-241920` `890:-60480` `891:-725760` `893:-60480` `894:-120960` `895:-120960` `896:-120960` `897:-120960` `898:-120960`|
|2|5|`899:-80640` `901:-20160` `903:-161280` `904:-40320` `905:-483840` `909:-967680` `910:-967680` `911:-483840` `912:-483840` `916:-120960`|
|2|6|`22:120960` `24:120960` `28:161280` `29:161280` `34:40320` `482:-1209600` `483:-1209600` `488:-1209600` `490:-1209600` `588:-604800` `589:-604800` `593:-604800`|
|2|6|`605:-604800` `606:-604800` `632:-604800` `633:-604800` `668:-1209600` `679:-1814400` `680:-1814400` `681:-1209600` `702:-1209600` `706:-3628800` `707:-1209600` `752:-241920`|
|2|6|`753:-120960` `760:-362880` `781:-483840` `782:-2177280` `783:-483840` `792:-362880` `796:-725760` `799:-1088640` `802:-725760` `814:-120960` `832:-725760` `837:-362880`|
|2|6|`840:-725760` `872:-40320` `878:-241920` `879:-120960` `880:-241920` `899:-40320` `905:-241920`|
|2|7|`27:33600` `29:302400` `32:302400` `33:33600` `34:20160` `36:181440` `37:181440` `225:-201600` `256:-1814400` `276:-1814400` `308:-201600` `377:-201600`|
|2|7|`379:-134400` `385:-1814400` `392:-1612800` `394:-1209600` `402:-403200` `404:-403200` `407:-1814400` `415:-1209600` `434:-134400` `450:-2419200` `451:-2419200` `456:-3628800`|
|2|7|`469:-2419200` `470:-1209600` `480:-3628800` `489:-806400` `490:-1209600` `514:-1209600` `516:-2419200` `529:-403200` `566:-907200` `568:-604800` `571:-201600` `575:-100800`|
|2|7|`576:-201600` `579:-302400` `580:-1814400` `587:-3628800` `589:-2419200` `592:-302400` `593:-604800` `605:-1209600` `606:-1209600` `609:-1814400` `610:-907200` `613:-604800`|
|2|7|`614:-1209600` `620:-604800` `621:-100800` `624:-302400` `625:-100800` `632:-907200` `633:-1814400` `636:-201600` `638:-907200` `642:-1814400` `648:-302400` `652:-907200`|
|2|7|`658:-907200` `661:-604800` `662:-201600` `677:-1814400` `679:-1814400` `681:-1209600` `687:-604800` `688:-1209600` `690:-1814400` `695:-907200` `696:-302400` `697:-604800`|
|2|7|`705:-302400` `706:-907200` `709:-907200` `713:-1814400` `716:-302400` `717:-604800` `726:-604800` `735:-20160` `738:-181440` `743:-362880` `744:-362880` `749:-362880`|
|2|7|`752:-120960` `760:-362880` `761:-362880` `762:-181440` `764:-40320` `768:-120960` `771:-60480` `772:-60480` `775:-362880` `779:-725760` `780:-725760` `782:-725760`|
|2|7|`786:-725760` `793:-362880` `799:-725760` `800:-725760` `802:-725760` `803:-1451520` `817:-362880` `819:-181440` `821:-362880` `825:-362880` `827:-362880`|
|2|8|`21:201600` `24:604800` `26:201600` `28:322560` `31:967680` `32:322560` `35:120960` `36:362880` `263:-1209600` `264:-1209600` `292:-3628800` `293:-3628800`|
|2|8|`310:-1209600` `311:-1209600` `397:-2016000` `398:-1209600` `402:-806400` `404:-806400` `406:-806400` `409:-3628800` `410:-3628800` `416:-1209600` `425:-2419200` `428:-2419200`|
|2|8|`432:-806400` `435:-806400` `456:-2419200` `463:-4838400` `468:-4838400` `469:-2419200` `472:-2419200` `473:-4032000` `479:-1612800` `484:-7257600` `485:-7257600` `487:-2419200`|
|2|8|`491:-2419200` `503:-4838400` `507:-2419200` `513:-2419200` `515:-4838400` `519:-2419200` `520:-2419200` `523:-806400` `524:-2419200` `525:-4032000` `528:-806400` `530:-1612800`|
|2|8|`571:-201600` `577:-604800` `578:-604800` `579:-201600` `580:-604800` `584:-1814400` `585:-1612800` `590:-3628800` `591:-3628800` `592:-1209600` `597:-604800` `602:-1209600`|
|2|8|`607:-3628800` `608:-3628800` `609:-604800` `610:-1814400` `613:-1814400` `617:-1209600` `618:-604800` `621:-604800` `622:-403200` `624:-3024000` `625:-1209600` `627:-1814400`|
|2|8|`628:-604800` `629:-1612800` `630:-1008000` `634:-3628800` `635:-3024000` `637:-1814400` `639:-3628800` `641:-3628800` `643:-1814400` `644:-1814400` `645:-1814400` `646:-1209600`|
|2|8|`647:-3024000` `650:-604800` `652:-1814400` `653:-1209600` `655:-604800` `656:-3628800` `657:-806400` `660:-1814400` `661:-604800` `665:-2419200` `672:-3628800` `673:-604800`|
|2|8|`674:-1209600` `675:-1209600` `684:-2419200` `687:-604800` `689:-3628800` `690:-5443200` `691:-604800` `692:-2419200` `693:-1814400` `694:-1814400` `695:-1814400` `696:-604800`|
|2|8|`697:-1814400` `700:-1209600` `701:-2419200` `703:-1209600` `704:-201600` `708:-1814400` `711:-3628800` `712:-1814400` `713:-1814400` `714:-1209600` `715:-1209600` `721:-1209600`|
|2|8|`724:-1209600` `726:-604800` `740:-120960` `741:-483840` `742:-120960` `744:-362880` `749:-362880` `756:-362880` `761:-725760` `766:-241920` `767:-483840` `768:-120960`|
|2|8|`771:-120960` `773:-241920` `774:-241920` `775:-725760` `776:-241920` `777:-241920` `779:-725760` `784:-1451520` `785:-2903040` `788:-725760` `789:-483840` `790:-1088640`|
|2|8|`791:-2177280` `792:-725760` `793:-725760` `794:-241920` `795:-483840` `797:-483840` `801:-2903040` `802:-725760` `803:-725760` `804:-725760` `807:-725760` `808:-725760`|
|2|8|`809:-725760` `813:-241920` `814:-120960` `816:-1451520` `817:-362880` `818:-2177280` `820:-1451520` `821:-362880` `823:-241920` `824:-120960` `825:-362880` `826:-725760`|
|2|8|`827:-362880`|
|2|9|`15:100800` `18:100800` `22:463680` `24:463680` `28:342720` `29:342720` `34:60480` `290:-604800` `294:-604800` `304:-2419200` `305:-604800` `312:-2419200`|
|2|9|`313:-604800` `402:-403200` `404:-403200` `411:-604800` `412:-604800` `413:-2419200` `414:-604800` `417:-2419200` `425:-1612800` `428:-1612800` `432:-403200` `435:-403200`|
|2|9|`440:-806400` `456:-1209600` `458:-1209600` `459:-1612800` `465:-1209600` `469:-806400` `470:-403200` `476:-3225600` `478:-806400` `482:-4838400` `483:-4838400` `488:-1209600`|
|2|9|`490:-1209600` `503:-4838400` `508:-4838400` `512:-1612800` `515:-3225600` `516:-1612800` `523:-2419200` `527:-403200` `528:-806400` `529:-403200` `564:-100800` `571:-201600`|
|2|9|`575:-302400` `576:-403200` `579:-302400` `580:-604800` `582:-1209600` `583:-806400` `588:-2419200` `589:-2419200` `593:-604800` `598:-302400` `601:-604800` `602:-806400`|
|2|9|`603:-403200` `605:-2419200` `606:-1814400` `609:-604800` `610:-907200` `611:-907200` `612:-1209600` `613:-604800` `614:-604800` `617:-1209600` `620:-1008000` `626:-2419200`|
|2|9|`628:-806400` `631:-100800` `632:-1512000` `633:-1512000` `636:-302400` `639:-7257600` `640:-1209600` `641:-2419200` `642:-1209600` `646:-1209600` `648:-403200` `665:-2419200`|
|2|9|`668:-2419200` `669:-806400` `674:-1209600` `679:-3628800` `680:-3628800` `681:-3628800` `684:-2419200` `688:-1209600` `691:-604800` `696:-302400` `702:-604800` `704:-201600`|
|2|9|`706:-1814400` `707:-604800` `714:-604800` `716:-302400` `717:-302400` `746:-60480` `749:-362880` `751:-241920` `752:-483840` `753:-483840` `754:-241920` `757:-181440`|
|2|9|`760:-362880` `761:-362880` `762:-181440` `767:-241920` `771:-60480` `781:-967680` `782:-2903040` `783:-967680` `785:-1451520` `791:-725760` `796:-725760` `797:-483840`|
|2|9|`799:-1088640` `801:-1451520` `802:-725760` `803:-725760` `816:-362880` `822:-120960` `823:-120960` `825:-181440`|
|2|10|`15:67200` `21:107520` `27:40320` `277:-403200` `309:-403200` `392:-268800` `408:-403200` `415:-672000` `446:-806400` `451:-806400` `488:-806400` `490:-806400`|
|2|10|`568:-201600` `583:-201600` `588:-604800` `589:-1209600` `593:-201600` `632:-604800` `633:-604800` `657:-201600` `661:-604800` `741:-120960` `744:-362880` `777:-120960`|
|2|10|`779:-362880`|
|2|11|`21:100800` `24:302400` `26:100800` `28:100800` `31:302400` `32:100800` `123:-302400` `134:-907200` `150:-302400` `223:-907200` `225:-302400` `230:-1814400`|
|2|11|`235:-907200` `242:-604800` `246:-302400` `256:-907200` `260:-1814400` `263:-1209600` `264:-604800` `282:-1814400` `283:-907200` `292:-3628800` `302:-907200` `303:-907200`|
|2|11|`310:-1209600` `313:-302400` `319:-604800` `377:-100800` `380:-403200` `385:-302400` `389:-302400` `391:-604800` `393:-2419200` `394:-1209600` `397:-1612800` `398:-201600`|
|2|11|`402:-403200` `403:-201600` `404:-403200` `406:-806400` `409:-3628800` `414:-302400` `416:-2419200` `417:-2419200` `425:-1209600` `427:-302400` `428:-1209600` `430:-1209600`|
|2|11|`431:-1209600` `432:-403200` `435:-806400` `437:-806400` `444:-1209600` `449:-302400` `456:-1209600` `457:-604800` `463:-2419200` `466:-604800` `467:-604800` `472:-604800`|
|2|11|`473:-806400` `474:-100800` `484:-1814400` `486:-302400` `487:-604800` `503:-2419200` `505:-604800` `507:-1209600` `511:-604800` `519:-1209600` `520:-1209600` `522:-302400`|
|2|11|`523:-403200` `525:-806400` `566:-302400` `570:-302400` `577:-302400` `578:-201600` `579:-100800` `584:-302400` `590:-907200` `591:-1209600` `592:-604800` `607:-1814400`|
|2|11|`610:-907200` `615:-604800` `617:-604800` `618:-604800` `619:-604800` `624:-302400` `627:-604800` `630:-201600` `631:-100800` `634:-907200` `635:-604800` `637:-907200`|
|2|11|`639:-1814400` `643:-907200` `644:-604800` `645:-907200` `646:-604800` `647:-1209600` `650:-302400` `651:-302400`|
|2|12|`15:201600` `18:201600` `22:403200` `24:403200` `28:201600` `29:201600` `138:-604800` `140:-604800` `152:-604800` `155:-604800` `236:-604800` `238:-1209600`|
|2|12|`239:-1209600` `240:-1209600` `244:-1209600` `245:-1209600` `248:-604800` `261:-1209600` `266:-1209600` `284:-1814400` `287:-1814400` `290:-2419200` `291:-604800` `294:-2419200`|
|2|12|`296:-3024000` `300:-1814400` `303:-604800` `304:-2419200` `307:-604800` `312:-2419200` `316:-1209600` `319:-604800` `383:-201600` `395:-403200` `399:-1612800` `400:-1008000`|
|2|12|`401:-201600` `402:-806400` `403:-201600` `404:-806400` `405:-1612800` `411:-2419200` `412:-2419200` `413:-2419200` `417:-2419200` `418:-2419200` `421:-2419200` `423:-604800`|
|2|12|`425:-1612800` `426:-806400` `427:-201600` `428:-1612800` `429:-6451200` `431:-806400` `432:-806400` `433:-806400` `435:-806400` `436:-1612800` `437:-806400` `438:-1209600`|
|2|12|`440:-1612800` `441:-403200` `454:-604800` `456:-2419200` `457:-604800` `458:-1209600` `459:-2016000` `461:-2419200` `465:-2419200` `466:-604800` `467:-604800` `469:-1612800`|
|2|12|`470:-806400` `476:-1612800` `477:-403200` `482:-1209600` `483:-1209600` `497:-1209600` `503:-2419200` `505:-604800` `508:-2419200` `511:-604800` `512:-806400` `515:-1612800`|
|2|12|`516:-806400` `573:-201600` `578:-201600` `582:-604800` `583:-604800` `588:-604800` `589:-604800` `594:-1209600` `597:-604800` `601:-1209600` `602:-403200` `603:-403200`|
|2|12|`605:-604800` `606:-604800` `608:-2419200` `610:-1814400` `611:-1814400` `612:-604800` `613:-1814400` `614:-604800` `615:-604800` `617:-604800` `619:-604800` `620:-201600`|
|2|12|`626:-1209600` `628:-403200` `629:-403200` `632:-604800` `633:-604800` `639:-3628800` `640:-604800` `641:-2419200` `642:-604800` `644:-604800` `646:-604800` `647:-604800`|
|2|12|`648:-201600`|
|2|13|`11:33600` `15:168000` `21:168000` `27:33600` `133:-100800` `149:-403200` `157:-100800` `225:-100800` `231:-201600` `235:-403200` `243:-806400` `246:-302400`|
|2|13|`253:-302400` `256:-302400` `257:-403200` `275:-1209600` `277:-1612800` `282:-1209600` `305:-302400` `309:-403200` `313:-302400` `377:-100800` `379:-134400` `382:-100800`|
|2|13|`385:-302400` `386:-806400` `390:-403200` `391:-403200` `392:-1075200` `394:-403200` `408:-1612800` `413:-1612800` `414:-100800` `415:-1075200` `417:-1612800` `432:-403200`|
|2|13|`434:-134400` `435:-403200` `441:-403200` `444:-1209600` `446:-1612800` `451:-1612800` `474:-100800` `486:-302400` `488:-403200` `490:-403200` `564:-100800` `566:-302400`|
|2|13|`567:-302400` `568:-403200` `585:-403200` `588:-1209600` `589:-1209600` `591:-1209600` `593:-403200` `625:-100800` `632:-302400` `633:-302400` `635:-302400` `636:-100800`|
|2|14|`15:60480` `18:60480` `22:100800` `24:100800` `46:-60480` `52:-60480` `113:-120960` `114:-302400` `118:-120960` `119:-302400` `122:-241920` `130:-362880`|
|2|14|`132:-362880` `138:-604800` `139:-120960` `146:-241920` `152:-604800` `154:-120960` `221:-60480` `222:-120960` `224:-604800` `229:-181440` `232:-907200` `234:-907200`|
|2|14|`236:-302400` `237:-60480` `238:-1209600` `240:-604800` `241:-302400` `244:-1209600` `245:-604800` `248:-302400` `249:-302400` `252:-362880` `255:-181440` `261:-604800`|
|2|14|`262:-120960` `271:-362880` `284:-907200` `285:-181440` `287:-907200` `289:-181440` `290:-604800` `291:-302400` `294:-604800` `296:-604800` `303:-302400` `378:-100800`|
|2|14|`387:-604800` `389:-302400` `399:-806400` `400:-201600` `401:-201600` `410:-604800` `411:-604800` `412:-604800` `418:-1209600` `419:-302400` `421:-1209600` `423:-302400`|
|2|14|`424:-302400` `425:-403200` `426:-403200` `427:-100800` `428:-403200` `429:-806400` `431:-403200`|
|2|15|`11:40320` `15:107520` `21:67200` `49:-40320` `54:-40320` `113:-80640` `117:-201600` `118:-80640` `120:-201600` `127:-241920` `130:-241920` `133:-403200`|
|2|15|`137:-241920` `146:-241920` `149:-403200` `222:-120960` `225:-201600` `228:-604800` `229:-120960` `231:-806400` `232:-604800` `235:-403200` `239:-604800` `243:-806400`|
|2|15|`245:-604800` `246:-201600` `251:-120960` `252:-362880` `253:-604800` `256:-604800` `257:-403200` `262:-120960` `271:-362880` `275:-604800` `277:-403200` `282:-604800`|
|2|15|`383:-201600` `386:-403200` `387:-604800` `390:-806400` `391:-201600` `392:-268800` `394:-806400` `398:-201600` `408:-403200` `410:-604800` `413:-806400` `414:-201600`|
|2|15|`415:-268800` `417:-806400`|
|2|16|`11:6720` `15:20160` `41:-6720` `42:-20160` `44:-20160` `47:-20160` `49:-40320` `111:-20160` `112:-120960` `113:-40320` `115:-120960` `117:-201600`|
|2|16|`118:-40320` `121:-20160` `125:-60480` `127:-120960` `130:-120960` `133:-100800` `221:-60480` `226:-181440` `228:-302400` `229:-60480` `231:-201600` `232:-302400`|
|2|16|`235:-100800`|
|3|3|`4:480` `21:-10080` `28:-13440` `30:60480` `31:80640` `35:25920` `658:-907200` `659:-1209600` `662:-302400` `743:-181440` `828:-725760` `829:-1451520`|
|3|3|`836:-241920` `877:-120960` `893:-60480` `908:-241920` `918:-20160` `929:-8640`|
|3|4|`21:40320` `24:120960` `26:40320` `28:134400` `31:403200` `32:134400` `35:40320` `36:120960` `458:-403200` `459:-1209600` `482:-1209600` `483:-3628800`|
|3|4|`488:-403200` `490:-1209600` `582:-403200` `583:-604800` `588:-604800` `589:-1814400` `593:-201600` `612:-604800` `655:-403200` `656:-3024000` `657:-403200` `660:-3628800`|
|3|4|`661:-1209600` `668:-1209600` `680:-3628800` `694:-1814400` `702:-604800` `703:-1209600` `704:-201600` `706:-1814400` `712:-3628800` `716:-604800` `740:-161280` `741:-604800`|
|3|4|`742:-80640` `744:-725760` `751:-80640` `753:-362880` `763:-40320` `764:-120960` `765:-40320` `774:-362880` `775:-362880` `776:-120960` `777:-120960` `781:-483840`|
|3|4|`782:-1451520` `789:-967680` `792:-1451520` `796:-241920` `797:-967680` `798:-80640` `802:-1451520` `813:-483840` `814:-362880` `815:-80640` `823:-120960` `824:-120960`|
|3|4|`830:-1451520` `831:-725760` `832:-725760` `833:-1451520` `834:-241920` `837:-725760` `848:-1451520` `856:-725760` `871:-13440` `872:-201600` `873:-40320` `874:-40320`|
|3|4|`875:-80640` `876:-13440` `878:-241920` `879:-241920` `880:-241920` `894:-120960` `899:-80640` `903:-161280` `904:-80640` `910:-483840` `913:-120960` `914:-120960`|
|3|5|`23:181440` `25:241920` `26:30240` `30:241920` `31:322560` `35:60480` `484:-1814400` `485:-2419200` `486:-302400` `487:-1814400` `491:-2419200` `590:-907200`|
|3|5|`591:-1209600` `592:-907200` `652:-1814400` `653:-2419200` `654:-302400` `658:-1814400` `659:-2419200` `672:-1814400` `689:-3628800` `690:-1814400` `691:-604800` `692:-1209600`|
|3|5|`693:-1814400` `708:-3628800` `709:-1814400` `710:-1209600` `738:-362880` `739:-483840` `743:-362880` `757:-181440` `784:-725760` `785:-1451520` `786:-725760` `787:-483840`|
|3|5|`788:-725760` `790:-725760` `791:-725760` `800:-725760` `804:-1451520` `805:-241920` `808:-725760` `810:-181440` `818:-725760` `819:-181440` `828:-725760` `829:-1451520`|
|3|5|`847:-1451520` `851:-725760` `869:-120960` `870:-241920` `877:-120960` `885:-60480` `891:-483840` `892:-241920` `893:-60480` `908:-241920`|
|3|6|`23:120960` `24:60480` `28:80640` `480:-1209600` `481:-604800` `489:-403200` `587:-604800` `677:-1814400` `678:-1209600` `705:-604800` `780:-725760` `790:-362880`|
|3|6|`800:-362880` `836:-120960` `877:-120960`|
|3|7|`21:67200` `24:201600` `26:67200` `28:40320` `31:120960` `32:40320` `257:-403200` `277:-1209600` `309:-403200` `386:-672000` `408:-1209600` `415:-403200`|
|3|7|`446:-1612800` `458:-806400` `459:-806400` `482:-2419200` `488:-1612800` `512:-1612800` `527:-268800` `567:-201600` `582:-604800` `583:-403200` `588:-2419200` `593:-403200`|
|3|7|`612:-604800` `632:-604800` `640:-1209600` `648:-201600` `655:-201600` `656:-1209600` `657:-604800` `668:-1209600` `669:-806400` `702:-604800` `703:-1209600` `704:-604800`|
|3|7|`740:-40320` `741:-241920` `742:-40320` `751:-80640` `753:-120960` `777:-120960` `781:-483840` `789:-241920` `796:-241920` `797:-967680` `798:-80640` `814:-120960`|
|3|7|`823:-120960` `824:-120960`|
|3|8|`16:100800` `19:134400` `20:16800` `23:362880` `25:483840` `26:60480` `30:181440` `31:241920` `290:-604800` `294:-1814400` `304:-806400` `305:-100800`|
|3|8|`312:-2419200` `313:-302400` `411:-604800` `412:-1814400` `413:-806400` `414:-100800` `417:-2419200` `440:-806400` `465:-2419200` `472:-1209600` `473:-2419200` `474:-201600`|
|3|8|`476:-1612800` `478:-268800` `479:-537600` `484:-3628800` `485:-4838400` `486:-604800` `508:-4838400` `519:-2419200` `523:-806400` `524:-1612800` `564:-100800` `584:-1209600`|
|3|8|`585:-1209600` `586:-100800` `590:-1814400` `591:-2419200` `598:-302400` `601:-604800` `603:-403200` `611:-1814400` `621:-302400` `622:-403200` `623:-50400` `624:-604800`|
|3|8|`625:-100800` `626:-1209600` `627:-1814400` `628:-201600` `629:-806400` `630:-403200` `631:-302400` `639:-3628800` `643:-1814400` `644:-1209600` `652:-907200` `653:-1209600`|
|3|8|`654:-151200` `672:-1814400` `674:-2419200` `689:-3628800` `690:-1814400` `691:-604800` `692:-1209600` `696:-604800` `700:-1209600` `735:-60480` `736:-80640` `737:-10080`|
|3|8|`738:-181440` `739:-241920` `746:-60480` `757:-181440` `766:-241920` `767:-967680` `768:-241920` `769:-80640` `770:-161280` `771:-120960` `773:-241920` `784:-725760`|
|3|8|`785:-1451520` `786:-725760` `787:-483840` `791:-725760` `808:-725760` `818:-725760` `819:-181440`|
|3|9|`17:201600` `19:100800` `23:322560` `24:161280` `28:120960` `292:-1209600` `293:-604800` `310:-1209600` `311:-604800` `409:-1209600` `410:-604800` `416:-1209600`|
|3|9|`456:-2419200` `457:-1209600` `463:-2419200` `468:-1612800` `469:-806400` `480:-2419200` `481:-1209600` `507:-2419200` `513:-1612800` `514:-806400` `580:-1209600` `581:-604800`|
|3|9|`587:-1209600` `597:-604800` `607:-1814400` `608:-1209600` `609:-604800` `610:-1814400` `613:-604800` `615:-1209600` `618:-604800` `637:-1814400` `638:-604800` `673:-604800`|
|3|9|`677:-1814400` `678:-1209600` `687:-1209600` `749:-725760` `750:-483840` `756:-362880` `761:-725760` `780:-725760` `790:-362880` `800:-362880`|
|3|10|`16:50400` `21:47040` `276:-302400` `308:-100800` `407:-302400` `450:-806400` `489:-268800` `587:-604800` `592:-302400` `662:-100800` `743:-181440`|
|3|11|`16:100800` `19:134400` `20:16800` `23:100800` `25:134400` `26:16800` `133:-302400` `149:-403200` `157:-50400` `231:-604800` `243:-806400` `246:-100800`|
|3|11|`253:-604800` `275:-1209600` `290:-1209600` `291:-604800` `304:-1612800` `305:-403200` `307:-403200` `382:-100800` `390:-806400` `411:-1209600` `413:-3225600` `414:-201600`|
|3|11|`426:-806400` `432:-268800` `433:-537600` `440:-806400` `441:-806400` `472:-604800` `473:-1612800` `474:-302400` `476:-1612800` `477:-403200` `478:-268800` `479:-537600`|
|3|11|`564:-100800` `584:-302400` `585:-806400` `586:-50400` `601:-604800` `603:-403200` `625:-100800` `626:-1209600` `627:-604800` `628:-201600` `629:-806400` `631:-100800`|
|3|12|`12:67200` `14:33600` `17:268800` `19:134400` `23:201600` `24:100800` `138:-201600` `140:-100800` `152:-604800` `155:-302400` `238:-403200` `239:-201600`|
|3|12|`244:-1209600` `245:-604800` `261:-604800` `263:-806400` `264:-604800` `266:-403200` `287:-1814400` `292:-2419200` `293:-1209600` `300:-1209600` `302:-604800` `383:-201600`|
|3|12|`397:-1612800` `398:-604800` `399:-806400` `401:-201600` `402:-403200` `403:-201600` `404:-268800` `405:-537600` `406:-806400` `409:-2419200` `410:-1209600` `421:-2419200`|
|3|12|`429:-1612800` `430:-806400` `454:-604800` `456:-1209600` `457:-604800` `463:-2419200` `466:-1209600` `468:-1612800` `469:-806400` `571:-201600` `572:-100800` `573:-201600`|
|3|12|`577:-604800` `578:-806400` `579:-201600` `580:-604800` `581:-302400` `597:-604800` `607:-1814400` `608:-1209600` `609:-604800` `613:-604800` `618:-604800`|
|3|13|`12:50400` `16:100800` `21:50400` `134:-151200` `150:-151200` `230:-302400` `242:-302400` `256:-604800` `260:-604800` `276:-604800` `283:-604800` `385:-604800`|
|3|13|`389:-302400` `393:-806400` `394:-806400` `407:-604800` `416:-806400` `449:-302400` `450:-806400` `566:-604800` `570:-302400` `587:-604800` `592:-302400`|
|3|14|`12:40320` `14:20160` `17:67200` `19:33600` `49:-40320` `54:-20160` `117:-201600` `120:-100800` `127:-241920` `137:-161280` `138:-403200` `139:-80640`|
|3|14|`140:-201600` `228:-604800` `238:-806400` `239:-806400` `241:-201600` `251:-120960` `261:-604800` `262:-241920` `263:-403200` `264:-403200` `266:-403200` `383:-201600`|
|3|14|`397:-403200` `398:-403200` `399:-806400` `401:-201600` `405:-537600` `406:-268800`|
|3|15|`8:10080` `12:47040` `16:50400` `46:-10080` `52:-30240` `114:-50400` `119:-151200` `122:-80640` `123:-100800` `132:-241920` `134:-302400` `221:-60480`|
|3|15|`223:-403200` `224:-201600` `225:-100800` `230:-604800` `234:-604800` `255:-181440` `256:-302400` `260:-604800` `377:-100800` `378:-100800` `380:-268800` `385:-302400`|
|3|15|`389:-302400` `393:-806400`|
|3|16|`8:3360` `12:10080` `42:-10080` `44:-13440` `46:-20160` `112:-80640` `114:-100800` `121:-20160` `122:-80640` `123:-50400` `221:-60480` `223:-100800`|
|3|16|`224:-201600`|
|4|4|`5:2240` `15:-33600` `21:-26880` `22:-40320` `24:-40320` `27:98560` `28:-20160` `29:893760` `31:-20160` `32:900480` `33:100800` `34:114240`|
|4|4|`36:1028160` `37:1028160` `38:33600` `39:302400` `40:2240` `236:-604800` `237:-201600` `284:-5443200` `285:-1814400` `320:-5443200` `321:-1814400` `329:-604800`|
|4|4|`395:-806400` `396:-268800` `402:-1612800` `403:-537600` `418:-3628800` `419:-1209600` `425:-2419200` `427:-806400` `438:-7257600` `439:-2419200` `442:-7257600` `456:-14515200`|
|4|4|`457:-4838400` `480:-14515200` `481:-4838400` `489:-1612800` `492:-7257600` `493:-2419200` `503:-14515200` `505:-4838400` `520:-3628800` `522:-1209600` `532:-1209600` `534:-7257600`|
|4|4|`535:-2419200` `542:-7257600` `543:-2419200` `548:-1209600` `561:-100800` `562:-33600` `563:-907200` `571:-806400` `572:-268800` `575:-201600` `576:-67200` `580:-7257600`|
|4|4|`581:-2419200` `587:-7257600` `594:-3628800` `595:-1209600` `602:-2419200` `604:-806400` `605:-3024000` `606:-604800` `610:-3628800` `615:-1209600` `620:-403200` `621:-403200`|
|4|4|`622:-537600` `623:-134400` `624:-1209600` `630:-403200` `632:-2419200` `633:-604800` `634:-907200` `636:-201600` `637:-907200` `638:-302400` `639:-3628800` `644:-1209600`|
|4|4|`645:-907200` `648:-201600` `650:-302400` `651:-302400` `652:-3628800` `653:-4838400` `654:-1209600` `658:-3628800` `659:-4838400` `662:-403200` `663:-7257600` `665:-14515200`|
|4|4|`666:-4838400` `675:-3628800` `676:-1209600` `677:-7257600` `678:-9676800` `679:-7257600` `681:-2419200` `685:-7257600` `690:-7257600` `692:-2419200` `695:-3628800` `697:-1209600`|
|4|4|`705:-1512000` `706:-5443200` `707:-1209600` `708:-3628800` `709:-3628800` `710:-2419200` `711:-3628800` `713:-10886400` `715:-2419200` `716:-604800` `717:-1209600` `719:-3628800`|
|4|4|`721:-1209600` `722:-1209600` `726:-1209600` `727:-907200` `728:-907200` `729:-302400` `731:-3628800` `732:-1209600` `733:-907200` `734:-604800` `735:-80640` `736:-107520`|
|4|4|`737:-26880` `738:-725760` `739:-967680` `743:-725760` `745:-362880` `747:-241920` `748:-80640` `749:-1451520` `750:-1935360` `752:-483840` `754:-161280` `759:-362880`|
|4|4|`764:-161280` `765:-53760` `768:-483840` `770:-161280` `775:-1451520` `776:-483840` `779:-1451520` `780:-3628800` `782:-2903040` `783:-967680` `785:-2903040` `786:-2903040`|
|4|4|`787:-1935360` `793:-1451520` `795:-483840` `796:-483840` `798:-161280` `799:-1451520` `800:-1451520` `803:-2903040` `804:-725760` `805:-241920` `807:-725760` `809:-725760`|
|4|4|`811:-241920` `812:-241920` `813:-483840` `815:-161280` `829:-4354560` `832:-1451520` `833:-2903040` `834:-1451520` `838:-725760` `839:-241920` `840:-1451520` `841:-1451520`|
|4|4|`843:-2903040` `844:-2903040` `845:-967680` `849:-2903040` `850:-967680` `851:-725760` `854:-2903040` `855:-725760` `865:-6720` `866:-2240` `868:-60480` `870:-362880`|
|4|4|`874:-40320` `875:-80640` `876:-40320` `880:-483840` `883:-120960` `886:-60480` `887:-60480` `888:-20160` `890:-120960` `892:-241920` `898:-241920`|
|4|5|`21:302400` `24:907200` `26:302400` `28:544320` `31:1632960` `32:544320` `35:262080` `36:786240` `38:20160` `290:-1814400` `291:-1814400` `322:-5443200`|
|4|5|`323:-5443200` `330:-1814400` `331:-1814400` `411:-1209600` `426:-1209600` `440:-2419200` `441:-2419200` `443:-7257600` `444:-7257600` `450:-2419200` `458:-2419200` `459:-2419200`|
|4|5|`472:-2419200` `473:-2419200` `482:-9676800` `483:-7257600` `484:-10886400` `485:-7257600` `487:-2419200` `488:-2419200` `490:-2419200` `491:-2419200` `499:-7257600` `504:-7257600`|
|4|5|`512:-2419200` `519:-3628800` `533:-3628800` `534:-3628800` `536:-3628800` `537:-3628800` `539:-7257600` `541:-7257600` `544:-3628800` `545:-1209600` `546:-6048000` `547:-1209600`|
|4|5|`549:-2419200` `564:-302400` `565:-302400` `566:-907200` `582:-1814400` `583:-1209600` `584:-1814400` `585:-1209600` `588:-4233600` `589:-3628800` `590:-4536000` `591:-3628800`|
|4|5|`592:-1209600` `593:-1209600` `601:-1209600` `603:-1209600` `611:-907200` `612:-604800` `621:-302400` `624:-1209600` `625:-302400` `627:-604800` `634:-907200` `635:-302400`|
|4|5|`640:-604800` `643:-1814400` `652:-1814400` `655:-1209600` `656:-9072000` `657:-2419200` `658:-907200` `660:-9072000` `661:-6048000` `664:-7257600` `665:-7257600` `668:-6048000`|
|4|5|`669:-2419200` `672:-7257600` `674:-3628800` `677:-3628800` `680:-9072000` `681:-3628800` `683:-10886400` `684:-7257600` `686:-5443200` `687:-1209600` `688:-6048000` `689:-3628800`|
|4|5|`690:-5443200` `693:-2721600` `694:-3628800` `695:-1814400` `696:-604800` `697:-1209600` `702:-1814400` `703:-4838400` `704:-604800` `706:-1814400` `708:-1814400` `709:-2721600`|
|4|5|`711:-3628800` `712:-10886400` `713:-3628800` `714:-604800` `715:-1209600` `716:-604800` `718:-3628800` `720:-907200` `723:-3628800` `726:-604800` `727:-1814400` `730:-1814400`|
|4|5|`731:-3628800` `732:-604800` `740:-241920` `741:-1209600` `742:-483840` `744:-1451520` `746:-120960` `747:-120960` `749:-725760` `751:-241920` `753:-483840` `754:-241920`|
|4|5|`757:-362880` `758:-362880` `760:-362880` `761:-725760` `762:-362880` `766:-241920` `767:-241920` `774:-241920` `775:-725760` `776:-241920` `777:-241920` `779:-725760`|
|4|5|`781:-1935360` `782:-1451520` `784:-2903040` `785:-1451520` `786:-725760` `788:-725760` `789:-1451520` `790:-725760` `791:-2177280` `792:-1451520` `793:-725760` `794:-241920`|
|4|5|`795:-483840` `796:-725760` `797:-2419200` `798:-241920` `801:-1451520` `802:-2177280` `803:-725760` `808:-725760` `814:-241920` `828:-1451520` `830:-4354560` `831:-725760`|
|4|5|`832:-2177280` `833:-2903040` `834:-241920` `837:-725760` `840:-1451520` `842:-2903040` `843:-2903040` `846:-1451520` `847:-2177280` `848:-2903040` `849:-4354560` `850:-483840`|
|4|5|`852:-725760` `853:-1451520` `854:-1451520` `869:-60480` `872:-80640` `874:-40320` `875:-40320` `878:-241920` `879:-120960` `880:-241920` `882:-120960` `885:-60480`|
|4|5|`886:-60480` `889:-120960` `891:-241920` `893:-60480` `894:-120960` `895:-120960` `896:-120960` `897:-241920` `898:-120960`|
|4|6|`15:100800` `18:100800` `22:383040` `24:383040` `28:248640` `29:248640` `34:20160` `324:-604800` `325:-2419200` `332:-1814400` `445:-1612800` `446:-3225600`|
|4|6|`451:-2419200` `458:-806400` `459:-806400` `482:-3628800` `483:-2419200` `500:-3628800` `512:-806400` `536:-7257600` `537:-2419200` `567:-201600` `568:-403200` `575:-100800`|
|4|6|`576:-100800` `582:-806400` `583:-403200` `588:-1814400` `589:-1209600` `605:-907200` `606:-302400` `612:-1209600` `620:-403200` `632:-302400` `640:-604800` `648:-100800`|
|4|6|`667:-806400` `668:-6048000` `669:-806400` `679:-5443200` `680:-7257600` `681:-3628800` `694:-907200` `706:-907200` `712:-1814400` `716:-302400` `751:-322560` `752:-604800`|
|4|6|`753:-362880` `754:-241920` `760:-362880` `781:-967680` `782:-1451520` `783:-483840` `789:-483840` `796:-241920` `797:-483840` `798:-80640` `830:-1451520` `831:-725760`|
|4|6|`832:-725760` `833:-1451520` `834:-241920` `871:-6720` `872:-40320` `873:-40320` `874:-20160` `875:-40320` `876:-6720` `889:-120960` `890:-60480`|
|4|7|`27:44800` `29:403200` `32:403200` `33:44800` `34:67200` `36:604800` `37:604800` `38:13440` `39:120960` `113:-134400` `130:-1209600` `163:-1209600`|
|4|7|`185:-134400` `222:-403200` `225:-268800` `232:-1209600` `235:-403200` `236:-403200` `240:-403200` `252:-3628800` `256:-2419200` `267:-3628800` `276:-2419200` `278:-2419200`|
|4|7|`282:-2419200` `284:-3628800` `296:-2419200` `303:-1209600` `308:-268800` `320:-3628800` `329:-806400` `331:-1209600` `344:-1209600` `349:-2419200` `363:-403200` `376:-134400`|
|4|7|`377:-268800` `379:-44800` `381:-1209600` `385:-2419200` `387:-2419200` `391:-806400` `392:-672000` `394:-806400` `395:-537600` `400:-806400` `402:-403200` `403:-268800`|
|4|7|`404:-268800` `407:-2419200` `415:-537600` `416:-403200` `417:-806400` `418:-1209600` `425:-806400` `428:-403200` `429:-806400` `431:-403200` `432:-134400` `434:-44800`|
|4|7|`438:-4838400` `442:-9676800` `444:-4838400` `450:-1612800` `451:-1612800` `456:-3628800` `457:-2419200` `461:-4838400` `467:-2419200` `469:-1612800` `470:-806400` `480:-3628800`|
|4|7|`481:-4838400` `489:-806400` `490:-1209600` `491:-806400` `492:-2419200` `493:-2419200` `498:-4838400` `503:-4838400` `513:-806400` `514:-806400` `515:-1612800` `516:-2419200`|
|4|7|`518:-1612800` `520:-1209600` `523:-403200` `525:-806400` `529:-268800` `532:-403200` `534:-2419200` `542:-2419200` `547:-403200` `548:-806400` `549:-806400` `555:-806400`|
|4|7|`556:-403200` `559:-134400` `561:-67200` `563:-1209600` `565:-201600` `566:-604800` `568:-201600` `571:-201600` `572:-134400` `574:-201600` `576:-67200` `579:-201600`|
|4|7|`580:-1814400` `581:-2419200` `587:-3628800` `589:-1209600` `591:-1209600` `593:-201600` `594:-1209600` `595:-1209600` `600:-604800` `602:-806400` `606:-604800` `608:-1209600`|
|4|7|`609:-1209600` `614:-604800` `615:-604800` `617:-604800` `619:-604800` `620:-201600` `622:-134400` `628:-201600` `630:-201600` `633:-1209600` `638:-604800` `642:-1209600`|
|4|7|`653:-1209600` `659:-2419200` `661:-604800` `663:-2419200` `665:-4838400` `675:-1209600` `678:-4838400` `681:-1209600` `684:-2419200` `685:-4838400` `688:-1209600` `691:-604800`|
|4|7|`692:-1209600` `697:-604800` `701:-1209600` `707:-604800` `710:-1209600` `715:-2419200` `722:-1209600` `725:-604800` `736:-26880` `739:-483840` `742:-40320` `745:-120960`|
|4|7|`747:-80640` `750:-967680` `754:-80640` `758:-120960` `759:-241920` `765:-13440` `769:-40320` `770:-80640` `776:-241920` `783:-483840` `787:-967680` `795:-483840`|
|4|7|`805:-120960` `812:-241920`|
|4|8|`21:403200` `24:1209600` `26:403200` `28:739200` `31:2217600` `32:739200` `35:322560` `36:967680` `38:40320` `138:-1209600` `139:-403200` `172:-3628800`|
|4|8|`173:-1209600` `188:-1209600` `189:-403200` `238:-1209600` `241:-403200` `261:-3628800` `262:-1209600` `263:-2419200` `264:-806400` `270:-10886400` `271:-3628800` `278:-3628800`|
|4|8|`290:-3628800` `291:-1209600` `292:-10886400` `293:-2419200` `295:-7257600` `297:-2419200` `302:-1209600` `304:-2419200` `307:-806400` `310:-2419200` `311:-806400` `322:-10886400`|
|4|8|`323:-3628800` `330:-3628800` `331:-1209600` `343:-3628800` `344:-1209600` `345:-7257600` `348:-2419200` `353:-3628800` `354:-1209600` `357:-3628800` `358:-1209600` `360:-1209600`|
|4|8|`362:-403200` `367:-2419200` `368:-806400` `383:-1209600` `384:-403200` `387:-3628800` `397:-3225600` `398:-806400` `399:-2419200` `401:-806400` `402:-403200` `404:-134400`|
|4|8|`406:-268800` `409:-8467200` `410:-2419200` `411:-1209600` `413:-806400` `416:-2419200` `421:-1209600` `425:-1209600` `426:-403200` `428:-403200` `430:-806400` `432:-403200`|
|4|8|`433:-268800` `435:-134400` `440:-4838400` `441:-1612800` `443:-14515200` `444:-4838400` `450:-4838400` `456:-2419200` `458:-1209600` `459:-403200` `460:-14515200` `461:-4838400`|
|4|8|`463:-12096000` `466:-2419200` `468:-1612800` `469:-806400` `472:-2419200` `473:-4838400` `474:-806400` `476:-4838400` `477:-1612800` `479:-537600` `480:-1209600` `482:-6048000`|
|4|8|`483:-1209600` `484:-10886400` `485:-12096000` `486:-2419200` `487:-2419200` `488:-1612800` `490:-403200` `491:-3225600` `492:-7257600` `495:-19353600` `497:-4838400` `499:-7257600`|
|4|8|`503:-4838400` `504:-9676800` `505:-2419200` `507:-3628800` `508:-2419200` `510:-9676800` `511:-2419200` `512:-806400` `513:-3225600` `514:-806400` `515:-1612800` `517:-1612800`|
|4|8|`518:-4838400` `519:-3628800` `520:-2419200` `523:-806400` `524:-2419200` `525:-2419200` `526:-806400` `527:-134400` `528:-268800` `530:-537600` `533:-3628800` `534:-2419200`|
|4|8|`536:-3628800` `537:-1209600` `539:-7257600` `540:-1209600` `541:-4838400` `542:-2419200` `544:-3628800` `545:-1209600` `546:-4838400` `547:-1209600` `549:-806400` `552:-2419200`|
|4|8|`554:-806400` `555:-2419200` `556:-403200` `557:-403200` `558:-268800` `559:-134400` `564:-604800` `565:-201600` `566:-1814400` `573:-604800` `574:-201600` `577:-604800`|
|4|8|`578:-201600` `582:-604800` `583:-201600` `584:-1209600` `585:-1612800` `586:-403200` `588:-1814400` `589:-604800` `590:-3628800` `591:-4838400` `592:-1209600` `593:-604800`|
|4|8|`594:-3628800` `597:-1814400` `599:-604800` `601:-1209600` `603:-1612800` `604:-403200` `607:-7257600` `608:-4838400` `609:-604800` `610:-1814400` `613:-604800` `615:-1814400`|
|4|8|`618:-604800` `619:-1814400` `624:-604800` `625:-201600` `626:-1209600` `627:-1814400` `629:-806400` `631:-201600` `634:-1814400` `635:-1209600` `637:-1814400` `639:-3628800`|
|4|8|`641:-1209600` `643:-1814400` `644:-604800` `645:-1814400` `646:-604800` `647:-1814400` `656:-1209600` `657:-403200` `660:-3628800` `661:-2419200` `664:-7257600` `665:-2419200`|
|4|8|`668:-2419200` `669:-806400` `672:-3628800` `674:-1209600` `677:-7257600` `680:-3628800` `681:-1209600` `682:-7257600` `683:-7257600` `684:-4838400` `686:-3628800` `687:-1814400`|
|4|8|`688:-3628800` `689:-10886400` `690:-1814400` `691:-1209600` `692:-1209600` `696:-604800` `699:-3628800` `700:-1209600` `701:-1209600` `702:-604800` `703:-1209600` `704:-604800`|
|4|8|`708:-5443200` `711:-3628800` `712:-1814400` `713:-1814400` `714:-2419200` `715:-1209600` `716:-604800` `717:-604800` `718:-3628800` `719:-7257600` `721:-1209600` `723:-1814400`|
|4|8|`724:-2419200` `725:-604800` `726:-604800` `741:-241920` `742:-80640` `744:-725760` `746:-120960` `747:-40320` `749:-1088640` `753:-241920` `754:-80640` `756:-362880`|
|4|8|`757:-362880` `758:-241920` `760:-362880` `761:-725760` `762:-362880` `767:-483840` `769:-80640` `770:-80640` `777:-120960` `778:-40320` `779:-362880` `785:-2903040`|
|4|8|`790:-362880` `791:-725760` `794:-362880` `797:-483840` `798:-161280` `801:-1451520` `802:-725760` `803:-725760` `804:-725760` `808:-725760` `809:-725760` `811:-241920`|
|4|8|`814:-120960` `815:-40320` `816:-362880` `818:-725760` `820:-725760` `821:-362880`|
|4|9|`15:403200` `18:403200` `22:1008000` `24:1008000` `28:725760` `29:725760` `34:120960` `182:-1209600` `183:-1209600` `192:-1209600` `193:-1209600` `274:-6048000`|
|4|9|`275:-3628800` `281:-3628800` `282:-3628800` `290:-2419200` `294:-2419200` `304:-2419200` `306:-2419200` `312:-2419200` `322:-3628800` `324:-3628800` `325:-8467200` `332:-3628800`|
|4|9|`340:-3628800` `347:-2419200` `350:-1209600` `353:-7257600` `359:-1209600` `364:-2419200` `365:-1209600` `390:-1612800` `391:-1209600` `394:-1209600` `402:-403200` `404:-403200`|
|4|9|`411:-2419200` `412:-2419200` `413:-2419200` `417:-2419200` `425:-403200` `426:-403200` `428:-403200` `440:-1612800` `443:-2419200` `445:-7257600` `446:-10483200` `451:-4838400`|
|4|9|`456:-2419200` `458:-2419200` `459:-3225600` `465:-3628800` `469:-1612800` `470:-806400` `475:-4838400` `476:-6451200` `480:-1209600` `482:-8467200` `483:-6048000` `488:-2419200`|
|4|9|`490:-2419200` `496:-2419200` `499:-14515200` `500:-8467200` `501:-4838400` `502:-3628800` `503:-9676800` `504:-2419200` `508:-9676800` `512:-6451200` `514:-403200` `515:-6451200`|
|4|9|`516:-4032000` `519:-1209600` `534:-1209600` `536:-7257600` `537:-3628800` `541:-2419200` `546:-1209600` `547:-403200` `549:-806400` `567:-604800` `568:-1209600` `582:-1209600`|
|4|9|`583:-1209600` `588:-2419200` `589:-2419200` `593:-1209600` `601:-3024000` `602:-1209600` `605:-604800` `606:-604800` `610:-3628800` `611:-3628800` `612:-1209600` `613:-2419200`|
|4|9|`614:-1209600` `617:-1814400` `620:-604800` `626:-1209600` `627:-604800` `632:-604800` `633:-604800` `639:-3628800` `640:-604800` `641:-1209600` `642:-604800` `664:-3628800`|
|4|9|`667:-2419200` `668:-12096000` `669:-2419200` `672:-1814400` `677:-1814400` `679:-5443200` `680:-7257600` `681:-3628800` `683:-3628800` `686:-1814400` `687:-604800` `688:-1209600`|
|4|9|`689:-3628800` `690:-1814400` `702:-1209600` `703:-1209600` `706:-3628800` `707:-1209600` `708:-1814400` `711:-3628800` `712:-1814400` `713:-1814400` `714:-604800` `715:-1209600`|
|4|9|`751:-241920` `752:-604800` `753:-362880` `754:-241920` `760:-362880` `766:-120960` `768:-120960` `781:-483840` `782:-1451520` `783:-483840` `784:-725760` `788:-362880`|
|4|9|`789:-241920` `790:-362880` `791:-725760` `792:-362880` `793:-362880` `794:-120960` `795:-241920` `796:-483840` `799:-725760` `807:-725760` `813:-241920` `816:-362880`|
|4|9|`817:-362880`|
|4|10|`11:44800` `15:201600` `21:215040` `27:40320` `165:-134400` `187:-403200` `257:-268800` `269:-806400` `277:-2419200` `325:-1209600` `332:-1209600` `379:-44800`|
|4|10|`386:-806400` `392:-806400` `408:-1209600` `415:-134400` `445:-806400` `446:-1612800` `451:-806400` `459:-403200` `482:-2419200` `483:-3628800` `490:-403200` `567:-201600`|
|4|10|`568:-201600` `582:-403200` `583:-201600` `588:-1209600` `589:-604800` `593:-201600` `605:-604800` `606:-604800` `656:-604800` `660:-1814400` `661:-604800` `740:-40320`|
|4|10|`741:-120960` `742:-40320` `774:-120960` `775:-362880` `776:-120960`|
|4|11|`21:100800` `24:302400` `26:100800` `28:268800` `31:806400` `32:268800` `35:100800` `36:302400` `46:-100800` `59:-302400` `73:-100800` `114:-201600`|
|4|11|`122:-806400` `123:-403200` `128:-2419200` `134:-1814400` `135:-1209600` `138:-806400` `140:-403200` `142:-806400` `150:-403200` `170:-604800` `171:-1209600` `172:-2419200`|
|4|11|`180:-604800` `181:-604800` `188:-806400` `195:-201600` `206:-403200` `221:-604800` `223:-1008000` `224:-806400` `225:-100800` `226:-1814400` `230:-2721600` `234:-302400`|
|4|11|`235:-302400` `238:-403200` `239:-201600` `242:-806400` `246:-100800` `256:-604800` `258:-4838400` `260:-3628800` `261:-2419200` `263:-806400` `264:-1209600` `266:-1612800`|
|4|11|`270:-7257600` `276:-302400` `278:-4838400` `280:-6048000` `282:-1209600` `283:-907200` `290:-604800` `292:-3628800` `293:-3024000` `295:-2419200` `297:-2419200` `300:-604800`|
|4|11|`301:-3024000` `302:-604800` `303:-604800` `304:-806400` `305:-100800` `310:-806400` `311:-1612800` `313:-201600` `317:-1612800` `319:-403200` `322:-1814400` `330:-604800`|
|4|11|`331:-302400` `334:-302400` `343:-1209600` `345:-2419200` `348:-604800` `349:-604800` `352:-604800` `353:-1209600` `357:-1209600` `358:-302400` `360:-403200` `367:-806400`|
|4|11|`368:-201600` `370:-100800` `378:-403200` `380:-268800` `383:-806400` `387:-4838400` `389:-1209600` `393:-3225600` `394:-806400` `397:-806400` `398:-806400` `399:-806400`|
|4|11|`401:-806400` `404:-268800` `405:-537600` `406:-806400` `409:-2419200` `410:-4838400` `416:-1612800` `417:-1612800` `419:-1209600` `424:-1209600` `428:-806400` `430:-806400`|
|4|11|`431:-806400` `435:-537600` `437:-537600` `440:-806400` `443:-2419200` `450:-2419200` `452:-3225600` `460:-4838400` `463:-2419200` `468:-4838400` `469:-806400` `471:-1612800`|
|4|11|`473:-806400` `476:-1612800` `478:-268800` `479:-537600` `485:-2419200` `491:-2419200` `492:-2419200` `495:-4838400` `504:-2419200` `510:-2419200` `513:-1612800` `515:-1612800`|
|4|11|`517:-1612800` `518:-4838400` `524:-806400` `525:-806400` `526:-806400` `528:-268800` `530:-1075200` `531:-268800` `564:-100800` `566:-604800` `570:-302400` `573:-201600`|
|4|11|`578:-403200` `585:-403200` `591:-2419200` `592:-302400` `594:-1209600` `597:-604800` `603:-403200` `608:-2419200` `613:-604800` `615:-1209600` `618:-604800` `619:-1209600`|
|4|11|`625:-100800` `629:-403200` `631:-100800` `635:-604800` `641:-1209600` `644:-1209600` `647:-1209600` `651:-302400`|
|4|12|`15:302400` `18:302400` `22:907200` `24:907200` `28:571200` `29:571200` `34:100800` `64:-302400` `65:-100800` `79:-302400` `80:-100800` `136:-3628800`|
|4|12|`137:-806400` `138:-1209600` `140:-403200` `141:-403200` `145:-2419200` `146:-806400` `152:-1209600` `155:-403200` `172:-1814400` `177:-1814400` `180:-1209600` `182:-2419200`|
|4|12|`183:-1411200` `184:-201600` `192:-2419200` `193:-806400` `197:-1209600` `204:-604800` `205:-403200` `207:-201600` `228:-2116800` `229:-604800` `232:-1814400` `236:-302400`|
|4|12|`238:-2419200` `239:-806400` `240:-403200` `241:-100800` `244:-2419200` `245:-806400` `248:-100800` `261:-1209600` `265:-4838400` `266:-2016000` `270:-1814400` `274:-9676800`|
|4|12|`275:-3024000` `281:-7257600` `282:-2419200` `284:-1814400` `287:-2721600` `290:-2419200` `291:-604800` `294:-2419200` `295:-14515200` `296:-4233600` `297:-604800` `300:-3024000`|
|4|12|`302:-302400` `303:-604800` `304:-3225600` `305:-806400` `306:-3225600` `307:-1209600` `312:-3225600` `313:-806400` `314:-4838400` `316:-2016000` `318:-2419200` `319:-1008000`|
|4|12|`320:-907200` `322:-3628800` `323:-907200` `324:-1814400` `325:-3628800` `332:-604800` `336:-1814400` `340:-3628800` `342:-907200` `344:-604800` `347:-2419200` `348:-604800`|
|4|12|`349:-604800` `350:-1209600` `353:-2419200` `354:-302400` `358:-302400` `359:-403200` `362:-100800` `364:-806400` `365:-403200` `368:-201600` `390:-2419200` `391:-806400`|
|4|12|`394:-2419200` `399:-4838400` `400:-806400` `405:-537600` `406:-268800` `411:-2419200` `412:-2419200` `413:-3225600` `414:-806400` `417:-3225600` `418:-6048000` `421:-6048000`|
|4|12|`423:-1209600` `425:-806400` `426:-2016000` `427:-403200` `428:-806400` `429:-5644800` `431:-2016000` `432:-268800` `433:-268800` `435:-268800` `436:-537600` `437:-268800`|
|4|12|`445:-2419200` `446:-3225600` `451:-806400` `460:-4838400` `463:-2419200` `468:-1612800` `469:-806400` `475:-4838400` `476:-4838400` `478:-537600` `479:-537600` `482:-2419200`|
|4|12|`483:-2419200` `488:-806400` `490:-806400` `492:-2419200` `495:-4838400` `499:-14515200` `500:-2419200` `501:-4838400` `502:-2419200` `503:-4838400` `504:-2419200` `508:-4838400`|
|4|12|`510:-2419200` `512:-1612800` `513:-1612800` `515:-4838400` `516:-1612800` `518:-1612800` `523:-1612800` `524:-806400` `525:-806400` `527:-268800` `528:-806400` `529:-268800`|
|4|12|`530:-537600` `567:-302400` `568:-403200` `577:-302400` `579:-100800` `588:-1209600` `589:-1209600` `593:-403200` `601:-1209600` `602:-403200` `607:-1814400` `610:-2721600`|
|4|12|`611:-1814400` `612:-604800` `613:-1209600` `614:-604800` `617:-1209600` `620:-201600` `626:-1209600` `627:-604800` `628:-403200` `630:-201600` `632:-302400` `633:-302400`|
|4|12|`636:-100800` `637:-907200` `639:-5443200` `640:-604800` `641:-1209600` `642:-604800` `643:-907200` `645:-907200` `646:-1209600` `647:-604800` `648:-201600` `650:-302400`|
|4|13|`11:100800` `15:369600` `21:369600` `27:100800` `72:-100800` `84:-100800` `133:-403200` `144:-806400` `149:-1209600` `160:-604800` `164:-604800` `165:-806400`|
|4|13|`183:-604800` `187:-806400` `193:-604800` `225:-100800` `231:-1411200` `235:-705600` `243:-806400` `253:-907200` `256:-604800` `257:-806400` `268:-907200` `269:-3628800`|
|4|13|`274:-2419200` `275:-604800` `276:-302400` `277:-5241600` `281:-2419200` `282:-604800` `291:-302400` `304:-2419200` `309:-806400` `312:-2419200` `323:-907200` `325:-1209600`|
|4|13|`331:-302400` `332:-1209600` `386:-1612800` `390:-806400` `392:-1075200` `394:-806400` `408:-1612800` `411:-1209600` `412:-1209600` `413:-806400` `415:-1075200` `417:-806400`|
|4|13|`425:-1209600` `428:-1209600` `440:-806400` `443:-2419200` `445:-2419200` `446:-3225600` `450:-806400` `451:-806400` `473:-806400` `482:-2419200` `483:-2419200` `485:-2419200`|
|4|13|`488:-806400` `490:-806400` `491:-806400` `582:-604800` `583:-604800` `584:-302400` `588:-604800` `589:-604800` `590:-907200` `592:-302400` `605:-604800` `606:-604800`|
|4|13|`624:-302400` `632:-604800` `633:-604800` `634:-907200` `635:-302400`|
|4|14|`15:40320` `18:40320` `22:201600` `24:201600` `28:134400` `29:134400` `45:-282240` `46:-80640` `50:-201600` `52:-80640` `59:-120960` `62:-120960`|
|4|14|`64:-201600` `65:-40320` `76:-80640` `79:-201600` `83:-40320` `112:-443520` `113:-40320` `114:-403200` `115:-403200` `118:-40320` `119:-403200` `122:-161280`|
|4|14|`124:-806400` `128:-241920` `130:-241920` `132:-362880` `135:-2419200` `136:-2016000` `137:-80640` `138:-403200` `139:-80640` `140:-403200` `141:-403200` `145:-1612800`|
|4|14|`146:-161280` `151:-806400` `152:-403200` `154:-80640` `155:-403200` `156:-403200` `163:-120960` `167:-241920` `172:-604800` `173:-120960` `177:-604800` `179:-120960`|
|4|14|`182:-403200` `183:-201600` `189:-40320` `192:-403200` `197:-403200` `199:-80640` `204:-201600` `224:-1411200` `228:-1209600` `232:-3024000` `234:-1814400` `238:-806400`|
|4|14|`239:-806400` `240:-201600` `241:-604800` `244:-806400` `245:-2016000` `248:-201600` `249:-604800` `258:-1209600` `260:-604800` `265:-1612800` `266:-806400` `274:-1209600`|
|4|14|`278:-604800` `280:-1209600` `281:-1209600` `295:-4838400` `296:-604800` `297:-604800` `300:-604800` `301:-604800` `304:-806400` `306:-806400` `307:-201600` `311:-201600`|
|4|14|`312:-806400` `314:-1612800` `316:-403200` `317:-403200` `318:-806400` `319:-201600` `380:-134400` `390:-403200` `393:-806400` `394:-806400` `399:-806400` `405:-537600`|
|4|14|`406:-268800` `413:-806400` `416:-403200` `417:-1612800` `418:-1209600` `421:-1209600` `426:-403200` `429:-2419200` `430:-403200` `431:-806400` `432:-134400` `433:-268800`|
|4|14|`435:-268800` `436:-537600` `437:-537600`|
|4|15|`11:40320` `15:215040` `21:201600` `27:44800` `49:-80640` `53:-201600` `54:-94080` `63:-120960` `65:-40320` `69:-120960` `72:-201600` `80:-40320`|
|4|15|`84:-67200` `113:-40320` `117:-806400` `118:-147840` `120:-134400` `127:-362880` `130:-241920` `133:-403200` `136:-1209600` `137:-80640` `139:-40320` `140:-403200`|
|4|15|`144:-1612800` `145:-1209600` `146:-80640` `149:-1075200` `155:-403200` `157:-134400` `159:-362880` `160:-604800` `163:-120960` `164:-604800` `165:-403200` `173:-120960`|
|4|15|`183:-201600` `187:-134400` `189:-40320` `193:-201600` `228:-604800` `231:-2016000` `232:-604800` `235:-604800` `238:-1209600` `239:-201600` `240:-403200` `243:-1075200`|
|4|15|`244:-1209600` `245:-201600` `246:-336000` `248:-403200` `261:-604800` `264:-201600` `269:-1209600` `270:-1814400` `274:-2419200` `275:-604800` `277:-1209600` `278:-604800`|
|4|15|`281:-2419200` `282:-604800` `293:-604800` `304:-806400` `305:-201600` `309:-268800` `311:-201600` `312:-806400` `313:-201600` `386:-403200` `392:-134400` `397:-403200`|
|4|15|`404:-134400` `408:-806400` `409:-1209600` `411:-1209600` `412:-1209600` `413:-806400` `415:-403200` `416:-403200` `417:-806400` `425:-403200` `428:-806400` `432:-268800`|
|4|15|`434:-44800` `435:-403200`|
|4|16|`11:2240` `15:26880` `21:33600` `41:-2240` `42:-60480` `44:-20160` `45:-80640` `47:-13440` `49:-26880` `50:-80640` `53:-134400` `54:-26880`|
|4|16|`57:-20160` `63:-40320` `66:-6720` `69:-40320` `72:-33600` `112:-80640` `114:-201600` `115:-80640` `117:-403200` `118:-26880` `119:-201600` `120:-134400`|
|4|16|`122:-80640` `128:-241920` `136:-403200` `137:-80640` `142:-80640` `144:-268800` `145:-403200` `146:-80640` `149:-134400` `223:-100800` `230:-302400` `231:-201600`|
|4|16|`238:-403200` `239:-201600` `242:-100800` `243:-268800` `244:-403200` `245:-201600` `246:-33600`|
|5|5|`4:2520` `16:189000` `19:302400` `20:37800` `21:-30240` `23:1028160` `24:-30240` `25:1451520` `26:181440` `28:-20160` `30:952560` `31:1270080`|
|5|5|`35:241920` `326:-1360800` `327:-7257600` `328:-1587600` `333:-7257600` `334:-2721600` `447:-3628800` `448:-9676800` `449:-2116800` `452:-9676800` `472:-3628800` `473:-4838400`|
|5|5|`474:-604800` `484:-18144000` `485:-19353600` `486:-2419200` `487:-3628800` `491:-4838400` `506:-10886400` `519:-2419200` `520:-1209600` `521:-2419200` `522:-302400` `538:-14515200`|
|5|5|`539:-14515200` `540:-6048000` `541:-9676800` `542:-4838400` `543:-1209600` `545:-3628800` `547:-2419200` `548:-1209600` `569:-453600` `570:-1209600` `584:-3628800` `585:-2419200`|
|5|5|`586:-302400` `590:-9072000` `591:-9676800` `592:-1814400` `616:-1814400` `621:-907200` `622:-1209600` `623:-151200` `624:-1814400` `627:-3628800` `630:-1209600` `631:-604800`|
|5|5|`634:-907200` `643:-1814400` `645:-907200` `649:-453600` `650:-302400` `651:-302400` `652:-4536000` `653:-6048000` `654:-756000` `658:-4536000` `659:-6048000` `662:-907200`|
|5|5|`670:-1814400` `671:-10886400` `672:-10886400` `673:-2419200` `674:-2419200` `675:-1209600` `676:-302400` `682:-14515200` `683:-21772800` `684:-9676800` `685:-4838400` `687:-3628800`|
|5|5|`689:-14515200` `690:-7257600` `691:-2419200` `692:-4838400` `693:-1814400` `698:-3628800` `708:-3628800` `709:-1814400` `710:-1209600` `711:-3628800` `718:-7257600` `719:-3628800`|
|5|5|`720:-3628800` `721:-1209600` `722:-1209600` `727:-2721600` `728:-907200` `729:-302400` `735:-181440` `736:-241920` `737:-30240` `738:-907200` `739:-1209600` `743:-907200`|
|5|5|`755:-362880` `756:-1088640` `757:-362880` `758:-241920` `759:-120960` `761:-1451520` `766:-725760` `767:-1451520` `768:-725760` `769:-241920` `770:-483840` `784:-4354560`|
|5|5|`785:-5806080` `786:-2903040` `787:-1935360` `788:-362880` `790:-725760` `800:-725760` `804:-725760` `805:-120960` `806:-1451520` `807:-725760` `808:-1451520` `809:-725760`|
|5|5|`810:-362880` `811:-241920` `812:-241920` `828:-1451520` `829:-2903040` `836:-362880` `842:-8709120` `843:-8709120` `844:-2903040` `845:-967680` `847:-1451520` `851:-725760`|
|5|5|`869:-60480` `870:-120960` `877:-241920` `884:-15120` `885:-241920` `886:-181440` `887:-60480` `888:-20160` `891:-483840` `892:-241920`|
|5|6|`17:302400` `19:151200` `23:544320` `24:272160` `28:262080` `322:-1814400` `323:-907200` `330:-1814400` `331:-907200` `443:-2419200` `444:-1209600` `450:-2419200`|
|5|6|`456:-2419200` `457:-1209600` `480:-2419200` `481:-1209600` `499:-3628800` `503:-2419200` `504:-2419200` `505:-604800` `533:-3628800` `534:-4838400` `535:-604800` `566:-302400`|
|5|6|`580:-1209600` `581:-604800` `587:-1209600` `610:-1814400` `611:-907200` `615:-1209600` `639:-1814400` `643:-907200` `644:-604800` `664:-3628800` `665:-4838400` `666:-604800`|
|5|6|`677:-5443200` `678:-3628800` `693:-907200` `708:-1814400` `709:-907200` `710:-604800` `749:-1088640` `750:-725760` `780:-725760` `784:-725760` `785:-1451520` `786:-725760`|
|5|6|`787:-483840` `828:-725760` `829:-1451520` `869:-60480` `870:-120960`|
|5|7|`21:201600` `24:604800` `26:201600` `28:302400` `31:907200` `32:302400` `35:60480` `36:181440` `133:-604800` `164:-1814400` `186:-604800` `231:-604800`|
|5|7|`253:-1814400` `257:-1209600` `268:-5443200` `276:-1814400` `277:-4838400` `281:-3628800` `290:-1814400` `294:-1814400` `309:-1209600` `322:-5443200` `330:-3628800` `332:-1814400`|
|5|7|`347:-3628800` `350:-1814400` `361:-604800` `365:-1209600` `382:-604800` `385:-1814400` `386:-1612800` `390:-1209600` `408:-4032000` `411:-604800` `412:-1209600` `415:-1209600`|
|5|7|`440:-2419200` `443:-14515200` `446:-4838400` `450:-2419200` `458:-1209600` `459:-2419200` `465:-3628800` `472:-604800` `480:-3628800` `482:-4838400` `483:-9676800` `484:-3628800`|
|5|7|`487:-1814400` `488:-2419200` `490:-1209600` `496:-7257600` `499:-3628800` `501:-2419200` `502:-4838400` `512:-2419200` `514:-1209600` `516:-2419200` `527:-403200` `533:-1814400`|
|5|7|`536:-1814400` `537:-1209600` `539:-3628800` `544:-1814400` `545:-1814400` `546:-3628800` `557:-403200` `564:-604800` `566:-907200` `567:-302400` `580:-1814400` `582:-604800`|
|5|7|`583:-1209600` `584:-302400` `588:-4233600` `589:-1814400` `590:-1814400` `592:-302400` `593:-604800` `598:-907200` `601:-604800` `606:-907200` `609:-1814400` `611:-907200`|
|5|7|`612:-1209600` `614:-907200` `632:-907200` `640:-1814400` `648:-302400` `656:-604800` `658:-907200` `660:-3628800` `661:-604800` `664:-3628800` `668:-1209600` `672:-1814400`|
|5|7|`677:-1814400` `680:-5443200` `683:-7257600` `686:-3628800` `687:-604800` `688:-1209600` `693:-907200` `694:-907200` `706:-907200` `709:-907200` `712:-3628800` `713:-1814400`|
|5|7|`716:-302400` `720:-907200` `723:-1814400` `726:-604800` `738:-181440` `741:-241920` `744:-362880` `746:-60480` `749:-362880` `753:-241920` `757:-362880` `760:-181440`|
|5|7|`761:-362880` `762:-181440` `775:-362880` `782:-725760` `786:-725760` `792:-725760` `793:-362880` `802:-725760` `810:-181440` `817:-362880` `819:-181440` `821:-362880`|
|5|8|`16:604800` `19:806400` `20:100800` `23:1512000` `25:2016000` `26:252000` `30:1088640` `31:1451520` `35:181440` `174:-1814400` `175:-1814400` `190:-2419200`|
|5|8|`191:-2419200` `194:-302400` `195:-302400` `272:-9072000` `273:-5443200` `279:-7257600` `280:-7257600` `283:-907200` `290:-3628800` `294:-3628800` `299:-3628800` `304:-4838400`|
|5|8|`305:-604800` `312:-4838400` `313:-604800` `322:-3628800` `326:-5443200` `327:-16329600` `328:-907200` `333:-7257600` `334:-907200` `339:-7257600` `346:-7257600` `347:-3628800`|
|5|8|`353:-2419200` `355:-1209600` `356:-6048000` `357:-2419200` `358:-604800` `360:-1209600` `364:-2419200` `366:-2419200` `367:-2419200` `368:-604800` `388:-2419200` `389:-1814400`|
|5|8|`393:-2419200` `411:-3628800` `412:-3628800` `413:-4838400` `414:-604800` `417:-4838400` `422:-604800` `440:-2419200` `443:-2419200` `447:-10886400` `448:-19353600` `449:-1209600`|
|5|8|`452:-9676800` `462:-7257600` `464:-9676800` `465:-2419200` `472:-5443200` `473:-9676800` `474:-907200` `476:-4838400` `478:-806400` `479:-1612800` `484:-16329600` `485:-14515200`|
|5|8|`486:-1512000` `487:-4233600` `491:-4838400` `494:-14515200` `495:-14515200` `496:-2419200` `499:-7257600` `501:-7257600` `504:-2419200` `506:-18144000` `507:-2419200` `508:-9676800`|
|5|8|`509:-19353600` `510:-7257600` `511:-1209600` `515:-4838400` `517:-4838400` `518:-4838400` `519:-7257600` `520:-2419200` `521:-1814400` `522:-604800` `523:-806400` `524:-2419200`|
|5|8|`525:-806400` `526:-806400` `538:-7257600` `539:-7257600` `540:-1209600` `541:-9676800` `542:-2419200` `543:-604800` `545:-1209600` `551:-2419200` `552:-2419200` `554:-806400`|
|5|8|`555:-1612800` `569:-907200` `570:-2116800` `584:-2721600` `585:-3628800` `586:-453600` `590:-4536000` `591:-6048000` `592:-1814400` `596:-2721600` `597:-1814400` `601:-1814400`|
|5|8|`607:-3628800` `611:-3628800` `613:-1814400` `616:-2721600` `617:-1209600` `618:-3024000` `619:-1209600` `626:-3628800` `627:-1814400` `628:-604800` `629:-1209600` `639:-3628800`|
|5|8|`643:-1814400` `644:-1209600` `652:-1814400` `653:-2419200` `654:-302400` `658:-1814400` `659:-2419200` `664:-3628800` `670:-3628800` `671:-18144000` `672:-9072000` `673:-1209600`|
|5|8|`674:-4838400` `675:-2419200` `676:-604800` `682:-7257600` `683:-14515200` `684:-4838400` `685:-2419200` `689:-10886400` `690:-1814400` `691:-604800` `692:-1209600` `693:-3628800`|
|5|8|`698:-7257600` `699:-3628800` `700:-1209600` `701:-2419200` `708:-3628800` `709:-1814400` `710:-1209600` `711:-3628800` `718:-7257600` `719:-3628800` `720:-907200` `721:-1209600`|
|5|8|`722:-1209600` `738:-362880` `739:-483840` `743:-362880` `755:-544320` `756:-1088640` `757:-907200` `758:-483840` `759:-241920` `761:-725760` `784:-2177280` `785:-1451520`|
|5|8|`786:-725760` `787:-483840` `790:-725760` `791:-725760` `800:-725760` `806:-725760` `807:-725760` `808:-2177280` `809:-725760` `810:-181440` `811:-241920` `812:-241920`|
|5|8|`818:-725760` `819:-181440`|
|5|9|`12:201600` `14:100800` `17:1108800` `19:554400` `23:1471680` `24:735840` `28:544320` `174:-604800` `175:-302400` `190:-2419200` `191:-1209600` `194:-604800`|
|5|9|`195:-302400` `263:-1209600` `264:-604800` `272:-3628800` `273:-907200` `279:-7257600` `280:-3628800` `283:-1814400` `290:-1814400` `291:-907200` `292:-6048000` `293:-2419200`|
|5|9|`294:-604800` `299:-1209600` `302:-302400` `303:-302400` `310:-1209600` `311:-604800` `322:-7257600` `323:-3628800` `330:-1814400` `331:-907200` `339:-7257600` `345:-4838400`|
|5|9|`346:-4838400` `347:-2419200` `348:-1209600` `349:-1209600` `355:-1814400` `360:-2419200` `361:-604800` `362:-302400` `363:-302400` `388:-1209600` `389:-302400` `393:-2419200`|
|5|9|`397:-2419200` `398:-604800` `402:-1209600` `403:-604800` `404:-403200` `406:-806400` `409:-6048000` `410:-2419200` `411:-1814400` `412:-604800` `416:-1209600` `422:-604800`|
|5|9|`425:-1209600` `426:-1209600` `427:-302400` `428:-403200` `430:-403200` `431:-403200` `440:-2419200` `441:-1209600` `443:-9676800` `444:-4838400` `450:-2419200` `456:-6048000`|
|5|9|`457:-3024000` `462:-3628800` `463:-7257600` `464:-4838400` `465:-1209600` `466:-604800` `467:-604800` `468:-3225600` `469:-1612800` `480:-6048000` `481:-3024000` `489:-1209600`|
|5|9|`494:-14515200` `495:-24192000` `496:-4838400` `497:-2419200` `498:-2419200` `499:-7257600` `501:-2419200` `503:-4838400` `504:-4838400` `505:-1209600` `507:-4838400` `509:-2419200`|
|5|9|`513:-3225600` `514:-1612800` `515:-1612800` `517:-1612800` `518:-1612800` `521:-604800` `533:-1814400` `534:-2419200` `535:-302400` `540:-1209600` `545:-1209600` `547:-806400`|
|5|9|`548:-403200` `564:-302400` `565:-151200` `566:-1209600` `571:-604800` `572:-302400` `577:-907200` `578:-604800` `579:-302400` `580:-3024000` `581:-1512000` `587:-3024000`|
|5|9|`596:-1814400` `597:-1209600` `598:-302400` `599:-151200` `600:-151200` `601:-1814400` `602:-1209600` `603:-1209600` `604:-302400` `607:-9072000` `608:-4838400` `609:-2419200`|
|5|9|`610:-907200` `611:-907200` `615:-604800` `616:-907200` `617:-604800` `618:-604800` `619:-604800` `637:-907200` `638:-302400` `664:-7257600` `665:-9676800` `666:-1209600`|
|5|9|`671:-3628800` `677:-5443200` `678:-3628800` `682:-7257600` `683:-7257600` `684:-4838400` `685:-2419200` `693:-907200` `705:-907200` `708:-1814400` `709:-907200` `710:-604800`|
|5|9|`746:-181440` `747:-241920` `748:-30240` `749:-1088640` `750:-725760` `755:-181440` `756:-362880` `757:-362880` `758:-241920` `759:-120960` `780:-1451520` `784:-725760`|
|5|9|`785:-1451520` `786:-725760` `787:-483840`|
|5|10|`12:100800` `16:252000` `21:181440` `164:-302400` `186:-302400` `256:-604800` `268:-907200` `276:-1512000` `281:-1209600` `282:-604800` `330:-1209600` `331:-604800`|
|5|10|`385:-907200` `394:-806400` `407:-604800` `412:-604800` `417:-806400` `443:-2419200` `444:-1209600` `450:-806400` `480:-2419200` `481:-1209600` `487:-604800` `491:-806400`|
|5|10|`566:-604800` `580:-604800` `581:-302400` `587:-604800` `590:-907200` `591:-1209600` `658:-907200` `659:-1209600` `738:-181440` `739:-241920`|
|5|11|`16:226800` `19:302400` `20:37800` `23:604800` `25:806400` `26:100800` `30:226800` `31:302400` `60:-226800` `74:-302400` `85:-37800` `129:-2721600`|
|5|11|`133:-907200` `143:-2419200` `149:-1209600` `150:-302400` `157:-151200` `164:-907200` `169:-1814400` `174:-1814400` `175:-907200` `190:-2419200` `193:-604800` `194:-604800`|
|5|11|`203:-604800` `204:-604800` `207:-151200` `227:-1587600` `230:-1814400` `231:-1814400` `243:-2419200` `246:-302400` `253:-907200` `259:-3628800` `268:-907200` `272:-7257600`|
|5|11|`273:-1814400` `275:-1814400` `279:-14515200` `281:-3628800` `283:-907200` `290:-1814400` `291:-907200` `294:-1814400` `299:-3628800` `304:-2419200` `305:-604800` `307:-604800`|
|5|11|`310:-1209600` `312:-4838400` `313:-302400` `315:-2419200` `318:-2419200` `319:-604800` `322:-1814400` `323:-1814400` `326:-1360800` `327:-3628800` `328:-680400` `339:-3628800`|
|5|11|`341:-907200` `342:-907200` `353:-1209600` `354:-302400` `355:-604800` `356:-1209600` `357:-1209600` `358:-907200` `388:-4233600` `390:-1209600` `393:-2419200` `409:-3628800`|
|5|11|`411:-1814400` `412:-3628800` `413:-4838400` `414:-302400` `417:-2419200` `422:-1814400` `426:-1209600` `428:-1209600` `430:-1209600` `431:-1209600` `432:-403200` `433:-806400`|
|5|11|`443:-2419200` `447:-1814400` `448:-4838400` `449:-302400` `462:-3628800` `464:-2419200` `465:-1209600` `484:-1814400` `485:-4838400` `486:-302400` `487:-604800` `494:-7257600`|
|5|11|`499:-3628800` `504:-2419200` `506:-3628800` `507:-1209600` `508:-4838400` `509:-4838400` `510:-4838400` `511:-604800` `519:-1209600` `521:-604800` `523:-403200` `524:-1612800`|
|5|11|`525:-806400` `526:-806400` `569:-453600` `570:-302400` `590:-2721600` `591:-1209600` `596:-907200` `607:-1814400` `611:-1814400` `616:-1814400` `617:-604800` `618:-604800`|
|5|11|`619:-604800` `634:-907200` `639:-1814400` `643:-2721600` `644:-604800` `645:-907200` `649:-226800` `650:-302400` `651:-302400`|
|5|12|`12:302400` `14:151200` `17:1108800` `19:554400` `23:1108800` `24:554400` `28:302400` `77:-302400` `81:-453600` `86:-151200` `138:-1209600` `140:-604800`|
|5|12|`147:-2419200` `148:-3628800` `151:-1209600` `152:-1209600` `155:-604800` `161:-1814400` `172:-1209600` `174:-2419200` `175:-1814400` `178:-3024000` `180:-302400` `181:-302400`|
|5|12|`190:-2419200` `191:-1209600` `197:-1209600` `201:-1209600` `203:-604800` `205:-302400` `206:-302400` `233:-1814400` `234:-2721600` `238:-2419200` `239:-1209600` `244:-2419200`|
|5|12|`245:-1209600` `254:-907200` `261:-1814400` `263:-2419200` `264:-1814400` `266:-1209600` `270:-1814400` `272:-10886400` `273:-4536000` `279:-7257600` `280:-3628800` `286:-7257600`|
|5|12|`287:-1814400` `288:-1814400` `290:-1814400` `291:-907200` `292:-7257600` `293:-3024000` `294:-1209600` `295:-4838400` `297:-604800` `298:-12096000` `299:-4838400` `300:-2419200`|
|5|12|`301:-1814400` `302:-1209600` `303:-604800` `310:-2419200` `311:-1209600` `314:-9676800` `315:-2419200` `316:-1209600` `317:-1209600` `322:-1814400` `323:-907200` `336:-1814400`|
|5|12|`339:-3628800` `341:-1814400` `345:-2419200` `346:-2419200` `347:-1209600` `348:-1814400` `349:-604800` `352:-604800` `388:-2419200` `389:-1209600` `393:-2419200` `397:-2419200`|
|5|12|`398:-1209600` `399:-2419200` `405:-1612800` `406:-806400` `409:-4838400` `410:-2419200` `416:-2419200` `420:-3628800` `421:-10886400` `422:-1209600` `423:-604800` `424:-604800`|
|5|12|`429:-4032000` `430:-2016000` `440:-2419200` `441:-1209600` `443:-2419200` `444:-1209600` `453:-2419200` `456:-2419200` `457:-1209600` `460:-4838400` `462:-7257600` `463:-7257600`|
|5|12|`464:-12096000` `465:-2419200` `466:-1209600` `467:-1209600` `468:-4838400` `469:-806400` `471:-1612800` `480:-2419200` `481:-1209600` `494:-7257600` `495:-14515200` `496:-2419200`|
|5|12|`497:-1209600` `498:-1209600` `501:-2419200` `507:-2419200` `509:-2419200` `513:-1612800` `514:-806400` `515:-1612800` `517:-1612800` `518:-1612800` `564:-302400` `565:-151200`|
|5|12|`566:-302400` `580:-1209600` `581:-604800` `587:-1209600` `596:-2721600` `597:-2419200` `598:-604800` `599:-302400` `600:-302400` `607:-5443200` `608:-2419200` `609:-1209600`|
|5|12|`611:-1814400` `616:-907200` `617:-604800` `618:-604800` `619:-604800` `637:-1814400` `638:-604800` `639:-1814400` `643:-907200` `644:-604800`|
|5|13|`8:37800` `12:252000` `16:478800` `21:252000` `60:-37800` `74:-151200` `85:-37800` `123:-151200` `129:-604800` `133:-302400` `134:-756000` `143:-1209600`|
|5|13|`150:-453600` `164:-1209600` `169:-1209600` `171:-604800` `186:-302400` `194:-302400` `195:-151200` `223:-604800` `225:-302400` `227:-453600` `230:-2419200` `231:-604800`|
|5|13|`235:-302400` `242:-302400` `253:-907200` `256:-1512000` `259:-1209600` `260:-1814400` `268:-3628800` `273:-907200` `276:-2419200` `279:-4838400` `280:-3628800` `281:-2419200`|
|5|13|`282:-1209600` `283:-302400` `308:-302400` `310:-1209600` `311:-604800` `328:-226800` `330:-604800` `331:-302400` `334:-302400` `377:-302400` `380:-403200` `382:-302400`|
|5|13|`385:-2721600` `388:-1209600` `389:-302400` `390:-1209600` `391:-604800` `393:-2419200` `394:-403200` `407:-1512000` `409:-2419200` `410:-1209600` `416:-403200` `443:-4838400`|
|5|13|`444:-2419200` `448:-2419200` `450:-1612800` `452:-3225600` `480:-1209600` `481:-604800` `487:-604800` `489:-403200` `491:-806400` `564:-302400` `565:-151200` `566:-302400`|
|5|13|`569:-226800` `570:-302400` `580:-1209600` `581:-604800` `587:-1209600` `590:-907200` `591:-1209600`|
|5|14|`12:60480` `14:30240` `17:302400` `19:151200` `23:201600` `24:100800` `49:-120960` `51:-302400` `54:-60480` `55:-151200` `58:-181440` `69:-120960`|
|5|14|`71:-120960` `77:-302400` `78:-60480` `80:-30240` `81:-151200` `83:-30240` `116:-604800` `117:-604800` `119:-302400` `120:-302400` `126:-181440` `127:-362880`|
|5|14|`131:-1814400` `137:-241920` `138:-604800` `139:-120960` `140:-302400` `145:-1209600` `147:-2419200` `148:-2419200` `152:-604800` `153:-604800` `155:-604800` `156:-302400`|
|5|14|`159:-362880` `161:-907200` `162:-362880` `172:-604800` `173:-241920` `174:-604800` `175:-604800` `178:-604800` `179:-120960` `180:-151200` `181:-151200` `228:-907200`|
|5|14|`233:-3628800` `234:-907200` `238:-1209600` `239:-1209600` `241:-302400` `244:-3628800` `245:-604800` `247:-604800` `248:-302400` `249:-302400` `254:-907200` `270:-1814400`|
|5|14|`272:-1814400` `273:-1814400` `286:-3628800` `287:-907200` `288:-907200` `292:-1209600` `293:-1209600` `295:-2419200` `297:-604800` `298:-2419200` `299:-1209600` `300:-1209600`|
|5|14|`301:-1209600` `302:-302400` `388:-1209600` `389:-302400` `409:-2419200` `410:-604800` `412:-604800` `420:-1814400` `421:-3628800` `422:-1209600` `423:-302400` `424:-302400`|
|5|14|`428:-403200` `429:-806400` `430:-806400` `431:-403200`|
|5|15|`8:30240` `12:181440` `16:252000` `21:100800` `46:-60480` `51:-151200` `52:-60480` `55:-151200` `58:-120960` `59:-60480` `60:-151200` `71:-120960`|
|5|15|`74:-151200` `76:-60480` `114:-302400` `116:-302400` `119:-604800` `122:-241920` `123:-302400` `126:-181440` `128:-241920` `129:-1814400` `131:-1209600` `132:-241920`|
|5|15|`133:-302400` `134:-907200` `135:-604800` `143:-1209600` `148:-1209600` `150:-302400` `151:-604800` `162:-181440` `164:-302400` `167:-241920` `169:-604800` `171:-302400`|
|5|15|`223:-604800` `224:-604800` `227:-907200` `230:-2116800` `233:-1209600` `234:-1209600` `242:-604800` `244:-1209600` `245:-604800` `253:-907200` `254:-907200` `256:-604800`|
|5|15|`258:-1209600` `259:-2419200` `260:-1814400` `268:-907200` `273:-907200` `276:-604800` `279:-2419200` `280:-2419200` `283:-604800` `382:-302400` `385:-907200` `388:-1814400`|
|5|15|`389:-604800` `393:-1612800` `407:-604800` `409:-1209600` `410:-604800` `412:-604800` `416:-806400` `417:-806400`|
|5|16|`8:2520` `12:30240` `16:37800` `42:-30240` `43:-37800` `44:-20160` `46:-30240` `48:-120960` `50:-60480` `51:-151200` `52:-30240` `56:-15120`|
|5|16|`57:-20160` `58:-60480` `59:-30240` `60:-37800` `112:-120960` `114:-151200` `116:-604800` `119:-302400` `126:-181440` `128:-241920` `129:-302400` `131:-604800`|
|5|16|`132:-120960` `134:-151200` `135:-302400` `227:-453600` `230:-604800` `233:-604800` `234:-604800`|
|6|6|`3:672` `15:-10080` `17:151200` `18:30240` `22:151200` `28:-6720` `320:-907200` `321:-181440` `329:-302400` `442:-1209600` `492:-2419200` `493:-604800`|
|6|6|`532:-806400` `637:-907200` `638:-302400` `663:-2419200` `705:-302400` `780:-725760`|
|6|7|`15:134400` `18:134400` `22:201600` `24:201600` `28:40320` `29:40320` `165:-403200` `187:-403200` `269:-2016000` `277:-1209600` `324:-1209600` `325:-3628800`|
|6|7|`359:-806400` `386:-537600` `392:-403200` `445:-4032000` `446:-2419200` `458:-403200` `459:-806400` `482:-3628800` `488:-403200` `500:-3628800` `512:-3225600` `527:-134400`|
|6|7|`567:-201600` `568:-201600` `582:-806400` `583:-201600` `605:-604800` `612:-604800` `620:-403200` `656:-604800` `667:-806400` `668:-3628800` `669:-806400` `702:-604800`|
|6|7|`703:-1209600` `704:-201600` `740:-40320` `751:-80640` `752:-120960` `753:-120960` `754:-80640` `774:-120960` `789:-241920` `813:-241920` `814:-120960` `815:-40320`|
|6|8|`12:134400` `14:67200` `17:604800` `19:302400` `23:645120` `24:322560` `28:120960` `182:-403200` `183:-201600` `192:-1209600` `193:-604800` `274:-2419200`|
|6|8|`275:-604800` `281:-3628800` `282:-1814400` `290:-1209600` `291:-604800` `304:-806400` `306:-806400` `307:-201600` `322:-3628800` `323:-1814400` `340:-3628800` `353:-4838400`|
|6|8|`354:-604800` `390:-806400` `391:-201600` `394:-1209600` `402:-403200` `403:-201600` `411:-1209600` `413:-806400` `425:-806400` `426:-1209600` `427:-201600` `432:-134400`|
|6|8|`433:-268800` `440:-1612800` `441:-806400` `443:-4838400` `444:-2419200` `456:-1209600` `457:-604800` `465:-1209600` `475:-2419200` `476:-3225600` `477:-403200` `499:-10886400`|
|6|8|`503:-7257600` `504:-7257600` `505:-1814400` `508:-2419200` `519:-1209600` `523:-403200` `524:-806400` `564:-201600` `565:-100800` `566:-604800` `571:-201600` `572:-100800`|
|6|8|`580:-604800` `581:-302400` `601:-2419200` `602:-1209600` `603:-1209600` `604:-302400` `610:-1814400` `615:-1209600` `626:-1209600` `627:-604800` `628:-201600` `629:-403200`|
|6|8|`664:-3628800` `665:-4838400` `666:-604800` `672:-1814400` `689:-3628800` `690:-1814400` `691:-604800` `692:-1209600` `746:-120960` `747:-161280` `748:-20160` `749:-362880`|
|6|8|`750:-241920` `766:-120960` `767:-241920` `768:-120960` `769:-40320` `770:-80640` `788:-362880` `804:-725760` `805:-120960`|
|6|9|`13:201600` `14:40320` `17:504000` `18:100800` `22:362880` `172:-604800` `173:-120960` `188:-604800` `189:-120960` `270:-1814400` `271:-362880` `278:-1814400`|
|6|9|`284:-1814400` `285:-362880` `295:-2419200` `296:-604800` `297:-604800` `320:-1814400` `321:-362880` `343:-2419200` `344:-1209600` `387:-604800` `418:-2419200` `419:-604800`|
|6|9|`421:-1209600` `429:-806400` `430:-403200` `438:-2419200` `439:-483840` `442:-2419200` `460:-4838400` `461:-2419200` `492:-7257600` `493:-1814400` `507:-1209600` `513:-806400`|
|6|9|`514:-403200` `563:-302400` `594:-3628800` `595:-907200` `607:-1814400` `608:-1209600` `609:-604800` `663:-2419200` `677:-1814400` `678:-1209600` `745:-241920` `749:-362880`|
|6|9|`750:-241920`|
|6|10|`12:40320` `15:73920` `163:-120960` `185:-40320` `267:-362880` `278:-604800` `329:-201600` `416:-403200` `442:-1209600` `489:-134400` `587:-604800`|
|6|11|`12:100800` `14:50400` `17:268800` `19:134400` `23:100800` `24:50400` `72:-100800` `84:-50400` `144:-806400` `149:-403200` `160:-604800` `182:-806400`|
|6|11|`183:-1209600` `184:-100800` `231:-604800` `235:-302400` `253:-302400` `274:-6048000` `275:-1814400` `290:-604800` `291:-604800` `304:-2419200` `305:-100800` `306:-2419200`|
|6|11|`307:-806400` `390:-806400` `391:-403200` `411:-1209600` `425:-806400` `426:-806400` `427:-201600` `440:-1612800` `441:-403200` `473:-806400` `475:-2419200` `476:-4838400`|
|6|11|`477:-403200` `478:-268800` `479:-537600` `564:-100800` `565:-50400` `584:-302400` `601:-604800` `602:-403200` `603:-403200` `604:-100800` `624:-302400` `627:-604800`|
|6|11|`630:-201600` `631:-100800`|
|6|12|`9:50400` `10:10080` `13:285600` `14:57120` `17:453600` `18:90720` `22:151200` `64:-50400` `65:-10080` `79:-151200` `80:-30240` `136:-806400`|
|6|12|`137:-80640` `138:-403200` `139:-80640` `140:-100800` `141:-100800` `145:-1209600` `146:-241920` `172:-1209600` `173:-241920` `177:-1209600` `180:-604800` `228:-604800`|
|6|12|`229:-60480` `232:-907200` `236:-302400` `237:-60480` `238:-806400` `239:-201600` `240:-201600` `241:-302400` `261:-1209600` `262:-241920` `265:-1612800` `266:-806400`|
|6|12|`270:-3628800` `271:-725760` `284:-907200` `285:-181440` `287:-907200` `295:-7257600` `296:-1814400` `297:-1814400` `300:-604800` `302:-302400` `383:-403200` `384:-80640`|
|6|12|`387:-1209600` `395:-403200` `396:-80640` `399:-3225600` `400:-604800` `401:-604800` `405:-537600` `406:-268800` `418:-2419200` `419:-604800` `438:-1209600` `439:-241920`|
|6|12|`460:-4838400` `461:-2419200` `463:-2419200` `468:-1612800` `469:-806400` `561:-50400` `562:-10080` `563:-151200` `573:-403200` `574:-201600` `577:-302400` `578:-201600`|
|6|12|`579:-100800` `594:-1209600` `595:-302400` `610:-907200` `615:-604800`|
|6|13|`9:30240` `12:110880` `15:110880` `59:-30240` `73:-30240` `128:-241920` `130:-241920` `135:-302400` `142:-241920` `163:-241920` `170:-302400` `226:-181440`|
|6|13|`232:-604800` `234:-302400` `252:-725760` `258:-1209600` `267:-725760` `278:-1814400` `283:-302400` `381:-241920` `387:-1814400` `393:-806400` `442:-1209600` `450:-806400`|
|6|13|`563:-302400` `566:-302400`|
|6|14|`9:20160` `10:4032` `13:100800` `14:20160` `17:67200` `18:13440` `53:-100800` `54:-20160` `63:-80640` `64:-100800` `65:-60480` `117:-201600`|
|6|14|`118:-40320` `127:-120960` `136:-2016000` `137:-241920` `138:-201600` `139:-80640` `140:-302400` `141:-302400` `228:-604800` `229:-120960` `238:-806400` `240:-201600`|
|6|14|`241:-201600` `261:-1209600` `262:-120960` `264:-201600` `265:-1612800` `266:-1209600` `383:-201600` `384:-40320` `397:-403200` `399:-806400` `400:-201600` `401:-201600`|
|6|14|`404:-134400` `406:-268800`|
|6|15|`7:4032` `9:32256` `12:73920` `15:40320` `45:-40320` `46:-20160` `50:-60480` `59:-60480` `62:-60480` `112:-80640` `113:-40320` `114:-100800`|
|6|15|`115:-120960` `122:-161280` `124:-201600` `128:-483840` `130:-120960` `132:-120960` `135:-907200` `221:-120960` `222:-120960` `224:-806400` `226:-362880` `232:-604800`|
|6|15|`252:-362880` `258:-1209600` `260:-604800` `376:-40320` `378:-201600` `380:-134400` `381:-120960` `387:-604800` `394:-403200`|
|6|16|`7:672` `9:8064` `12:10080` `42:-10080` `44:-6720` `45:-100800` `46:-10080` `112:-80640` `114:-100800` `122:-161280` `124:-201600` `221:-60480`|
|6|16|`223:-100800` `224:-201600`|
|7|7|`5:1120` `11:-11200` `15:-33600` `18:-16800` `21:-6720` `22:-10080` `24:-20160` `26:-3360` `27:11200` `29:100800` `32:100800` `33:11200`|
|7|7|`34:33600` `36:302400` `37:302400` `38:16800` `39:151200` `40:1120` `41:-11200` `47:-100800` `66:-100800` `88:-11200` `111:-100800` `113:-134400`|
|7|7|`115:-201600` `118:-134400` `125:-907200` `130:-1209600` `142:-403200` `146:-806400` `154:-201600` `158:-907200` `163:-1209600` `185:-201600` `189:-403200` `199:-403200`|
|7|7|`210:-67200` `220:-100800` `222:-403200` `225:-67200` `226:-907200` `229:-604800` `232:-604800` `235:-201600` `236:-100800` `237:-100800` `240:-201600` `242:-100800`|
|7|7|`245:-403200` `246:-67200` `248:-100800` `249:-100800` `250:-907200` `252:-3628800` `256:-604800` `267:-5443200` `271:-3628800` `276:-604800` `278:-1209600` `282:-1209600`|
|7|7|`284:-907200` `285:-907200` `289:-907200` `296:-1209600` `303:-604800` `308:-100800` `311:-403200` `313:-201600` `316:-403200` `317:-403200` `319:-403200` `320:-907200`|
|7|7|`321:-2721600` `329:-302400` `331:-604800` `337:-1814400` `344:-604800` `349:-1209600` `362:-100800` `363:-302400` `368:-403200` `370:-100800` `373:-33600` `375:-11200`|
|7|7|`376:-134400` `377:-67200` `381:-1814400` `384:-403200` `385:-604800` `387:-1209600` `391:-403200` `395:-134400` `396:-134400` `400:-403200` `403:-134400` `407:-907200`|
|7|7|`410:-1209600` `414:-201600` `419:-302400` `423:-302400` `424:-302400` `427:-201600` `438:-1209600` `439:-3628800` `442:-3628800` `444:-2419200` `455:-604800` `457:-1209600`|
|7|7|`461:-2419200` `467:-1209600` `481:-3628800` `486:-604800` `493:-1814400` `497:-1209600` `498:-3628800` `505:-1209600` `511:-1209600` `522:-302400` `535:-907200` `543:-1814400`|
|7|7|`553:-302400` `561:-16800` `562:-50400` `563:-453600` `565:-100800` `572:-67200` `574:-100800` `581:-1814400` `586:-100800` `595:-907200` `599:-151200` `600:-453600`|
|7|7|`604:-201600` `623:-16800` `654:-453600` `666:-1814400` `676:-907200` `737:-10080` `748:-30240`|
|7|8|`21:134400` `24:403200` `26:134400` `28:403200` `31:1209600` `32:403200` `35:201600` `36:604800` `38:13440` `49:-134400` `69:-403200` `89:-134400`|
|7|8|`117:-268800` `127:-1209600` `133:-806400` `138:-806400` `145:-1612800` `149:-537600` `152:-806400` `159:-3628800` `163:-1209600` `164:-2419200` `172:-2419200` `186:-806400`|
|7|8|`188:-1612800` `193:-806400` `197:-1612800` `204:-806400` `209:-268800` `214:-537600` `228:-1209600` `231:-403200` `238:-403200` `243:-268800` `244:-806400` `251:-1209600`|
|7|8|`252:-3628800` `253:-2419200` `257:-403200` `261:-2419200` `263:-403200` `268:-7257600` `270:-14515200` `275:-2419200` `276:-2419200` `277:-2016000` `278:-2419200` `281:-2419200`|
|7|8|`287:-3628800` `290:-1209600` `291:-1209600` `292:-2419200` `294:-1209600` `295:-2419200` `304:-806400` `309:-537600` `310:-1209600` `312:-1612800` `314:-1612800` `318:-806400`|
|7|8|`320:-3628800` `322:-3628800` `323:-7257600` `330:-2419200` `331:-1209600` `332:-1209600` `336:-7257600` `342:-3628800` `343:-1209600` `344:-1209600` `345:-2419200` `347:-2419200`|
|7|8|`349:-2419200` `350:-1209600` `353:-1209600` `357:-1209600` `360:-1209600` `361:-403200` `364:-806400` `365:-1209600` `367:-2419200` `372:-268800` `374:-134400` `376:-134400`|
|7|8|`382:-806400` `383:-1612800` `385:-2419200` `386:-403200` `387:-2419200` `390:-806400` `397:-403200` `399:-806400` `408:-1209600` `409:-2419200` `413:-806400` `415:-403200`|
|7|8|`416:-403200` `421:-1209600` `426:-403200` `438:-4838400` `440:-1612800` `441:-3225600` `443:-9676800` `444:-4838400` `446:-1612800` `450:-1612800` `454:-2419200` `457:-2419200`|
|7|8|`459:-806400` `460:-4838400` `461:-4838400` `463:-2419200` `465:-2419200` `467:-2419200` `473:-806400` `476:-1612800` `480:-3628800` `483:-4838400` `485:-4838400` `488:-403200`|
|7|8|`490:-806400` `491:-806400` `492:-2419200` `495:-9676800` `496:-4838400` `502:-2419200` `504:-4838400` `507:-1209600` `508:-2419200` `510:-4838400` `512:-806400` `513:-806400`|
|7|8|`514:-806400` `516:-1612800` `518:-1612800` `519:-1209600` `534:-2419200` `537:-1209600` `540:-1209600` `541:-4838400` `542:-2419200` `546:-2419200` `547:-403200` `549:-806400`|
|7|8|`552:-2419200` `555:-806400` `556:-403200` `561:-67200` `564:-403200` `565:-201600` `566:-604800` `572:-134400` `573:-201600` `574:-201600` `580:-1814400` `583:-806400`|
|7|8|`585:-806400` `589:-1209600` `591:-1209600` `594:-1209600` `597:-1209600` `598:-604800` `603:-806400` `606:-604800` `608:-1209600` `609:-1209600` `614:-604800` `615:-604800`|
|7|8|`619:-604800` `653:-1209600` `657:-403200` `661:-604800` `665:-4838400` `669:-806400` `673:-604800` `674:-2419200` `675:-1209600` `681:-1209600` `684:-2419200` `688:-1209600`|
|7|8|`691:-604800` `692:-1209600` `697:-604800` `701:-1209600` `736:-26880` `742:-40320` `747:-80640` `754:-80640` `758:-120960` `765:-13440` `769:-40320` `770:-80640`|
|7|9|`15:201600` `18:201600` `22:604800` `24:604800` `28:302400` `29:302400` `34:20160` `72:-201600` `90:-201600` `144:-806400` `160:-1814400` `164:-1814400`|
|7|9|`165:-1209600` `182:-1209600` `187:-1612800` `192:-3628800` `212:-806400` `215:-403200` `231:-201600` `253:-1814400` `256:-1814400` `269:-4838400` `274:-8467200` `277:-4032000`|
|7|9|`281:-4838400` `290:-604800` `294:-1209600` `322:-7257600` `324:-1814400` `325:-9676800` `330:-604800` `332:-2419200` `340:-7257600` `347:-6048000` `350:-3024000` `359:-1209600`|
|7|9|`361:-201600` `365:-403200` `377:-201600` `386:-1209600` `390:-1209600` `392:-1209600` `394:-1209600` `411:-1209600` `412:-604800` `440:-2419200` `445:-4838400` `446:-4838400`|
|7|9|`451:-2419200` `456:-3628800` `459:-1209600` `465:-3628800` `469:-2419200` `470:-1209600` `475:-2419200` `482:-6048000` `483:-1209600` `484:-1814400` `499:-10886400` `500:-6048000`|
|7|9|`501:-2419200` `502:-1209600` `512:-2419200` `533:-1814400` `536:-5443200` `537:-1209600` `539:-3628800` `544:-1814400` `545:-604800` `546:-1209600` `567:-302400` `568:-604800`|
|7|9|`571:-201600` `576:-100800` `579:-302400` `582:-604800` `584:-302400` `588:-604800` `589:-604800` `601:-604800` `605:-907200` `610:-907200` `611:-907200` `612:-1209600`|
|7|9|`613:-604800` `614:-302400` `620:-604800` `652:-907200` `656:-604800` `668:-1209600` `672:-1814400` `679:-1814400` `680:-1814400` `681:-1209600` `690:-1814400` `693:-907200`|
|7|9|`694:-907200` `695:-907200` `696:-302400` `697:-604800` `735:-20160` `752:-120960` `760:-181440` `764:-40320` `768:-120960` `771:-60480` `772:-60480`|
|7|10|`11:44800` `15:134400` `21:67200` `27:4480` `91:-44800` `165:-940800` `187:-268800` `257:-672000` `269:-1209600` `277:-806400` `309:-134400` `324:-1209600`|
|7|10|`325:-1209600` `379:-44800` `386:-403200` `445:-806400` `446:-806400` `458:-1209600` `459:-806400` `482:-1209600` `488:-403200` `575:-67200` `576:-67200` `582:-201600`|
|7|10|`583:-201600` `655:-201600` `656:-604800` `657:-201600` `763:-13440` `764:-40320` `765:-13440`|
|7|11|`21:16800` `24:50400` `26:16800` `28:100800` `31:302400` `32:100800` `35:100800` `36:302400` `38:16800` `42:-16800` `44:-100800` `46:-50400`|
|7|11|`49:-67200` `50:-100800` `52:-50400` `54:-33600` `57:-302400` `59:-151200` `66:-100800` `69:-201600` `73:-100800` `76:-100800` `80:-50400` `83:-50400`|
|7|11|`89:-67200` `97:-16800` `104:-33600` `112:-201600` `114:-50400` `117:-67200` `119:-100800` `120:-33600` `121:-302400` `122:-403200` `123:-50400` `125:-907200`|
|7|11|`127:-604800` `128:-2419200` `132:-604800` `133:-100800` `134:-302400` `135:-302400` `137:-403200` `138:-201600` `139:-201600` `140:-100800` `142:-403200` `145:-403200`|
|7|11|`149:-134400` `150:-151200` `151:-201600` `152:-201600` `155:-201600` `156:-100800` `157:-16800` `159:-1814400` `163:-1209600` `164:-302400` `167:-1209600` `170:-151200`|
|7|11|`171:-302400` `172:-604800` `173:-1209600` `179:-604800` `180:-151200` `181:-151200` `186:-100800` `188:-403200` `189:-403200` `193:-201600` `195:-151200` `197:-403200`|
|7|11|`199:-403200` `204:-201600` `205:-100800` `206:-302400` `207:-50400` `209:-67200` `214:-134400` `218:-33600` `219:-16800` `220:-100800` `221:-604800` `223:-100800`|
|7|11|`224:-201600` `226:-907200` `228:-302400` `230:-604800` `234:-302400` `239:-201600` `241:-100800` `242:-100800` `251:-604800` `252:-3628800` `253:-302400` `255:-907200`|
|7|11|`258:-1209600` `260:-604800` `261:-604800` `262:-1209600` `264:-201600` `266:-403200` `268:-907200` `270:-3628800` `271:-3628800` `275:-604800` `276:-604800` `278:-1209600`|
|7|11|`280:-2419200` `283:-302400` `285:-907200` `287:-907200` `289:-907200` `291:-302400` `293:-1209600` `297:-1209600` `300:-604800` `301:-1209600` `302:-302400` `305:-100800`|
|7|11|`307:-201600` `311:-403200` `317:-403200` `320:-907200` `323:-1814400` `331:-604800` `334:-302400` `336:-1814400` `342:-907200` `344:-604800` `348:-1209600` `349:-1209600`|
|7|11|`352:-604800` `354:-302400` `358:-604800` `362:-100800` `368:-403200` `370:-100800` `376:-134400` `378:-100800` `382:-100800` `383:-403200` `384:-403200` `385:-604800`|
|7|11|`387:-1209600` `389:-604800` `396:-134400` `398:-403200` `401:-403200` `410:-1209600` `419:-302400` `424:-302400` `438:-1209600` `441:-806400` `444:-2419200` `449:-302400`|
|7|11|`454:-604800` `457:-1209600` `461:-2419200` `466:-1209600` `467:-1209600` `474:-201600` `477:-403200` `486:-604800` `497:-1209600` `505:-1209600` `511:-1209600` `522:-302400`|
|7|11|`561:-16800` `565:-100800` `572:-67200` `574:-100800` `586:-100800` `599:-151200` `604:-201600` `623:-16800`|
|7|12|`15:67200` `18:67200` `22:403200` `24:403200` `28:403200` `29:403200` `34:67200` `53:-134400` `63:-403200` `64:-201600` `69:-403200` `72:-268800`|
|7|12|`79:-604800` `84:-67200` `90:-268800` `102:-134400` `107:-67200` `117:-67200` `127:-1209600` `130:-1209600` `136:-3628800` `138:-201600` `144:-537600` `145:-2016000`|
|7|12|`149:-134400` `152:-403200` `160:-2419200` `164:-2419200` `165:-403200` `172:-3024000` `177:-3024000` `182:-806400` `183:-1814400` `187:-672000` `188:-201600` `192:-2419200`|
|7|12|`193:-1008000` `197:-2419200` `204:-1209600` `209:-67200` `212:-537600` `214:-134400` `215:-268800` `222:-403200` `228:-1209600` `232:-1209600` `238:-806400` `244:-403200`|
|7|12|`253:-2419200` `256:-2419200` `261:-2419200` `265:-1612800` `269:-1209600` `274:-4838400` `275:-2419200` `277:-1209600` `281:-2419200` `282:-2419200` `284:-3628800` `287:-3628800`|
|7|12|`291:-1209600` `292:-1209600` `295:-7257600` `296:-2419200` `303:-1209600` `304:-1612800` `306:-1612800` `312:-806400` `314:-1612800` `318:-806400` `322:-3628800` `325:-2419200`|
|7|12|`332:-1209600` `340:-3628800` `343:-1209600` `345:-2419200` `347:-2419200` `350:-1209600` `353:-3628800` `357:-1209600` `359:-403200` `360:-403200` `364:-806400` `365:-403200`|
|7|12|`367:-806400` `377:-268800` `386:-403200` `390:-806400` `391:-806400` `392:-403200` `394:-806400` `395:-537600` `397:-403200` `399:-806400` `400:-806400` `403:-268800`|
|7|12|`413:-806400` `417:-806400` `418:-1209600` `421:-1209600` `426:-403200` `429:-806400` `431:-403200` `440:-1612800` `446:-1612800` `451:-1612800` `456:-3628800` `459:-806400`|
|7|12|`463:-2419200` `465:-2419200` `469:-1612800` `470:-806400` `473:-806400` `476:-1612800` `488:-403200` `490:-403200` `503:-4838400` `507:-1209600` `508:-2419200` `512:-806400`|
|7|12|`515:-1612800` `516:-806400` `519:-1209600` `520:-1209600` `523:-403200` `525:-806400` `568:-201600` `571:-201600` `576:-67200` `579:-201600` `593:-201600` `602:-806400`|
|7|12|`617:-604800` `620:-201600` `622:-134400` `628:-201600` `630:-201600`|
|7|13|`11:33600` `15:201600` `21:201600` `27:33600` `72:-403200` `90:-100800` `91:-134400` `133:-806400` `144:-806400` `160:-302400` `164:-604800` `165:-2016000`|
|7|13|`182:-1209600` `186:-100800` `187:-806400` `192:-1209600` `225:-201600` `231:-806400` `257:-1209600` `269:-1209600` `274:-1209600` `277:-1612800` `281:-1209600` `290:-2419200`|
|7|13|`294:-1814400` `322:-1814400` `324:-1814400` `325:-2419200` `330:-604800` `332:-604800` `379:-134400` `386:-403200` `392:-403200` `402:-403200` `404:-403200` `408:-403200`|
|7|13|`411:-604800` `412:-604800` `458:-1209600` `459:-1209600` `472:-604800` `482:-1209600` `483:-1209600` `484:-1814400` `487:-604800` `575:-100800` `576:-100800` `605:-302400`|
|7|13|`606:-302400` `621:-100800` `624:-302400` `625:-100800`|
|7|14|`15:3360` `18:3360` `22:50400` `24:50400` `28:100800` `29:100800` `34:33600` `42:-3360` `44:-100800` `45:-221760` `46:-10080` `47:-100800`|
|7|14|`50:-120960` `52:-20160` `53:-33600` `54:-6720` `59:-181440` `62:-181440` `63:-201600` `64:-50400` `65:-110880` `69:-201600` `72:-33600` `73:-10080`|
|7|14|`76:-141120` `79:-151200` `80:-60480` `83:-70560` `84:-16800` `90:-33600` `97:-3360` `102:-33600` `104:-6720` `107:-16800` `111:-100800` `112:-201600`|
|7|14|`114:-100800` `115:-201600` `119:-50400` `122:-403200` `124:-201600` `127:-604800` `130:-1209600` `132:-604800` `134:-151200` `135:-907200` `136:-806400` `137:-403200`|
|7|14|`139:-201600` `140:-201600` `141:-201600` `145:-403200` `146:-806400` `151:-201600` `154:-201600` `155:-100800` `156:-100800` `160:-302400` `164:-302400` `170:-151200`|
|7|14|`171:-302400` `172:-604800` `177:-604800` `180:-453600` `181:-151200` `183:-403200` `184:-50400` `193:-201600` `195:-50400` `197:-403200` `204:-201600` `205:-100800`|
|7|14|`206:-100800` `207:-50400` `222:-403200` `223:-100800` `224:-201600` `228:-302400` `229:-604800` `232:-604800` `234:-302400` `237:-100800` `239:-201600` `241:-100800`|
|7|14|`245:-403200` `249:-100800` `253:-302400` `256:-604800` `260:-604800` `261:-604800` `264:-201600` `266:-403200` `275:-604800` `282:-1209600` `283:-302400` `284:-907200`|
|7|14|`287:-907200` `291:-302400` `296:-1209600` `300:-604800` `302:-302400` `303:-604800` `305:-100800` `307:-201600` `313:-201600` `316:-403200` `319:-403200` `377:-67200`|
|7|14|`391:-403200` `395:-134400` `400:-403200` `403:-134400` `414:-201600` `423:-302400` `427:-201600`|
|7|15|`11:4480` `15:67200` `21:134400` `27:44800` `49:-161280` `53:-134400` `63:-40320` `64:-201600` `69:-80640` `72:-403200` `79:-201600` `84:-67200`|
|7|15|`89:-13440` `90:-67200` `91:-44800` `113:-134400` `117:-336000` `133:-806400` `136:-403200` `138:-1008000` `144:-537600` `145:-403200` `149:-672000` `152:-806400`|
|7|15|`165:-403200` `172:-604800` `182:-806400` `183:-201600` `187:-268800` `188:-201600` `192:-806400` `193:-201600` `225:-268800` `231:-403200` `235:-403200` `236:-403200`|
|7|15|`238:-403200` `240:-403200` `243:-268800` `244:-403200` `257:-403200` `263:-403200` `277:-806400` `290:-1209600` `292:-1209600` `294:-1209600` `304:-806400` `309:-134400`|
|7|15|`310:-403200` `312:-806400` `379:-44800` `392:-268800` `402:-403200` `404:-268800` `415:-134400` `425:-806400` `428:-403200` `432:-134400`|
|7|16|`15:3360` `21:16800` `27:11200` `41:-11200` `42:-20160` `45:-20160` `46:-60480` `49:-67200` `50:-20160` `52:-50400` `53:-33600` `54:-40320`|
|7|16|`59:-30240` `64:-50400` `65:-10080` `72:-33600` `73:-10080` `79:-50400` `80:-10080` `84:-16800` `113:-134400` `114:-50400` `117:-67200` `118:-134400`|
|7|16|`119:-50400` `120:-33600` `123:-50400` `133:-100800` `134:-151200` `138:-201600` `140:-100800` `149:-134400` `150:-50400` `152:-201600` `155:-100800` `157:-16800`|
|7|16|`225:-67200` `235:-201600` `236:-100800` `240:-201600` `246:-67200` `248:-100800`|
|8|8|`4:3360` `12:-33600` `16:554400` `17:-67200` `19:772800` `20:100800` `21:-10080` `23:1955520` `24:-20160` `25:2661120` `26:332640` `30:1512000`|
|8|8|`31:2016000` `35:362880` `77:-604800` `78:-201600` `92:-806400` `93:-268800` `96:-100800` `97:-33600` `147:-2419200` `153:-806400` `161:-5443200` `162:-1814400`|
|8|8|`166:-7257600` `167:-2419200` `170:-907200` `174:-7257600` `175:-2419200` `190:-14515200` `191:-3225600` `194:-1209600` `195:-403200` `196:-4838400` `198:-1612800` `203:-1612800`|
|8|8|`208:-806400` `209:-268800` `211:-3225600` `213:-1075200` `216:-1612800` `217:-940800` `218:-134400` `233:-604800` `247:-201600` `254:-5443200` `255:-1814400` `258:-7257600`|
|8|8|`272:-29030400` `273:-7257600` `279:-33868800` `280:-9676800` `283:-3628800` `286:-10886400` `288:-3628800` `290:-3628800` `294:-1209600` `298:-2419200` `299:-2419200` `304:-4838400`|
|8|8|`305:-604800` `312:-1612800` `313:-201600` `315:-2419200` `322:-7257600` `326:-5443200` `327:-21772800` `328:-2721600` `330:-1209600` `333:-12096000` `334:-3628800` `335:-21772800`|
|8|8|`336:-7257600` `339:-36288000` `341:-7257600` `343:-3628800` `346:-24192000` `347:-2419200` `348:-4838400` `351:-9676800` `352:-4233600` `353:-4838400` `355:-3024000` `356:-7257600`|
|8|8|`357:-4838400` `358:-1209600` `360:-806400` `361:-604800` `364:-1612800` `366:-4838400` `367:-3225600` `368:-403200` `369:-1612800` `370:-100800` `371:-268800` `372:-806400`|
|8|8|`373:-33600` `378:-604800` `388:-7257600` `389:-2419200` `393:-9676800` `411:-3628800` `412:-1209600` `413:-4838400` `414:-604800` `417:-1612800` `420:-1814400` `422:-604800`|
|8|8|`447:-7257600` `448:-19353600` `449:-3628800` `452:-12902400` `453:-7257600` `454:-2419200` `460:-14515200` `462:-14515200` `464:-19353600` `466:-4838400` `468:-9676800` `471:-4838400`|
|8|8|`484:-7257600` `485:-9676800` `486:-1209600` `487:-2419200` `491:-3225600` `494:-50803200` `495:-9676800` `496:-2419200` `499:-14515200` `501:-4838400` `506:-18144000` `507:-4838400`|
|8|8|`508:-9676800` `509:-19353600` `510:-9676800` `511:-2419200` `515:-3225600` `517:-3225600` `518:-3225600` `521:-1209600` `533:-5443200` `538:-29030400` `539:-7257600` `540:-3628800`|
|8|8|`541:-9676800` `542:-2419200` `543:-604800` `545:-3628800` `547:-1612800` `548:-806400` `550:-7257600` `551:-4838400` `552:-9676800` `553:-302400` `554:-806400` `555:-1612800`|
|8|8|`569:-907200` `570:-2419200` `573:-604800` `578:-1209600` `590:-3628800` `591:-4838400` `592:-1209600` `596:-4536000` `597:-1209600` `598:-302400` `607:-7257600` `611:-3628800`|
|8|8|`613:-1209600` `616:-3628800` `617:-2419200` `618:-3628800` `619:-2419200` `658:-907200` `659:-1209600` `662:-302400` `664:-3628800` `671:-10886400` `673:-1209600` `674:-2419200`|
|8|8|`682:-7257600` `683:-7257600` `684:-4838400` `685:-2419200` `687:-2419200` `689:-7257600` `693:-1814400` `696:-604800` `698:-3628800` `699:-3628800` `700:-2419200` `701:-2419200`|
|8|8|`743:-181440` `746:-60480` `756:-362880` `761:-725760` `767:-483840` `771:-120960` `773:-241920`|
|8|9|`12:403200` `14:201600` `17:1612800` `19:806400` `23:1814400` `24:907200` `28:645120` `94:-403200` `95:-604800` `98:-201600` `168:-6048000` `169:-5443200`|
|8|9|`171:-1814400` `174:-2419200` `175:-1209600` `182:-2419200` `183:-1209600` `190:-4032000` `191:-1209600` `192:-3225600` `193:-1209600` `200:-4032000` `203:-403200` `204:-403200`|
|8|9|`211:-3225600` `212:-806400` `213:-403200` `214:-403200` `259:-4838400` `260:-5443200` `263:-1209600` `264:-604800` `272:-10886400` `273:-3628800` `274:-10886400` `275:-3628800`|
|8|9|`279:-9676800` `280:-3628800` `281:-8467200` `282:-3628800` `290:-3628800` `291:-1814400` `292:-3628800` `293:-604800` `294:-1209600` `299:-4838400` `302:-604800` `303:-604800`|
|8|9|`304:-2419200` `306:-2419200` `307:-604800` `310:-403200` `312:-806400` `315:-806400` `318:-806400` `319:-201600` `322:-7257600` `323:-3628800` `330:-3628800` `331:-1814400`|
|8|9|`338:-10886400` `339:-29030400` `340:-7257600` `341:-1814400` `342:-1814400` `345:-9676800` `346:-12096000` `347:-6048000` `348:-2419200` `349:-2419200` `353:-4838400` `354:-604800`|
|8|9|`356:-1209600` `360:-403200` `364:-806400` `366:-806400` `367:-806400` `368:-201600` `380:-403200` `388:-2419200` `389:-1209600` `390:-2419200` `391:-1209600` `393:-2419200`|
|8|9|`394:-2419200` `397:-1209600` `398:-604800` `409:-1209600` `410:-604800` `440:-4838400` `441:-2419200` `443:-9676800` `444:-4838400` `450:-4838400` `456:-2419200` `457:-1209600`|
|8|9|`462:-18144000` `463:-16934400` `464:-9676800` `465:-4838400` `466:-2419200` `467:-2419200` `468:-4838400` `469:-2419200` `475:-7257600` `476:-9676800` `477:-1209600` `480:-2419200`|
|8|9|`481:-1209600` `494:-14515200` `495:-9676800` `496:-2419200` `497:-1209600` `498:-1209600` `499:-14515200` `503:-7257600` `504:-7257600` `505:-1814400` `506:-7257600` `507:-1209600`|
|8|9|`508:-2419200` `509:-2419200` `510:-2419200` `511:-604800` `533:-3628800` `534:-4838400` `535:-604800` `538:-7257600` `539:-7257600` `540:-1209600` `541:-4838400` `542:-2419200`|
|8|9|`543:-604800` `564:-604800` `565:-302400` `566:-1209600` `577:-604800` `578:-403200` `579:-201600` `580:-1209600` `581:-604800` `587:-1209600` `596:-1814400` `597:-2419200`|
|8|9|`598:-604800` `599:-302400` `600:-302400` `601:-1814400` `602:-1209600` `603:-1209600` `604:-302400` `607:-1814400` `608:-1209600` `609:-604800` `610:-1814400` `615:-1209600`|
|8|9|`664:-3628800` `665:-4838400` `666:-604800` `670:-1814400` `671:-3628800` `672:-5443200` `673:-604800` `674:-2419200` `675:-1209600` `676:-302400` `677:-3628800` `678:-2419200`|
|8|9|`689:-3628800` `690:-1814400` `691:-604800` `692:-1209600` `749:-362880` `750:-241920` `766:-120960` `767:-241920` `768:-120960` `769:-40320` `770:-80640`|
|8|10|`8:33600` `12:201600` `16:352800` `21:154560` `72:-33600` `90:-100800` `133:-201600` `144:-268800` `149:-134400` `160:-302400` `164:-1512000` `192:-806400`|
|8|10|`193:-403200` `225:-100800` `231:-604800` `235:-201600` `243:-268800` `246:-33600` `253:-907200` `256:-1209600` `268:-1814400` `274:-1209600` `275:-604800` `281:-1209600`|
|8|10|`282:-604800` `294:-604800` `312:-806400` `313:-100800` `322:-3628800` `323:-1814400` `377:-134400` `382:-201600` `385:-907200` `390:-806400` `391:-403200` `411:-604800`|
|8|10|`413:-806400` `414:-100800` `440:-806400` `441:-403200` `443:-2419200` `444:-1209600` `456:-2419200` `457:-1209600` `484:-1814400` `485:-2419200` `486:-302400` `564:-201600`|
|8|10|`565:-100800` `571:-67200` `572:-33600` `580:-604800` `581:-302400` `584:-302400` `585:-403200` `586:-50400` `652:-907200` `653:-1209600` `654:-151200` `735:-20160`|
|8|10|`736:-26880` `737:-3360`|
|8|11|`16:100800` `19:134400` `20:16800` `23:604800` `25:806400` `26:100800` `30:604800` `31:806400` `35:100800` `51:-201600` `58:-604800` `60:-302400`|
|8|11|`67:-806400` `73:-100800` `74:-604800` `75:-403200` `77:-403200` `81:-201600` `85:-50400` `92:-537600` `96:-134400` `98:-134400` `103:-134400` `105:-134400`|
|8|11|`109:-33600` `116:-100800` `126:-1814400` `128:-2419200` `129:-3024000` `131:-2419200` `133:-302400` `143:-3628800` `147:-806400` `148:-806400` `149:-403200` `150:-403200`|
|8|11|`157:-50400` `161:-3628800` `164:-604800` `166:-9676800` `169:-3628800` `170:-604800` `174:-1814400` `175:-1814400` `178:-2419200` `186:-100800` `188:-806400` `190:-4032000`|
|8|11|`191:-4032000` `193:-403200` `194:-705600` `195:-201600` `196:-1612800` `198:-1612800` `201:-403200` `202:-2016000` `203:-403200` `204:-403200` `206:-403200` `207:-100800`|
|8|11|`208:-268800` `211:-1075200` `213:-134400` `214:-134400` `216:-537600` `217:-403200` `219:-16800` `221:-604800` `227:-1814400` `230:-2419200` `231:-604800` `233:-604800`|
|8|11|`243:-806400` `246:-100800` `254:-7257600` `258:-4838400` `259:-2419200` `270:-7257600` `272:-5443200` `273:-7257600` `279:-16934400` `280:-4838400` `281:-2419200` `283:-907200`|
|8|11|`286:-3628800` `288:-3628800` `293:-2419200` `294:-1209600` `297:-2419200` `298:-2419200` `299:-3628800` `301:-2419200` `310:-806400` `312:-3225600` `313:-201600` `315:-1612800`|
|8|11|`318:-1612800` `319:-403200` `327:-1814400` `330:-1209600` `333:-4838400` `334:-302400` `335:-7257600` `339:-7257600` `343:-1209600` `346:-12096000` `347:-1209600` `351:-4838400`|
|8|11|`352:-604800` `355:-604800` `356:-1209600` `360:-403200` `364:-806400` `366:-1612800` `367:-2419200` `368:-201600` `369:-806400` `370:-100800` `378:-403200` `383:-806400`|
|8|11|`388:-3628800` `389:-1209600` `393:-2419200` `398:-806400` `401:-806400` `409:-2419200` `412:-2419200` `417:-1612800` `422:-1209600` `428:-806400` `430:-806400` `431:-806400`|
|8|11|`443:-2419200` `448:-4838400` `452:-3225600` `453:-2419200` `460:-4838400` `464:-9676800` `468:-3225600` `471:-1612800` `485:-2419200` `487:-1209600` `491:-806400` `495:-4838400`|
|8|11|`501:-2419200` `504:-2419200` `509:-4838400` `510:-2419200` `515:-1612800` `517:-1612800` `518:-1612800` `521:-604800` `524:-806400` `525:-806400` `526:-806400` `564:-100800`|
|8|11|`570:-302400` `573:-201600` `578:-403200` `585:-403200` `592:-302400` `597:-604800` `603:-403200` `613:-604800` `618:-604800` `625:-100800` `629:-403200` `631:-100800`|
|8|12|`12:201600` `14:100800` `17:1276800` `19:638400` `23:1612800` `24:806400` `28:604800` `70:-1814400` `71:-403200` `75:-604800` `76:-201600` `77:-604800`|
|8|12|`81:-504000` `82:-201600` `86:-100800` `92:-403200` `94:-806400` `95:-873600` `98:-134400` `100:-403200` `103:-369600` `105:-100800` `106:-67200` `108:-33600`|
|8|12|`109:-33600` `131:-4233600` `132:-1209600` `135:-1814400` `138:-604800` `140:-302400` `147:-4838400` `148:-4032000` `151:-806400` `152:-201600` `153:-201600` `155:-100800`|
|8|12|`161:-1814400` `166:-1209600` `168:-9676800` `169:-6652800` `171:-1209600` `172:-1209600` `174:-2419200` `175:-1814400` `176:-7257600` `178:-4838400` `180:-302400` `181:-302400`|
|8|12|`182:-1209600` `183:-604800` `188:-201600` `190:-4838400` `191:-2016000` `192:-1209600` `193:-201600` `194:-806400` `195:-403200` `196:-9676800` `197:-2016000` `198:-806400`|
|8|12|`200:-4838400` `201:-3225600` `202:-1612800` `203:-1814400` `204:-403200` `205:-504000` `206:-604800` `209:-67200` `211:-1075200` `212:-268800` `213:-403200` `214:-134400`|
|8|12|`217:-134400` `218:-33600` `224:-1209600` `233:-3628800` `234:-3024000` `238:-1209600` `239:-604800` `244:-403200` `245:-201600` `259:-7257600` `260:-6048000` `272:-7257600`|
|8|12|`273:-3628800` `274:-3628800` `275:-1814400` `279:-9676800` `280:-4838400` `281:-1209600` `282:-604800` `283:-2419200` `286:-21772800` `287:-3628800` `292:-2419200` `293:-1209600`|
|8|12|`295:-12096000` `298:-14515200` `299:-7257600` `300:-5443200` `301:-3024000` `302:-1209600` `310:-806400` `311:-403200` `314:-3225600` `315:-806400` `316:-403200` `317:-403200`|
|8|12|`322:-3628800` `323:-1814400` `330:-1209600` `331:-604800` `335:-7257600` `338:-10886400` `339:-29030400` `340:-3628800` `341:-1814400` `342:-1814400` `343:-1209600` `345:-4838400`|
|8|12|`346:-12096000` `347:-3628800` `348:-1209600` `349:-1209600` `351:-2419200` `352:-604800` `355:-1209600` `356:-1209600` `360:-2016000` `361:-403200` `362:-201600` `363:-201600`|
|8|12|`364:-806400` `366:-806400` `367:-806400` `368:-201600` `380:-806400` `388:-2419200` `389:-1209600` `390:-1209600` `391:-604800` `393:-3225600` `394:-403200` `399:-2419200`|
|8|12|`405:-1612800` `406:-806400` `409:-2419200` `410:-1209600` `416:-806400` `420:-3628800` `421:-6048000` `422:-1209600` `423:-604800` `424:-604800` `429:-806400` `430:-403200`|
|8|12|`443:-4838400` `444:-2419200` `450:-1612800` `462:-10886400` `463:-7257600` `464:-4838400` `465:-3628800` `466:-1209600` `467:-1209600` `468:-1612800` `469:-806400` `480:-1209600`|
|8|12|`481:-604800` `489:-403200` `494:-14515200` `495:-9676800` `496:-2419200` `497:-1209600` `498:-1209600` `499:-3628800` `506:-3628800` `507:-3628800` `508:-4838400` `509:-2419200`|
|8|12|`510:-2419200` `511:-604800` `513:-1612800` `514:-806400` `519:-1209600` `523:-403200` `524:-806400` `566:-604800` `577:-604800` `578:-403200` `579:-201600` `587:-604800`|
|8|12|`601:-604800` `607:-1814400` `608:-1209600` `609:-604800` `626:-1209600` `627:-604800` `628:-201600` `629:-403200`|
|8|13|`8:50400` `12:352800` `16:604800` `21:352800` `60:-151200` `68:-504000` `72:-201600` `74:-554400` `90:-201600` `95:-201600` `98:-100800` `123:-151200`|
|8|13|`129:-3024000` `133:-604800` `134:-1360800` `143:-1612800` `144:-806400` `149:-403200` `150:-50400` `160:-1814400` `164:-3024000` `168:-2419200` `169:-1814400` `171:-302400`|
|8|13|`175:-302400` `186:-604800` `190:-2419200` `191:-1612800` `192:-806400` `193:-403200` `195:-50400` `223:-604800` `227:-907200` `230:-1209600` `253:-3628800` `256:-2419200`|
|8|13|`259:-2419200` `260:-1209600` `268:-3628800` `272:-5443200` `273:-1814400` `274:-3628800` `275:-1814400` `276:-2419200` `279:-4838400` `280:-1209600` `281:-1209600` `282:-604800`|
|8|13|`283:-302400` `292:-3628800` `293:-1814400` `322:-3628800` `323:-1814400` `327:-1814400` `330:-1209600` `331:-604800` `333:-2419200` `334:-302400` `377:-201600` `382:-604800`|
|8|13|`385:-1814400` `388:-1209600` `389:-604800` `397:-806400` `398:-403200` `407:-604800` `409:-1209600` `410:-604800` `440:-2419200` `441:-1209600` `443:-2419200` `444:-1209600`|
|8|13|`447:-1814400` `448:-2419200` `449:-302400` `456:-2419200` `457:-1209600` `480:-2419200` `481:-1209600` `484:-1814400` `485:-2419200` `486:-302400` `571:-201600` `572:-100800`|
|8|13|`580:-604800` `581:-302400` `584:-302400` `585:-403200` `586:-50400`|
|8|14|`12:13440` `14:6720` `17:201600` `19:100800` `23:403200` `24:201600` `28:134400` `48:-443520` `49:-40320` `50:-201600` `51:-403200` `54:-20160`|
|8|14|`55:-201600` `58:-120960` `61:-604800` `67:-80640` `69:-80640` `70:-1008000` `71:-161280` `75:-1209600` `77:-201600` `78:-40320` `80:-20160` `81:-302400`|
|8|14|`82:-201600` `83:-20160` `86:-201600` `87:-100800` `89:-13440` `92:-134400` `93:-53760` `94:-134400` `95:-134400` `100:-134400` `101:-26880` `103:-33600`|
|8|14|`104:-6720` `105:-33600` `112:-403200` `116:-806400` `117:-201600` `119:-403200` `120:-100800` `131:-6652800` `135:-1209600` `145:-2419200` `147:-1612800` `148:-4838400`|
|8|14|`151:-806400` `152:-403200` `153:-1209600` `155:-806400` `156:-604800` `166:-1209600` `168:-1209600` `169:-1209600` `176:-2419200` `178:-1209600` `188:-201600` `190:-806400`|
|8|14|`191:-806400` `196:-3225600` `197:-403200` `198:-403200` `200:-806400` `201:-806400` `202:-806400` `203:-201600` `205:-100800` `206:-201600` `224:-806400` `228:-1209600`|
|8|14|`233:-2419200` `234:-1814400` `239:-806400` `241:-403200` `244:-2419200` `245:-403200` `247:-403200` `248:-201600` `249:-201600` `259:-2419200` `260:-604800` `279:-4838400`|
|8|14|`280:-1209600` `281:-1209600` `286:-3628800` `295:-2419200` `298:-4838400` `299:-2419200` `300:-604800` `301:-604800` `310:-806400` `311:-201600` `312:-806400` `314:-3225600`|
|8|14|`315:-1612800` `316:-403200` `317:-403200` `318:-806400` `319:-201600` `380:-134400` `390:-403200` `393:-806400` `399:-806400` `405:-537600` `406:-268800` `413:-806400`|
|8|14|`416:-403200` `421:-1209600` `426:-403200` `429:-806400` `430:-403200` `432:-134400` `433:-268800`|
|8|15|`8:10080` `12:154560` `16:352800` `21:201600` `46:-30240` `51:-604800` `52:-110880` `55:-100800` `58:-181440` `59:-60480` `60:-151200` `67:-80640`|
|8|15|`68:-806400` `70:-604800` `71:-40320` `72:-100800` `73:-10080` `74:-504000` `75:-302400` `76:-20160` `78:-20160` `81:-201600` `85:-50400` `86:-100800`|
|8|15|`90:-33600` `93:-26880` `95:-67200` `97:-3360` `98:-33600` `114:-453600` `116:-604800` `119:-252000` `129:-3024000` `131:-1209600` `134:-907200` `135:-604800`|
|8|15|`143:-1612800` `147:-2419200` `148:-1612800` `150:-504000` `151:-201600` `152:-806400` `155:-403200` `160:-907200` `161:-907200` `164:-907200` `166:-1209600` `168:-2419200`|
|8|15|`169:-1814400` `170:-151200` `171:-302400` `175:-302400` `186:-201600` `190:-806400` `191:-806400` `194:-201600` `195:-151200` `223:-604800` `227:-907200` `230:-1814400`|
|8|15|`233:-1209600` `234:-604800` `238:-1209600` `239:-604800` `242:-201600` `244:-403200` `245:-201600` `253:-907200` `256:-302400` `268:-1814400` `272:-5443200` `273:-1814400`|
|8|15|`276:-907200` `279:-4838400` `280:-1209600` `283:-302400` `292:-1209600` `293:-604800` `294:-604800` `308:-100800` `310:-806400` `311:-403200` `312:-806400` `313:-100800`|
|8|15|`377:-100800` `385:-604800` `397:-806400` `398:-403200` `407:-302400` `409:-1209600` `410:-604800` `411:-604800` `413:-806400` `414:-100800`|
|8|16|`12:10080` `16:50400` `21:33600` `42:-43680` `43:-50400` `48:-80640` `50:-40320` `51:-504000` `52:-20160` `55:-201600` `58:-60480` `67:-80640`|
|8|16|`68:-100800` `70:-201600` `71:-40320` `73:-10080` `74:-50400` `75:-100800` `76:-20160` `114:-201600` `116:-403200` `117:-268800` `119:-403200` `120:-134400`|
|8|16|`129:-604800` `143:-806400` `147:-806400` `148:-806400` `150:-100800` `151:-201600` `223:-100800` `230:-302400` `231:-201600` `238:-403200` `239:-201600` `242:-100800`|
|8|16|`243:-268800` `244:-403200` `245:-201600` `246:-33600`|
|9|9|`3:2016` `9:100800` `10:20160` `12:-20160` `13:705600` `14:141120` `15:-30240` `17:1411200` `18:282240` `22:887040` `23:-20160` `24:-10080`|
|9|9|`77:-100800` `78:-20160` `92:-403200` `93:-80640` `96:-100800` `97:-20160` `138:-1209600` `139:-241920` `147:-806400` `152:-403200` `153:-201600` `154:-40320`|
|9|9|`161:-907200` `162:-181440` `166:-3628800` `167:-725760` `170:-907200` `172:-4838400` `173:-967680` `188:-1209600` `189:-241920` `196:-3225600` `197:-1612800` `198:-806400`|
|9|9|`199:-161280` `208:-806400` `209:-604800` `210:-40320` `233:-604800` `236:-907200` `237:-181440` `238:-2419200` `240:-604800` `241:-604800` `244:-806400` `247:-201600`|
|9|9|`248:-100800` `249:-100800` `254:-907200` `255:-181440` `258:-3628800` `261:-3628800` `262:-725760` `270:-14515200` `271:-2903040` `278:-3628800` `284:-4536000` `285:-907200`|
|9|9|`286:-3628800` `287:-1814400` `288:-907200` `289:-181440` `295:-9676800` `296:-2419200` `297:-2419200` `298:-2419200` `314:-3225600` `315:-806400` `316:-403200` `317:-403200`|
|9|9|`320:-4536000` `321:-907200` `329:-907200` `335:-14515200` `336:-10886400` `337:-725760` `343:-6048000` `344:-3024000` `355:-604800` `360:-806400` `361:-201600` `362:-100800`|
|9|9|`363:-100800` `378:-100800` `383:-1209600` `384:-241920` `387:-4838400` `395:-1209600` `396:-241920` `399:-4838400` `400:-1209600` `401:-1209600` `418:-1209600` `419:-302400`|
|9|9|`420:-1814400` `421:-2419200` `422:-604800` `423:-302400` `424:-302400` `438:-6048000` `439:-1209600` `442:-6048000` `453:-2419200` `454:-1814400` `455:-120960` `460:-29030400`|
|9|9|`461:-14515200` `492:-9676800` `493:-2419200` `494:-7257600` `495:-9676800` `496:-2419200` `497:-1209600` `498:-1209600` `532:-1209600` `533:-1814400` `534:-2419200` `535:-302400`|
|9|9|`561:-151200` `562:-30240` `563:-756000` `573:-1411200` `574:-705600` `594:-6048000` `595:-1512000` `596:-907200` `597:-1209600` `598:-302400` `599:-151200` `600:-151200`|
|9|9|`663:-4838400` `664:-3628800` `665:-4838400` `666:-604800` `745:-120960` `746:-60480` `747:-80640` `748:-10080`|
|9|10|`9:40320` `12:161280` `15:181440` `69:-40320` `89:-40320` `130:-241920` `145:-403200` `146:-80640` `159:-362880` `163:-604800` `188:-403200` `189:-80640`|
|9|10|`232:-604800` `244:-403200` `245:-201600` `252:-1088640` `267:-725760` `270:-1814400` `271:-362880` `278:-604800` `310:-403200` `311:-201600` `320:-1814400` `321:-362880`|
|9|10|`381:-241920` `387:-1209600` `409:-1209600` `410:-604800` `438:-1209600` `439:-241920` `442:-1209600` `480:-1209600` `481:-604800` `563:-302400` `580:-604800` `581:-302400`|
|9|11|`12:100800` `14:50400` `17:604800` `19:302400` `23:604800` `24:302400` `28:100800` `68:-907200` `72:-302400` `74:-302400` `84:-151200` `90:-100800`|
|9|11|`94:-403200` `95:-604800` `106:-50400` `107:-50400` `129:-2116800` `134:-907200` `144:-2419200` `149:-1209600` `160:-907200` `164:-302400` `168:-8467200` `169:-2419200`|
|9|11|`174:-604800` `175:-604800` `182:-1209600` `183:-1814400` `184:-151200` `190:-3225600` `192:-2822400` `193:-806400` `194:-100800` `200:-3225600` `203:-1008000` `204:-1008000`|
|9|11|`207:-50400` `223:-604800` `231:-1814400` `235:-907200` `259:-3628800` `260:-1814400` `272:-9072000` `273:-907200` `274:-7257600` `275:-1814400` `281:-3628800` `282:-1814400`|
|9|11|`292:-3628800` `294:-1814400` `299:-3628800` `302:-907200` `303:-907200` `304:-2419200` `306:-2419200` `307:-604800` `322:-3628800` `323:-907200` `327:-1814400` `338:-5443200`|
|9|11|`339:-10886400` `340:-3628800` `341:-907200` `342:-907200` `353:-3628800` `354:-302400` `355:-604800` `356:-1209600` `357:-1209600` `358:-302400` `380:-403200` `388:-604800`|
|9|11|`389:-302400` `390:-1209600` `391:-604800` `394:-1209600` `397:-1209600` `404:-403200` `406:-806400` `411:-1814400` `425:-1209600` `426:-1209600` `427:-302400` `443:-2419200`|
|9|11|`444:-1209600` `447:-1814400` `462:-3628800` `463:-2419200` `464:-2419200` `465:-1209600` `466:-604800` `467:-604800` `484:-3628800` `499:-3628800` `503:-2419200` `504:-2419200`|
|9|11|`505:-604800` `506:-3628800` `519:-2419200` `520:-1209600` `521:-604800` `522:-302400` `566:-302400` `577:-302400` `578:-201600` `579:-100800` `584:-302400` `610:-907200`|
|9|11|`615:-604800` `624:-302400` `627:-604800` `630:-201600` `631:-100800`|
|9|12|`9:100800` `10:20160` `13:705600` `14:141120` `17:1209600` `18:241920` `22:705600` `64:-302400` `65:-60480` `70:-1008000` `71:-120960` `75:-604800`|
|9|12|`76:-120960` `77:-403200` `78:-80640` `79:-403200` `80:-60480` `81:-100800` `82:-100800` `83:-20160` `92:-403200` `93:-80640` `100:-403200` `102:-100800`|
|9|12|`103:-201600` `104:-20160` `131:-2419200` `132:-362880` `135:-1814400` `136:-3628800` `137:-483840` `138:-1814400` `139:-362880` `140:-302400` `141:-302400` `145:-2822400`|
|9|12|`146:-483840` `147:-1612800` `148:-403200` `152:-806400` `153:-604800` `154:-80640` `155:-100800` `156:-100800` `161:-3628800` `162:-725760` `166:-3628800` `167:-725760`|
|9|12|`172:-3024000` `173:-604800` `176:-4838400` `177:-2419200` `178:-3024000` `179:-241920` `180:-604800` `188:-1209600` `189:-241920` `196:-6451200` `197:-3628800` `198:-1612800`|
|9|12|`199:-322560` `201:-403200` `203:-201600` `205:-100800` `206:-100800` `224:-604800` `228:-1814400` `229:-362880` `232:-1814400` `254:-3628800` `255:-725760` `258:-3628800`|
|9|12|`261:-5443200` `262:-1088640` `265:-4838400` `266:-2419200` `270:-9072000` `271:-1814400` `278:-3628800` `284:-1814400` `285:-362880` `286:-18144000` `287:-7257600` `288:-3628800`|
|9|12|`289:-725760` `295:-16934400` `296:-3628800` `297:-3628800` `298:-2419200` `299:-1209600` `300:-604800` `301:-604800` `320:-1814400` `321:-362880` `335:-7257600` `336:-5443200`|
|9|12|`337:-362880` `339:-3628800` `343:-2419200` `344:-1209600` `345:-2419200` `346:-2419200` `347:-1209600` `348:-604800` `349:-604800` `378:-403200` `383:-1814400` `384:-362880`|
|9|12|`387:-3024000` `399:-4032000` `400:-1008000` `401:-1008000` `418:-2419200` `419:-604800` `438:-2419200` `439:-483840` `442:-2419200` `453:-4838400` `454:-3628800` `455:-241920`|
|9|12|`460:-9676800` `461:-4838400` `462:-3628800` `463:-2419200` `464:-2419200` `465:-1209600` `466:-604800` `467:-604800` `492:-4838400` `493:-1209600` `499:-3628800` `503:-2419200`|
|9|12|`504:-2419200` `505:-604800` `563:-302400` `573:-403200` `574:-201600` `594:-1209600` `595:-302400` `601:-604800` `602:-403200` `603:-403200` `604:-100800`|
|9|13|`7:10080` `9:100800` `12:312480` `15:312480` `46:-30240` `49:-40320` `51:-50400` `52:-10080` `58:-60480` `59:-120960` `67:-241920` `69:-161280`|
|9|13|`73:-90720` `75:-201600` `76:-40320` `89:-40320` `96:-50400` `97:-10080` `113:-120960` `114:-151200` `116:-100800` `117:-201600` `118:-40320` `119:-50400`|
|9|13|`122:-241920` `126:-181440` `127:-362880` `128:-1693440` `130:-604800` `131:-604800` `132:-120960` `135:-604800` `142:-241920` `145:-806400` `146:-161280` `148:-403200`|
|9|13|`151:-201600` `159:-1451520` `163:-967680` `166:-2419200` `167:-483840` `170:-151200` `185:-120960` `188:-806400` `189:-161280` `194:-100800` `195:-50400` `221:-241920`|
|9|13|`222:-362880` `224:-604800` `226:-725760` `228:-907200` `229:-181440` `232:-302400` `233:-604800` `234:-302400` `251:-362880` `252:-3265920` `254:-907200` `255:-181440`|
|9|13|`258:-2419200` `267:-1814400` `270:-7257600` `271:-1451520` `278:-1814400` `279:-2419200` `280:-1209600` `320:-907200` `321:-181440` `329:-302400` `330:-604800` `331:-302400`|
|9|13|`376:-161280` `378:-302400` `381:-604800` `383:-806400` `384:-161280` `387:-1814400` `388:-604800` `389:-302400` `438:-2419200` `439:-483840` `442:-2419200` `443:-2419200`|
|9|13|`444:-1209600` `561:-50400` `562:-10080` `563:-151200` `564:-100800` `565:-50400`|
|9|14|`9:10080` `10:2016` `13:151200` `14:30240` `17:302400` `18:60480` `22:100800` `51:-302400` `52:-60480` `53:-302400` `54:-60480` `58:-60480`|
|9|14|`63:-120960` `64:-151200` `65:-90720` `69:-40320` `70:-1411200` `71:-161280` `77:-100800` `78:-40320` `79:-352800` `80:-40320` `81:-201600` `82:-201600`|
|9|14|`83:-50400` `114:-302400` `117:-604800` `118:-120960` `131:-1814400` `132:-362880` `136:-2419200` `137:-241920` `140:-302400` `141:-302400` `145:-1209600` `146:-241920`|
|9|14|`147:-2419200` `152:-1209600` `153:-604800` `154:-120960` `161:-1814400` `162:-181440` `172:-1209600` `173:-120960` `175:-302400` `176:-2419200` `177:-1209600` `178:-1814400`|
|9|14|`179:-120960` `180:-453600` `181:-151200` `224:-604800` `228:-907200` `229:-181440` `232:-907200` `238:-2419200` `240:-604800` `241:-604800` `254:-907200` `255:-181440`|
|9|14|`270:-1814400` `271:-362880` `272:-1814400` `286:-3628800` `287:-1814400` `288:-907200` `289:-181440` `292:-1209600` `294:-604800` `295:-2419200` `296:-604800` `297:-604800`|
|9|14|`299:-1209600` `302:-302400` `303:-302400` `378:-100800` `387:-604800` `399:-806400` `400:-201600` `401:-201600` `411:-604800` `418:-1209600` `419:-302400` `425:-403200`|
|9|14|`426:-403200` `427:-100800`|
|9|15|`7:4032` `9:64512` `12:181440` `15:161280` `45:-181440` `46:-60480` `48:-161280` `49:-40320` `50:-262080` `51:-100800` `52:-20160` `58:-241920`|
|9|15|`59:-120960` `61:-302400` `62:-120960` `67:-241920` `69:-40320` `71:-40320` `73:-60480` `75:-403200` `76:-100800` `112:-362880` `115:-241920` `122:-483840`|
|9|15|`124:-604800` `126:-725760` `127:-362880` `128:-1693440` `130:-241920` `131:-3024000` `132:-483840` `135:-1209600` `142:-483840` `145:-1209600` `146:-241920` `159:-362880`|
|9|15|`163:-241920` `166:-1209600` `167:-241920` `169:-604800` `170:-302400` `171:-302400` `221:-604800` `224:-604800` `226:-725760` `228:-604800` `229:-120960` `232:-604800`|
|9|15|`251:-362880` `252:-1088640` `254:-1814400` `255:-362880` `258:-1209600` `259:-1209600` `260:-604800` `267:-725760` `270:-1814400` `271:-362880` `278:-1209600` `281:-1209600`|
|9|15|`282:-604800` `376:-40320` `381:-241920` `383:-403200` `384:-80640` `387:-604800` `390:-403200` `391:-201600`|
|9|16|`9:6048` `12:30240` `15:20160` `42:-50400` `45:-120960` `48:-120960` `50:-60480` `51:-151200` `52:-30240` `58:-120960` `59:-60480` `61:-151200`|
|9|16|`62:-60480` `112:-241920` `114:-151200` `115:-120960` `117:-201600` `118:-40320` `126:-181440` `128:-241920` `129:-302400` `131:-604800` `132:-120960` `134:-151200`|
|9|16|`135:-302400` `221:-60480` `226:-181440` `228:-302400` `229:-60480` `231:-201600` `232:-302400` `235:-100800`|
|10|10|`2:224` `8:6720` `11:20160` `15:-6720` `21:-3360` `66:-6720` `88:-2240` `142:-80640` `158:-60480` `185:-26880` `242:-100800` `267:-362880`|
|10|10|`308:-33600` `407:-302400`|
|10|11|`8:16800` `12:100800` `16:100800` `21:16800` `72:-235200` `84:-33600` `133:-504000` `144:-268800` `149:-268800` `157:-16800` `160:-604800` `182:-806400`|
|10|11|`183:-403200` `225:-100800` `253:-907200` `274:-1209600` `275:-604800` `290:-1814400` `291:-604800` `304:-806400` `305:-100800` `377:-67200` `382:-100800` `402:-268800`|
|10|11|`403:-134400` `440:-806400` `441:-403200` `472:-604800` `473:-806400` `474:-100800` `571:-67200` `572:-33600` `621:-100800` `622:-134400` `623:-16800`|
|10|12|`7:6720` `9:60480` `12:161280` `15:127680` `49:-26880` `53:-33600` `54:-6720` `63:-40320` `69:-201600` `79:-100800` `80:-20160` `113:-40320`|
|10|12|`117:-201600` `118:-26880` `120:-33600` `127:-362880` `130:-483840` `136:-403200` `137:-80640` `145:-403200` `146:-80640` `152:-201600` `155:-100800` `159:-725760`|
|10|12|`172:-1209600` `173:-241920` `222:-161280` `228:-604800` `229:-120960` `238:-403200` `239:-201600` `251:-241920` `252:-1088640` `261:-604800` `262:-120960` `270:-1814400`|
|10|12|`271:-362880` `284:-1814400` `285:-362880` `292:-1209600` `293:-604800` `376:-67200` `381:-120960` `383:-403200` `384:-80640` `395:-134400` `396:-26880` `397:-403200`|
|10|12|`398:-201600` `438:-1209600` `439:-241920` `456:-1209600` `457:-604800` `561:-33600` `562:-6720` `571:-67200` `572:-33600`|
|10|13|`7:3360` `8:23520` `11:40320` `47:-13440` `50:-20160` `57:-20160` `66:-33600` `73:-20160` `115:-80640` `119:-50400` `125:-181440` `128:-241920`|
|10|13|`142:-80640` `150:-50400` `158:-120960` `163:-241920` `226:-362880` `230:-302400` `250:-120960` `252:-362880` `267:-362880` `276:-302400` `381:-241920` `385:-302400`|
|10|14|`7:1344` `9:20160` `12:40320` `15:13440` `49:-67200` `53:-33600` `54:-13440` `63:-80640` `64:-100800` `65:-20160` `113:-40320` `127:-362880`|
|10|14|`136:-403200` `137:-80640` `138:-604800` `139:-80640` `140:-100800` `222:-80640` `236:-201600` `237:-40320` `251:-120960` `261:-604800` `262:-120960` `263:-403200`|
|10|14|`264:-201600` `376:-13440` `395:-134400` `396:-26880` `402:-134400` `403:-67200`|
|10|15|`6:224` `7:4032` `8:16800` `11:22400` `41:-2240` `42:-10080` `44:-20160` `45:-20160` `47:-26880` `50:-20160` `52:-10080` `57:-40320`|
|10|15|`59:-60480` `111:-26880` `112:-80640` `114:-50400` `121:-40320` `122:-80640` `125:-181440` `128:-241920` `130:-241920` `134:-151200` `220:-33600` `221:-120960`|
|10|15|`222:-40320` `223:-100800` `250:-60480` `252:-362880` `256:-302400` `375:-2240` `376:-26880` `377:-33600`|
|10|16|`7:672` `8:3360` `11:2240` `41:-2240` `44:-20160` `45:-20160` `46:-30240` `111:-13440` `113:-26880` `121:-20160` `122:-80640` `123:-50400`|
|10|16|`220:-6720` `222:-40320` `225:-33600`|
|11|11|`4:1050` `8:-6300` `12:-33600` `14:-8400` `16:-6300` `17:-16800` `19:-8400` `23:100800` `25:134400` `26:16800` `30:226800` `31:302400`|
|11|11|`35:100800` `43:-6300` `48:-201600` `51:-100800` `55:-75600` `56:-226800` `57:-302400` `58:-604800` `60:-75600` `67:-1209600` `71:-403200` `73:-100800`|
|11|11|`74:-201600` `75:-201600` `77:-100800` `78:-100800` `81:-100800` `85:-31500` `86:-75600` `87:-25200` `89:-67200` `92:-134400` `93:-403200` `96:-50400`|
|11|11|`97:-16800` `98:-67200` `101:-134400` `103:-67200` `104:-33600` `105:-67200` `108:-8400` `109:-42000` `110:-1050` `116:-100800` `121:-302400` `126:-2721600`|
|11|11|`128:-2419200` `129:-604800` `131:-1209600` `143:-1209600` `148:-806400` `150:-100800` `153:-201600` `159:-1814400` `161:-907200` `162:-2721600` `166:-3628800` `167:-1209600`|
|11|11|`169:-1209600` `170:-151200` `173:-1209600` `175:-604800` `178:-1209600` `179:-604800` `186:-100800` `188:-403200` `191:-2419200` `194:-201600` `195:-100800` `198:-1209600`|
|11|11|`201:-403200` `202:-1209600` `203:-403200` `206:-201600` `209:-67200` `213:-403200` `214:-134400` `217:-403200` `218:-33600` `219:-16800` `221:-604800` `227:-680400`|
|11|11|`230:-604800` `251:-604800` `254:-2721600` `255:-907200` `258:-1209600` `262:-1209600` `268:-907200` `270:-3628800` `273:-5443200` `280:-2419200` `288:-2721600` `293:-1209600`|
|11|11|`297:-1209600` `301:-1209600` `323:-1814400` `328:-680400` `334:-302400` `336:-1814400` `341:-2721600` `342:-907200` `348:-1209600` `352:-604800` `354:-302400` `358:-604800`|
|11|11|`378:-100800` `382:-100800` `383:-403200` `389:-604800` `398:-403200` `401:-403200` `441:-806400` `449:-302400` `454:-604800` `466:-1209600` `474:-201600` `477:-403200`|
|11|12|`12:16800` `14:8400` `17:268800` `19:134400` `23:604800` `24:302400` `28:268800` `51:-50400` `58:-604800` `59:-302400` `68:-1008000` `70:-1814400`|
|11|12|`72:-100800` `74:-453600` `75:-504000` `77:-100800` `81:-100800` `84:-50400` `90:-67200` `92:-672000` `94:-268800` `95:-1008000` `96:-16800` `98:-168000`|
|11|12|`100:-672000` `103:-201600` `105:-201600` `106:-33600` `107:-33600` `109:-8400` `122:-806400` `129:-2419200` `131:-2419200` `134:-1209600` `135:-1209600` `144:-806400`|
|11|12|`147:-1612800` `148:-403200` `149:-403200` `161:-3628800` `168:-4838400` `169:-3628800` `171:-1209600` `172:-2419200` `175:-1209600` `176:-2419200` `178:-2419200` `180:-604800`|
|11|12|`181:-604800` `190:-4032000` `191:-403200` `192:-1612800` `193:-403200` `196:-4838400` `200:-3225600` `201:-403200` `202:-403200` `203:-403200` `204:-403200` `208:-268800`|
|11|12|`211:-2150400` `212:-268800` `213:-134400` `214:-134400` `216:-537600` `217:-134400` `223:-806400` `224:-806400` `231:-604800` `233:-604800` `234:-302400` `235:-302400`|
|11|12|`259:-2419200` `260:-3628800` `261:-2419200` `264:-806400` `266:-1612800` `272:-5443200` `279:-2419200` `280:-1209600` `281:-2419200` `282:-1209600` `286:-3628800` `292:-2419200`|
|11|12|`294:-1209600` `295:-2419200` `298:-2419200` `299:-3628800` `300:-604800` `301:-604800` `302:-604800` `303:-604800` `322:-1814400` `327:-1814400` `330:-604800` `331:-302400`|
|11|12|`339:-7257600` `345:-2419200` `346:-2419200` `347:-1209600` `348:-604800` `349:-604800` `353:-1209600` `355:-604800` `356:-1209600` `357:-1209600` `358:-302400` `380:-268800`|
|11|12|`393:-806400` `394:-806400` `397:-806400` `399:-806400` `404:-268800` `405:-537600` `406:-806400` `440:-806400` `450:-806400` `463:-2419200` `468:-1612800` `469:-806400`|
|11|12|`473:-806400` `476:-1612800` `478:-268800` `479:-537600`|
|11|13|`8:6300` `12:100800` `16:226800` `21:100800` `60:-302400` `68:-403200` `72:-504000` `74:-100800` `84:-50400` `85:-6300` `90:-100800` `94:-403200`|
|11|13|`95:-201600` `123:-302400` `129:-1209600` `133:-907200` `160:-907200` `164:-1209600` `168:-1209600` `169:-604800` `174:-2419200` `175:-907200` `182:-1209600` `183:-604800`|
|11|13|`190:-806400` `192:-403200` `193:-201600` `194:-100800` `223:-302400` `225:-302400` `227:-226800` `253:-907200` `256:-907200` `263:-1209600` `264:-604800` `268:-907200`|
|11|13|`272:-1814400` `273:-907200` `290:-1814400` `291:-907200` `322:-1814400` `323:-907200` `326:-1360800` `327:-1814400` `328:-226800` `377:-100800` `385:-302400` `397:-403200`|
|11|13|`398:-201600` `402:-403200` `403:-201600` `456:-1209600` `457:-604800` `472:-604800` `473:-806400` `474:-100800`|
|11|14|`17:16800` `19:8400` `23:100800` `24:50400` `28:100800` `44:-100800` `48:-201600` `50:-100800` `51:-100800` `55:-25200` `58:-604800` `59:-151200`|
|11|14|`61:-151200` `68:-100800` `69:-201600` `70:-403200` `71:-403200` `74:-100800` `75:-403200` `76:-100800` `78:-100800` `80:-50400` `81:-100800` `82:-100800`|
|11|14|`83:-50400` `86:-50400` `87:-25200` `92:-134400` `95:-134400` `96:-16800` `98:-67200` `100:-134400` `103:-134400` `105:-67200` `106:-16800` `108:-8400`|
|11|14|`109:-16800` `112:-201600` `116:-100800` `119:-50400` `122:-403200` `127:-604800` `129:-604800` `131:-1209600` `132:-604800` `134:-151200` `135:-302400` `137:-403200`|
|11|14|`139:-201600` `145:-403200` `148:-806400` `151:-201600` `153:-201600` `155:-100800` `156:-100800` `161:-907200` `164:-302400` `169:-1209600` `171:-302400` `172:-604800`|
|11|14|`175:-604800` `178:-1209600` `180:-151200` `181:-151200` `193:-201600` `194:-201600` `195:-50400` `197:-403200` `201:-403200` `203:-403200` `204:-201600` `205:-100800`|
|11|14|`206:-100800` `207:-50400` `223:-100800` `224:-201600` `228:-302400` `234:-302400` `239:-201600` `241:-100800` `253:-302400` `260:-604800` `261:-604800` `264:-201600`|
|11|14|`266:-403200` `275:-604800` `283:-302400` `287:-907200` `291:-302400` `300:-604800` `302:-302400` `305:-100800` `307:-201600`|
|11|15|`12:16800` `16:100800` `21:100800` `46:-100800` `51:-252000` `60:-302400` `68:-201600` `70:-201600` `72:-100800` `74:-252000` `75:-100800` `77:-504000`|
|11|15|`81:-201600` `90:-67200` `92:-134400` `94:-268800` `95:-201600` `96:-16800` `98:-33600` `114:-201600` `116:-100800` `123:-403200` `129:-604800` `133:-302400`|
|11|15|`134:-604800` `138:-806400` `140:-403200` `143:-403200` `147:-806400` `148:-403200` `164:-604800` `174:-1814400` `175:-604800` `186:-100800` `190:-1612800` `191:-403200`|
|11|15|`194:-100800` `223:-201600` `225:-100800` `230:-302400` `238:-403200` `239:-201600` `256:-604800` `263:-806400` `264:-403200` `276:-302400` `290:-604800` `292:-1209600`|
|11|15|`293:-604800` `304:-806400` `305:-100800`|
|11|16|`16:6300` `21:16800` `42:-16800` `43:-6300` `46:-50400` `49:-67200` `51:-100800` `52:-50400` `54:-33600` `55:-50400` `60:-75600` `74:-100800`|
|11|16|`77:-100800` `81:-100800` `85:-12600` `86:-25200` `114:-50400` `117:-67200` `119:-50400` `120:-33600` `123:-50400` `133:-100800` `134:-151200` `138:-201600`|
|11|16|`140:-100800` `149:-134400` `150:-50400` `152:-201600` `155:-100800` `157:-16800`|
|12|12|`3:1680` `9:15120` `10:5040` `12:-26880` `13:378000` `14:75600` `15:-10080` `17:991200` `18:201600` `19:-8400` `22:705600` `61:-1058400`|
|12|12|`62:-483840` `64:-302400` `65:-60480` `70:-3225600` `71:-483840` `75:-1008000` `76:-161280` `77:-403200` `78:-80640` `79:-302400` `80:-20160` `81:-201600`|
|12|12|`82:-302400` `83:-40320` `86:-25200` `87:-25200` `92:-537600` `93:-107520` `96:-134400` `97:-26880` `99:-1612800` `100:-2016000` `101:-80640` `102:-201600`|
|12|12|`103:-470400` `104:-40320` `106:-16800` `108:-8400` `109:-8400` `124:-1209600` `131:-7257600` `132:-1451520` `135:-2419200` `136:-2419200` `137:-483840` `145:-806400`|
|12|12|`146:-161280` `161:-3628800` `162:-725760` `166:-4838400` `167:-967680` `170:-1209600` `172:-2419200` `173:-483840` `176:-19353600` `177:-6652800` `178:-8467200` `179:-846720`|
|12|12|`180:-1209600` `188:-806400` `189:-161280` `196:-8064000` `197:-2419200` `198:-1209600` `199:-241920` `200:-1612800` `201:-403200` `202:-403200` `208:-537600` `209:-403200`|
|12|12|`210:-26880` `211:-1075200` `212:-268800` `213:-134400` `214:-134400` `224:-2419200` `228:-1814400` `229:-362880` `232:-604800` `254:-3628800` `255:-725760` `258:-4838400`|
|12|12|`265:-4838400` `266:-2419200` `270:-7257600` `271:-1451520` `278:-2419200` `286:-14515200` `287:-7257600` `288:-3628800` `289:-725760` `295:-4838400` `296:-1209600` `297:-1209600`|
|12|12|`320:-907200` `321:-181440` `329:-302400` `335:-7257600` `336:-5443200` `337:-362880` `338:-5443200` `339:-7257600` `340:-3628800` `341:-907200` `342:-907200` `343:-2419200`|
|12|12|`344:-1209600` `353:-2419200` `354:-302400` `378:-403200` `387:-2419200` `399:-3225600` `400:-806400` `401:-806400` `442:-1209600` `460:-4838400` `461:-2419200` `475:-2419200`|
|12|12|`476:-3225600` `477:-403200`|
|12|13|`7:5040` `9:85680` `12:262080` `15:262080` `46:-30240` `49:-80640` `51:-151200` `52:-20160` `53:-100800` `54:-20160` `55:-25200` `58:-423360`|
|12|13|`59:-211680` `63:-241920` `67:-241920` `69:-403200` `70:-403200` `71:-80640` `75:-100800` `76:-20160` `79:-100800` `80:-20160` `81:-50400` `86:-25200`|
|12|13|`89:-80640` `92:-403200` `93:-80640` `122:-483840` `126:-725760` `127:-1451520` `128:-967680` `130:-967680` `131:-1209600` `132:-241920` `136:-1209600` `137:-241920`|
|12|13|`145:-403200` `146:-80640` `147:-806400` `148:-403200` `159:-1451520` `161:-1814400` `162:-362880` `163:-967680` `166:-1209600` `167:-241920` `172:-3024000` `173:-604800`|
|12|13|`188:-403200` `189:-80640` `190:-806400` `191:-403200` `221:-423360` `222:-241920` `226:-181440` `251:-725760` `252:-2177280` `254:-1814400` `255:-362880` `261:-3024000`|
|12|13|`262:-604800` `267:-725760` `270:-3628800` `271:-725760` `272:-1814400` `273:-907200` `284:-1814400` `285:-362880` `320:-1814400` `321:-362880` `322:-1814400` `323:-907200`|
|12|13|`376:-80640` `381:-241920` `383:-403200` `384:-80640` `395:-403200` `396:-80640` `438:-1209600` `439:-241920` `440:-806400` `441:-403200`|
|12|14|`13:33600` `14:6720` `17:201600` `18:40320` `22:201600` `45:-201600` `51:-403200` `52:-80640` `53:-100800` `54:-20160` `61:-604800` `62:-120960`|
|12|14|`70:-2016000` `71:-120960` `75:-403200` `76:-80640` `79:-403200` `80:-20160` `81:-302400` `82:-302400` `83:-60480` `92:-268800` `93:-26880` `95:-67200`|
|12|14|`99:-537600` `100:-672000` `101:-26880` `102:-33600` `103:-100800` `104:-6720` `105:-33600` `114:-403200` `117:-201600` `118:-40320` `124:-806400` `131:-1814400`|
|12|14|`132:-362880` `135:-1209600` `136:-1612800` `140:-403200` `141:-403200` `145:-806400` `146:-161280` `147:-1612800` `152:-806400` `153:-403200` `154:-80640` `166:-1209600`|
|12|14|`167:-241920` `168:-1209600` `176:-2419200` `177:-604800` `178:-1209600` `179:-120960` `188:-201600` `189:-40320` `190:-806400` `192:-403200` `196:-1612800` `197:-806400`|
|12|14|`198:-403200` `199:-80640` `200:-806400` `203:-201600` `204:-201600` `224:-604800` `232:-604800` `238:-806400` `240:-201600` `241:-201600` `258:-1209600` `265:-1612800`|
|12|14|`266:-806400` `274:-1209600` `278:-604800` `295:-2419200` `296:-604800` `297:-604800` `304:-806400` `306:-806400` `307:-201600`|
|12|15|`9:20160` `12:127680` `15:161280` `45:-181440` `48:-241920` `50:-100800` `58:-604800` `59:-181440` `61:-302400` `62:-60480` `63:-120960` `67:-322560`|
|12|15|`69:-120960` `70:-1209600` `71:-201600` `73:-100800` `75:-302400` `76:-40320` `79:-201600` `80:-40320` `89:-26880` `92:-134400` `93:-26880` `95:-67200`|
|12|15|`96:-33600` `97:-6720` `98:-33600` `112:-362880` `115:-40320` `122:-483840` `126:-725760` `127:-362880` `128:-1451520` `130:-120960` `131:-1209600` `132:-241920`|
|12|15|`136:-1209600` `137:-241920` `142:-161280` `145:-403200` `146:-80640` `159:-725760` `161:-1814400` `162:-362880` `163:-362880` `166:-1209600` `167:-241920` `168:-1209600`|
|12|15|`169:-604800` `172:-604800` `173:-120960` `185:-40320` `188:-403200` `189:-80640` `192:-403200` `193:-201600` `221:-241920` `222:-120960` `226:-362880` `228:-604800`|
|12|15|`229:-120960` `252:-725760` `261:-1209600` `262:-241920` `267:-362880` `270:-1814400` `271:-362880` `274:-1209600` `275:-604800`|
|12|16|`12:10080` `15:26880` `42:-36960` `45:-80640` `48:-80640` `50:-120960` `51:-100800` `52:-20160` `53:-134400` `54:-26880` `58:-60480` `67:-80640`|
|12|16|`68:-100800` `70:-201600` `71:-40320` `73:-10080` `74:-50400` `75:-100800` `76:-20160` `112:-80640` `115:-80640` `117:-134400` `118:-26880` `122:-80640`|
|12|16|`128:-241920` `136:-403200` `137:-80640` `142:-80640` `144:-268800` `145:-403200` `146:-80640` `149:-134400`|
|13|13|`2:420` `6:420` `7:8400` `8:39900` `11:67200` `12:-10080` `15:-10080` `16:-6300` `41:-6720` `42:-20160` `43:-6300` `44:-40320`|
|13|13|`47:-33600` `48:-40320` `50:-80640` `55:-25200` `56:-15120` `57:-221760` `66:-73920` `67:-161280` `73:-20160` `85:-6300` `88:-6720` `89:-40320`|
|13|13|`111:-60480` `112:-241920` `115:-40320` `116:-100800` `121:-141120` `125:-786240` `126:-181440` `128:-967680` `142:-161280` `143:-403200` `158:-302400` `159:-725760`|
|13|13|`163:-241920` `185:-40320` `186:-100800` `220:-100800` `221:-362880` `226:-181440` `227:-226800` `250:-302400` `251:-120960` `252:-1451520` `267:-725760` `268:-907200`|
|13|13|`375:-6720` `376:-80640` `381:-120960` `382:-100800`|
|13|14|`9:10080` `12:60480` `15:60480` `46:-60480` `49:-120960` `58:-241920` `63:-120960` `64:-151200` `65:-30240` `69:-161280` `70:-201600` `71:-40320`|
|13|14|`77:-403200` `78:-60480` `79:-50400` `80:-10080` `81:-50400` `113:-120960` `122:-241920` `126:-181440` `127:-362880` `130:-362880` `138:-1209600` `139:-241920`|
|13|14|`159:-362880` `161:-907200` `162:-181440` `172:-604800` `173:-120960` `174:-604800` `175:-302400` `221:-60480` `222:-120960` `236:-302400` `237:-60480` `252:-362880`|
|13|14|`261:-604800` `262:-120960` `284:-907200` `285:-181440` `290:-604800` `291:-302400`|
|13|15|`7:3360` `8:23520` `11:40320` `44:-120960` `45:-60480` `47:-53760` `48:-80640` `50:-20160` `51:-50400` `56:-60480` `57:-161280` `58:-120960`|
|13|15|`59:-60480` `66:-53760` `67:-80640` `69:-120960` `73:-20160` `74:-50400` `111:-40320` `121:-262080` `122:-241920` `125:-423360` `126:-362880` `127:-241920`|
|13|15|`128:-241920` `129:-302400` `130:-241920` `158:-120960` `159:-362880` `163:-241920` `164:-302400` `220:-60480` `222:-120960` `250:-120960` `251:-241920` `252:-362880`|
|13|15|`253:-302400`|
|13|16|`8:2520` `11:6720` `41:-6720` `44:-40320` `46:-30240` `47:-20160` `49:-40320` `56:-15120` `57:-20160` `58:-60480` `59:-30240` `60:-37800`|
|13|16|`111:-20160` `113:-40320` `121:-20160` `125:-60480` `127:-120960` `130:-120960` `133:-100800`|
|14|14|`3:336` `7:-1008` `9:-10080` `10:-1008` `12:-6720` `13:-8400` `14:-3360` `17:16800` `18:3360` `22:50400` `45:-201600` `51:-100800`|
|14|14|`52:-20160` `61:-151200` `62:-181440` `63:-201600` `65:-100800` `68:-100800` `70:-403200` `75:-201600` `76:-40320` `79:-100800` `81:-100800` `82:-100800`|
|14|14|`83:-20160` `90:-33600` `95:-134400` `96:-16800` `97:-3360` `100:-134400` `102:-33600` `103:-67200` `104:-6720` `106:-16800` `107:-16800` `114:-100800`|
|14|14|`124:-201600` `135:-604800` `136:-806400` `140:-201600` `141:-201600` `160:-302400` `170:-151200` `177:-604800` `180:-302400` `183:-403200` `184:-50400`|
|14|15|`12:13440` `15:40320` `45:-80640` `46:-80640` `48:-40320` `49:-40320` `58:-120960` `59:-120960` `64:-201600` `65:-40320` `67:-80640` `69:-80640`|
|14|15|`70:-201600` `71:-40320` `77:-201600` `78:-40320` `89:-13440` `92:-134400` `93:-26880` `94:-134400` `95:-67200` `112:-40320` `113:-40320` `122:-161280`|
|14|15|`128:-241920` `130:-241920` `136:-403200` `137:-80640` `138:-403200` `139:-80640` `163:-120960` `172:-604800` `173:-120960` `182:-403200` `183:-201600`|
|14|16|`15:3360` `42:-3360` `45:-20160` `46:-10080` `50:-20160` `53:-33600` `54:-6720` `59:-30240` `64:-50400` `65:-10080` `72:-33600` `73:-10080`|
|14|16|`79:-50400` `80:-10080` `84:-16800`|
|15|15|`2:224` `7:-672` `8:3360` `9:-4032` `11:20160` `12:-10080` `44:-120960` `47:-13440` `56:-60480` `57:-161280` `58:-241920` `63:-120960`|
|15|15|`66:-33600` `67:-80640` `68:-100800` `69:-80640` `88:-2240` `89:-26880` `90:-33600` `111:-40320` `121:-80640` `125:-241920` `127:-483840` `158:-60480`|
|15|15|`159:-362880` `160:-302400`|
|15|16|`11:2240` `41:-2240` `44:-20160` `47:-13440` `49:-26880` `57:-20160` `63:-40320` `66:-6720` `69:-40320` `72:-33600`|
|16|16|`2:28` `6:-28` `7:-672` `8:-840` `9:-1008`|

### Block 14695 — matching-occupancy

Order 1, 2 nonzero entries. Comparison result **PASS**. The coefficient SHA-256 is `810c6d79efe76268a48fa8948c08ebd2574052b42371978f747c19ead412276b`.

Actual label

```json
{"family":"matching-occupancy","caps":[6,9,3],"anchor":63,"lower":448,"upper":262080,"capacity":3,"occupied":2}
```

| Row $i$ | Column $j$ | Variable indices and coefficients |
|---:|---:|---|
|0|0|`11:840` `88:-840`|

### Block 14876 — coordinate-spectrum

Order 1, 2 nonzero entries. Comparison result **PASS**. The coefficient SHA-256 is `87eb768c43ac82580fd1caf7109c152d92088601af5edd23fc64b91e42c755c9`.

Actual label

```json
{"family":"coordinate-spectrum","inside":0,"outside":7,"capacity":13,"roots":[13],"sign":-1,"positive_scale_removed":1336608}
```

| Row $i$ | Column $j$ | Variable indices and coefficients |
|---:|---:|---|
|0|0|`0:5746` `1:-11`|

### Block 16249 — occupancy-matrix

Order 4, 2238 nonzero entries. Comparison result **PASS**. The coefficient SHA-256 is `6844492a76dd7efaf0b47c6ee8953db902c98cb6e9fbf67aa9414b9bdf61a812`.

Actual label

```json
{"family":"occupancy-matrix","core":2,"complement":false,"capacity":20,"polynomial":[0,-1,1],"kind":["gap",0],"classes":[0,1,2],"kept_rows":[0,1,2,3],"positive_scale_removed":15315300}
```

| Row $i$ | Column $j$ | Variable indices and coefficients |
|---:|---:|---|
|0|0|`4:9` `5:16`|
|0|1|`8:252` `11:1344` `12:4032` `15:12096` `16:7560` `17:10080` `18:2016` `21:16128` `22:18144` `23:24192` `24:12096` `27:2688`|
|0|1|`28:32256` `29:8064` `30:9072` `31:12096` `34:3456` `35:8640` `36:10368` `38:1008` `39:1296` `40:16`|
|0|2|`14:1008` `18:4032` `19:10080` `24:36288` `25:16128` `26:4032` `29:16128` `31:48384` `32:16128` `36:25920` `37:10368` `39:2592`|
|0|3|`4:18` `5:32` `20:630` `26:6048` `32:12096` `33:1344` `37:5184`|
|1|1|`8:252` `11:1344` `12:4032` `15:12096` `16:7560` `17:10080` `18:2016` `21:16128` `22:18144` `23:24192` `24:12096` `27:2688`|
|1|1|`28:32256` `29:8064` `30:9072` `31:12096` `34:3456` `35:8640` `36:10368` `38:1008` `39:1296` `40:16` `41:5376` `42:32256`|
|1|1|`43:7560` `44:64512` `45:145152` `46:96768` `47:16128` `48:96768` `49:129024` `50:48384` `51:241920` `52:48384` `53:120960` `54:48384`|
|1|1|`55:30240` `56:18144` `57:24192` `58:290304` `59:72576` `60:136080` `61:90720` `62:36288` `63:145152` `64:181440` `65:72576` `68:241920`|
|1|1|`69:96768` `70:483840` `71:96768` `72:322560` `74:120960` `77:241920` `78:48384` `79:120960` `80:24192` `81:120960` `82:60480` `83:24192`|
|1|1|`84:60480` `85:7560` `90:80640` `91:26880` `94:161280` `95:161280` `106:10080` `107:10080` `111:48384` `112:387072` `113:193536` `114:483840`|
|1|1|`115:96768` `116:362880` `117:645120` `118:129024` `119:241920` `120:120960` `121:120960` `122:774144` `123:302400` `124:362880` `125:145152` `126:653184`|
|1|1|`127:1161216` `128:580608` `129:2177280` `130:580608` `131:2903040` `132:580608` `133:1209600` `134:725760` `135:725760` `136:2903040` `137:580608` `138:1935360`|
|1|1|`139:387072` `140:725760` `141:362880` `143:483840` `144:1290240` `145:967680` `146:193536` `147:1935360` `148:967680` `149:967680` `150:120960` `152:967680`|
|1|1|`153:483840` `154:96768` `155:241920` `156:120960` `157:60480` `159:435456` `160:1451520` `161:2177280` `162:435456` `164:1451520` `165:967680` `168:2903040`|
|1|1|`169:1451520` `172:1451520` `173:290304` `174:2177280` `175:1451520` `176:1451520` `177:725760` `178:1451520` `179:145152` `180:362880` `181:181440` `182:1935360`|
|1|1|`183:1451520` `184:90720` `186:120960` `187:645120` `190:1935360` `192:1935360` `193:483840` `194:241920` `200:967680` `203:483840` `204:483840` `207:60480`|
|1|1|`212:161280` `215:80640` `220:56448` `221:725760` `222:580608` `223:1209600` `224:1935360` `225:564480` `226:435456` `227:1088640` `228:2903040` `229:580608`|
|1|1|`230:1451520` `231:2419200` `232:1451520` `233:2177280` `234:1451520` `235:967680` `236:725760` `237:145152` `238:3870720` `239:1451520` `240:967680` `241:967680`|
|1|1|`242:120960` `243:967680` `244:1935360` `245:483840` `246:241920` `247:241920` `248:241920` `249:241920` `250:72576` `251:725760` `252:1741824` `253:3628800`|
|1|1|`254:3265920` `255:653184` `256:2177280` `257:1128960` `258:1451520` `259:4354560` `260:2903040` `261:5806080` `262:1161216` `263:2419200` `264:1935360` `265:2903040`|
|1|1|`266:2903040` `268:2177280` `269:2903040` `270:4354560` `271:870912` `272:13063680` `273:4354560` `274:11612160` `275:4354560` `276:725760` `277:3870720` `279:5806080`|
|1|1|`281:5806080` `282:1451520` `283:725760` `284:2177280` `285:435456` `286:8709120` `287:4354560` `288:2177280` `289:435456` `290:7257600` `291:2903040` `292:5806080`|
|1|1|`293:1451520` `294:4354560` `295:5806080` `296:1451520` `297:1451520` `298:2903040` `299:5806080` `300:1451520` `301:1451520` `302:1451520` `303:1451520` `304:5806080`|
|1|1|`305:725760` `306:2903040` `307:1451520` `308:40320` `309:483840` `310:967680` `312:1935360` `313:241920` `315:967680` `318:967680` `319:483840` `322:8709120`|
|1|1|`323:2177280` `324:2177280` `325:5806080` `326:1632960` `327:4354560` `328:544320` `330:1451520` `332:1451520` `338:3265920` `339:8709120` `340:4354560` `341:1088640`|
|1|1|`342:1088640` `347:2903040` `350:1451520` `353:2903040` `354:362880` `355:725760` `356:1451520` `357:1451520` `358:725760` `359:725760` `361:241920` `365:483840`|
|1|1|`375:2688` `376:225792` `377:564480` `378:604800` `379:161280` `380:806400` `381:435456` `382:725760` `383:2419200` `384:483840` `385:2177280` `386:2257920`|
|1|1|`387:2903040` `388:5806080` `389:2177280` `390:4838400` `391:1935360` `392:1612800` `393:3870720` `394:2903040` `395:967680` `396:193536` `397:4838400` `398:1935360`|
|1|1|`399:7741440` `400:1935360` `401:1935360` `402:2257920` `403:967680` `404:1612800` `405:1935360` `406:2580480` `407:362880` `408:1935360` `409:5806080` `410:1451520`|
|1|1|`411:7257600` `412:4354560` `413:5806080` `414:725760` `415:967680` `416:967680` `417:1935360` `418:2903040` `419:725760` `420:2177280` `421:5806080` `422:2177280`|
|1|1|`423:725760` `424:725760` `425:3870720` `426:3870720` `427:967680` `428:1935360` `429:1935360` `430:1935360` `431:1935360` `432:967680` `433:967680` `434:53760`|
|1|1|`435:322560` `437:322560` `438:2903040` `439:580608` `440:9676800` `441:3870720` `443:11612160` `444:2903040` `445:5806080` `446:7741440` `447:4354560` `448:5806080`|
|1|1|`449:725760` `450:1935360` `451:1935360` `453:2903040` `454:2177280` `455:145152` `456:8709120` `457:2903040` `458:3386880` `459:4838400` `460:5806080` `461:2903040`|
|1|1|`462:13063680` `463:11612160` `464:11612160` `465:8709120` `466:2903040` `467:2903040` `468:3870720` `469:3870720` `470:1451520` `471:967680` `472:4354560` `473:7741440`|
|1|1|`474:967680` `475:5806080` `476:11612160` `477:1451520` `478:967680` `479:1935360` `480:2903040` `482:11612160` `483:5806080` `484:13063680` `485:5806080` `486:725760`|
|1|1|`487:2903040` `488:2903040` `489:322560` `490:967680` `494:8709120` `496:2903040` `499:17418240` `500:5806080` `501:2903040` `502:2903040` `503:5806080` `504:5806080`|
|1|1|`505:1451520` `506:8709120` `507:2903040` `508:5806080` `509:5806080` `510:5806080` `511:1451520` `512:5806080` `514:967680` `516:1935360` `519:5806080` `520:2903040`|
|1|1|`521:1451520` `522:725760` `523:967680` `524:1935360` `525:1935360` `526:967680` `527:483840` `529:161280` `533:2177280` `536:4354560` `537:1451520` `539:4354560`|
|1|1|`544:2177280` `545:1451520` `546:2903040` `557:241920` `561:141120` `562:28224` `563:544320` `564:1451520` `565:604800` `566:2177280` `567:846720` `568:1209600`|
|1|1|`569:680400` `570:1088640` `571:1128960` `572:483840` `573:1209600` `574:604800` `575:362880` `576:564480` `577:1814400` `578:1935360` `579:1209600` `580:4354560`|
|1|1|`581:1451520` `582:3386880` `583:2419200` `584:4354560` `585:3870720` `586:483840` `587:1451520` `588:5806080` `589:2903040` `590:6531840` `591:2903040` `592:1451520`|
|1|1|`593:1451520` `594:2903040` `595:725760` `596:4354560` `597:4354560` `598:1451520` `599:544320` `600:544320` `601:7257600` `602:3870720` `603:3870720` `604:967680`|
|1|1|`605:3628800` `606:2177280` `607:8709120` `608:2903040` `609:2903040` `610:6531840` `611:6531840` `612:3628800` `613:2903040` `614:2177280` `615:2903040` `616:3265920`|
|1|1|`617:2903040` `618:2903040` `619:2903040` `620:1935360` `621:967680` `622:967680` `623:120960` `624:3628800` `625:967680` `626:4354560` `627:5806080` `628:1451520`|
|1|1|`629:2903040` `630:1935360` `631:967680` `632:2177280` `633:725760` `634:2177280` `635:725760` `636:241920` `637:1088640` `638:362880` `639:4354560` `640:2177280`|
|1|1|`642:725760` `643:4354560` `644:1451520` `645:2177280` `646:725760` `647:1451520` `648:725760` `649:272160` `650:725760` `651:725760` `652:4354560` `653:2903040`|
|1|1|`654:362880` `655:967680` `656:7257600` `657:1935360` `658:2177280` `660:4354560` `661:1451520` `662:241920` `664:8709120` `665:5806080` `666:725760` `667:1451520`|
|1|1|`668:11612160` `669:2903040` `670:2177280` `671:8709120` `672:13063680` `673:1451520` `674:5806080` `675:2903040` `676:725760` `677:4354560` `679:4354560` `680:8709120`|
|1|1|`681:2903040` `683:8709120` `686:4354560` `687:1451520` `688:2903040` `689:8709120` `690:8709120` `691:1451520` `692:2903040` `693:4354560` `694:3265920` `695:3265920`|
|1|1|`696:1451520` `697:2903040` `698:2177280` `699:2177280` `700:1451520` `701:2903040` `702:2177280` `703:4354560` `704:1451520` `705:362880` `706:2177280` `709:2177280`|
|1|1|`712:4354560` `713:4354560` `716:725760` `717:725760` `720:1088640` `723:2177280` `726:1451520` `735:193536` `736:193536` `737:24192` `738:870912` `739:580608`|
|1|1|`740:387072` `741:1451520` `742:387072` `743:435456` `744:870912` `745:145152` `746:435456` `747:483840` `748:60480` `749:2612736` `750:1161216` `751:677376`|
|1|1|`752:1451520` `753:1451520` `754:774144` `755:544320` `756:1306368` `757:1741824` `758:870912` `759:435456` `760:1306368` `761:1741824` `762:653184` `763:161280`|
|1|1|`764:677376` `765:193536` `766:870912` `767:2322432` `768:1451520` `769:387072` `770:774144` `771:580608` `772:362880` `773:435456` `774:1741824` `775:2612736`|
|1|1|`776:580608` `777:1161216` `778:129024` `779:870912` `780:870912` `781:2322432` `782:3483648` `783:580608` `784:5225472` `785:3483648` `786:3483648` `787:1161216`|
|1|1|`788:1741824` `789:2903040` `790:1741824` `791:3483648` `792:2612736` `793:2612736` `794:580608` `795:1161216` `796:1741824` `797:3483648` `798:580608` `799:870912`|
|1|1|`800:870912` `802:1741824` `803:1741824` `804:1741824` `805:290304` `806:870912` `807:1741824` `808:3483648` `809:1741824` `810:653184` `811:580608` `812:580608`|
|1|1|`813:2322432` `814:1161216` `815:387072` `816:870912` `817:1741824` `818:1741824` `819:870912` `820:870912` `821:1741824` `822:145152` `823:870912` `824:435456`|
|1|1|`825:435456` `827:435456` `828:1741824` `830:3483648` `831:2612736` `832:1741824` `833:3483648` `834:580608` `835:435456` `836:290304` `837:1741824` `838:1741824`|
|1|1|`839:290304` `842:1741824` `846:870912` `847:1741824` `848:1741824` `849:3483648` `850:580608` `856:870912` `857:435456` `858:870912` `863:2688` `864:64512`|
|1|1|`865:48384` `866:5376` `867:145152` `868:145152` `869:290304` `870:290304` `871:32256` `872:483840` `873:290304` `874:193536` `875:387072` `876:64512`|
|1|1|`877:290304` `878:580608` `879:870912` `880:580608` `881:145152` `882:290304` `883:290304` `884:9072` `885:290304` `886:217728` `887:72576` `888:24192`|
|1|1|`889:870912` `890:290304` `891:1161216` `892:290304` `893:217728` `894:435456` `895:435456` `896:290304` `897:580608` `898:580608` `899:387072` `900:241920`|
|1|1|`901:64512` `902:193536` `903:290304` `904:290304` `905:580608` `906:580608` `907:290304` `908:290304` `910:1161216` `911:580608` `913:580608` `914:290304`|
|1|1|`915:580608` `916:290304` `917:580608` `918:24192` `919:145152` `920:145152` `927:20736` `928:41472` `929:10368` `930:62208` `931:82944` `932:10368`|
|1|1|`933:82944` `934:82944` `935:82944` `936:20736` `938:82944` `939:41472`|
|1|2|`47:16128` `50:72576` `52:24192` `54:8064` `55:30240` `57:48384` `59:108864` `62:36288` `66:16128` `67:145152` `69:145152` `71:48384`|
|1|2|`73:24192` `74:120960` `75:241920` `76:72576` `79:120960` `80:24192` `81:60480` `84:10080` `86:30240` `89:16128` `90:40320` `92:161280`|
|1|2|`93:32256` `95:80640` `98:40320` `100:80640` `102:20160` `103:40320` `104:4032` `115:96768` `118:32256` `119:181440` `120:40320` `125:181440`|
|1|2|`128:1161216` `130:580608` `132:290304` `134:544320` `135:1088640` `137:96768` `140:120960` `142:193536` `143:725760` `145:1451520` `146:387072` `148:967680`|
|1|2|`149:322560` `150:120960` `151:483840` `152:483840` `155:362880` `158:72576` `159:870912` `162:108864` `163:580608` `164:1088640` `166:2177280` `167:580608`|
|1|2|`169:1451520` `170:181440` `171:725760` `172:2177280` `173:435456` `175:362880` `177:725760` `178:725760` `180:362880` `183:241920` `186:120960` `187:161280`|
|1|2|`188:483840` `189:96768` `190:1935360` `191:1451520` `192:967680` `193:725760` `195:60480` `196:1935360` `197:1451520` `198:483840` `199:96768` `200:483840`|
|1|2|`201:483840` `203:241920` `205:120960` `206:120960` `211:645120` `212:161280` `213:80640` `214:80640` `226:544320` `229:145152` `230:1451520` `232:1451520`|
|1|2|`234:1088640` `235:362880` `239:483840` `240:241920` `242:241920` `243:483840` `244:1451520` `245:1209600` `246:80640` `248:241920` `250:72576` `252:2177280`|
|1|2|`255:217728` `256:1451520` `258:2903040` `260:2177280` `262:145152` `264:241920` `266:483840` `267:1306368` `268:2177280` `270:8709120` `271:2177280` `273:2177280`|
|1|2|`275:1451520` `276:1451520` `277:1451520` `278:2903040` `279:8709120` `280:7257600` `281:4354560` `282:3628800` `283:725760` `284:2177280` `285:435456` `287:2177280`|
|1|2|`288:544320` `291:362880` `292:4354560` `293:2903040` `294:1451520` `295:8709120` `296:2903040` `297:2177280` `298:2903040` `299:1451520` `300:2177280` `302:725760`|
|1|2|`304:1935360` `306:483840` `307:241920` `309:161280` `310:967680` `311:725760` `312:3870720` `313:241920` `314:3870720` `315:967680` `316:967680` `317:483840`|
|1|2|`319:241920` `320:2177280` `321:435456` `322:6531840` `323:4354560` `325:1451520` `327:2177280` `330:1451520` `331:1088640` `332:1451520` `333:2903040` `334:362880`|
|1|2|`335:4354560` `336:4354560` `337:217728` `339:8709120` `340:2177280` `341:1088640` `343:1451520` `344:725760` `345:5806080` `346:8709120` `347:4354560` `348:2177280`|
|1|2|`349:1451520` `352:362880` `353:4354560` `354:362880` `356:1451520` `359:241920` `360:483840` `364:1451520` `365:241920` `366:967680` `367:967680` `368:241920`|
|1|2|`381:435456` `384:96768` `385:1814400` `387:3628800` `389:1088640` `391:725760` `392:645120` `393:3870720` `394:1935360` `396:16128` `398:483840` `400:483840`|
|1|2|`401:241920` `403:80640` `404:322560` `405:645120` `406:322560` `407:1088640` `408:967680` `409:5806080` `410:4354560` `412:2177280` `413:2903040` `414:241920`|
|1|2|`415:645120` `416:1935360` `417:5806080` `418:2903040` `419:725760` `421:4354560` `422:362880` `423:725760` `425:1451520` `426:483840` `427:241920` `428:2419200`|
|1|2|`429:4838400` `430:1451520` `432:322560` `433:322560` `435:645120` `436:645120` `437:322560` `438:3628800` `439:870912` `441:967680` `442:4354560` `443:11612160`|
|1|2|`444:8709120` `446:2903040` `448:4354560` `449:362880` `450:3870720` `451:2903040` `452:5806080` `454:725760` `456:5806080` `457:3628800` `459:967680` `460:11612160`|
|1|2|`461:7257600` `463:8709120` `464:5806080` `465:2903040` `466:2177280` `468:7741440` `469:3870720` `473:1935360` `476:3870720` `477:483840` `479:645120` `480:5806080`|
|1|2|`481:3628800` `482:4354560` `483:7257600` `484:6531840` `485:14515200` `486:1088640` `487:1451520` `488:967680` `490:1935360` `491:3870720` `492:5806080` `493:1451520`|
|1|2|`494:13063680` `495:29030400` `496:5806080` `497:3628800` `498:2177280` `499:13063680` `500:1451520` `501:5806080` `503:14515200` `504:11612160` `505:2903040` `506:4354560`|
|1|2|`507:2903040` `508:11612160` `509:8709120` `511:725760` `512:2903040` `513:2903040` `514:967680` `515:9676800` `516:1935360` `517:3870720` `518:3870720` `519:4354560`|
|1|2|`523:967680` `524:2903040` `525:967680` `528:322560` `530:645120` `533:2177280` `534:4354560` `535:362880` `536:4354560` `537:2903040` `538:8709120` `539:8709120`|
|1|2|`540:1451520` `541:11612160` `542:2903040` `543:725760` `545:725760` `546:1451520` `547:483840` `549:967680` `551:1451520` `552:1451520` `554:483840` `555:967680`|
|1|2|`562:4032` `563:544320` `565:181440` `566:1814400` `568:483840` `570:725760` `572:80640` `574:120960` `576:80640` `578:483840` `579:241920` `580:3628800`|
|1|2|`581:2540160` `583:725760` `585:1451520` `586:120960` `587:4354560` `588:2903040` `589:5080320` `590:4354560` `591:10160640` `592:1451520` `593:967680` `594:3628800`|
|1|2|`595:1088640` `597:2177280` `598:362880` `599:362880` `602:1451520` `603:967680` `604:362880` `605:1451520` `606:2177280` `607:8709120` `608:8709120` `609:3628800`|
|1|2|`610:4354560` `611:3265920` `612:725760` `613:3628800` `615:3628800` `616:1088640` `617:3628800` `618:2903040` `620:967680` `622:161280` `624:1088640` `625:241920`|
|1|2|`626:2177280` `627:1451520` `628:483840` `629:1451520` `630:725760` `632:1451520` `633:2177280` `634:3265920` `635:1814400` `637:2177280` `638:725760` `639:13063680`|
|1|2|`640:1451520` `641:5080320` `642:1088640` `643:4354560` `644:3628800` `646:1451520` `647:1451520` `648:241920` `650:725760` `652:2177280` `653:4354560` `654:362880`|
|1|2|`656:2177280` `657:483840` `658:2177280` `659:4354560` `660:6531840` `661:3628800` `663:2903040` `664:8709120` `665:17418240` `666:2177280` `668:4354560` `669:1935360`|
|1|2|`671:6531840` `672:6531840` `673:725760` `674:5806080` `676:181440` `677:8709120` `678:7257600` `679:6531840` `680:10886400` `681:8709120` `682:13063680` `683:21772800`|
|1|2|`684:17418240` `685:4354560` `687:2177280` `688:2903040` `689:21772800` `690:10886400` `691:2177280` `692:7257600` `693:2177280` `694:2177280` `696:725760` `697:725760`|
|1|2|`698:4354560` `700:2903040` `701:1451520` `702:725760` `703:2903040` `704:241920` `706:4354560` `707:1451520` `708:8709120` `709:3265920` `710:2177280` `711:13063680`|
|1|2|`712:10886400` `713:4354560` `714:2177280` `715:4354560` `716:725760` `718:10886400` `719:4354560` `720:1088640` `721:2903040` `722:1451520` `723:1088640` `724:1451520`|
|1|2|`726:725760` `727:1088640` `730:1088640` `731:2177280` `732:362880` `736:48384` `737:4032` `738:544320` `739:1161216` `741:580608` `742:193536` `743:653184`|
|1|2|`744:1741824` `745:145152` `747:145152` `748:24192` `749:2177280` `750:2032128` `752:580608` `753:435456` `754:483840` `756:870912` `757:653184` `758:725760`|
|1|2|`760:1306368` `761:2612736` `764:145152` `765:32256` `767:870912` `768:435456` `769:96768` `770:387072` `771:145152` `772:145152` `773:290304` `774:580608`|
|1|2|`775:3048192` `776:870912` `777:580608` `779:3048192` `780:2612736` `781:1161216` `782:6096384` `783:2322432` `784:3483648` `785:12192768` `786:5225472` `787:4064256`|
|1|2|`788:870912` `789:870912` `790:1741824` `791:5225472` `792:2177280` `794:580608` `795:580608` `796:1161216` `797:2903040` `798:290304` `799:2612736` `800:2177280`|
|1|2|`801:6967296` `802:6096384` `803:2612736` `804:2612736` `805:290304` `806:435456` `807:2612736` `808:4354560` `810:217728` `811:725760` `812:145152` `813:1451520`|
|1|2|`814:290304` `815:96768` `816:2177280` `817:2612736` `818:5225472` `819:653184` `820:870912` `821:435456` `823:435456` `824:145152` `825:870912` `826:870912`|
|1|2|`827:435456` `828:1741824` `829:5225472` `830:5225472` `831:2612736` `832:4354560` `833:12192768` `834:1161216` `837:1741824` `838:870912` `840:4354560` `841:1741824`|
|1|2|`842:5225472` `843:12192768` `844:1741824` `845:1161216` `847:3483648` `848:5225472` `849:3483648` `850:580608` `851:1741824` `852:1741824` `853:5225472` `854:5225472`|
|1|2|`856:1306368` `858:870912` `859:435456` `860:870912` `864:16128` `865:16128` `867:145152` `868:290304` `869:181440` `870:580608` `872:193536` `873:96768`|
|1|2|`874:96768` `875:290304` `876:32256` `877:435456` `878:1161216` `879:1016064` `880:2612736` `882:435456` `883:217728` `885:145152` `886:217728` `888:24192`|
|1|2|`889:1016064` `890:435456` `891:2032128` `892:1161216` `893:145152` `894:435456` `896:362880` `897:1016064` `898:290304` `899:193536` `900:145152` `902:96768`|
|1|2|`903:290304` `904:193536` `905:2032128` `906:1161216` `907:870912` `908:870912` `909:1741824` `910:4644864` `911:870912` `912:2322432` `913:725760` `914:725760`|
|1|2|`915:290304` `916:725760` `917:870912` `919:145152` `921:1741824` `922:1161216` `923:1161216` `924:580608` `926:290304` `927:20736` `928:82944` `929:10368`|
|1|2|`930:82944` `931:414720` `933:145152` `934:186624` `935:62208` `936:31104` `937:414720` `938:290304` `939:82944` `940:82944` `941:10368`|
|1|3|`8:504` `11:2688` `12:8064` `15:24192` `16:15120` `17:20160` `18:4032` `21:32256` `22:36288` `23:48384` `24:24192` `27:5376`|
|1|3|`28:64512` `29:16128` `30:18144` `31:24192` `34:6912` `35:17280` `36:20736` `38:2016` `39:2592` `40:32` `66:12096` `73:36288`|
|1|3|`85:11340` `88:1344` `89:24192` `96:30240` `97:6048` `142:145152` `150:181440` `157:30240` `158:36288` `163:435456` `170:272160` `185:72576`|
|1|3|`186:181440` `188:725760` `189:145152` `194:362880` `195:181440` `208:161280` `209:120960` `210:8064` `242:181440` `246:120960` `267:653184` `276:1088640`|
|1|3|`278:2177280` `283:1088640` `305:362880` `308:181440` `309:241920` `310:1451520` `311:725760` `313:362880` `320:1632960` `321:326592` `328:408240` `329:544320`|
|1|3|`330:2177280` `331:1088640` `334:544320` `343:2177280` `344:1088640` `355:1088640` `360:1451520` `361:362880` `362:181440` `363:181440` `407:544320` `414:362880`|
|1|3|`415:483840` `416:1451520` `432:483840` `434:80640` `435:483840` `442:2177280` `449:544320` `450:2903040` `474:362880` `478:483840` `480:4354560` `481:2177280`|
|1|3|`486:1088640` `487:2177280` `488:1451520` `489:1451520` `490:1451520` `491:2903040` `492:4354560` `493:1088640` `507:4354560` `513:2903040` `514:1451520` `521:1088640`|
|1|3|`523:1451520` `527:241920` `528:483840` `529:241920` `532:725760` `533:3265920` `534:4354560` `535:544320` `540:2177280` `545:2177280` `547:1451520` `548:725760`|
|1|3|`586:181440` `587:2177280` `592:1088640` `593:725760` `623:30240` `625:362880` `628:725760` `631:362880` `632:1088640` `633:1088640` `635:1088640` `636:362880`|
|1|3|`637:1632960` `638:544320` `646:1088640` `648:362880` `649:136080` `654:272160` `657:725760` `658:1632960` `659:2177280` `661:2177280` `662:544320` `663:1451520`|
|1|3|`673:1088640` `677:6531840` `678:4354560` `687:2177280` `691:2177280` `693:3265920` `696:1088640` `702:1088640` `704:725760` `705:1632960` `706:3265920` `707:1088640`|
|1|3|`708:6531840` `709:3265920` `710:2177280` `714:2177280` `716:1088640` `717:1088640` `720:1632960` `727:1632960` `728:544320` `729:181440` `737:6048` `742:145152`|
|1|3|`743:326592` `765:48384` `769:145152` `771:217728` `776:435456` `777:435456` `778:145152` `779:1306368` `780:1306368` `790:1306368` `794:435456` `796:870912`|
|1|3|`798:290304` `799:1306368` `800:1306368` `805:217728` `810:326592` `814:435456` `815:145152` `816:1306368` `819:653184` `822:217728` `823:435456` `824:217728`|
|1|3|`825:653184` `828:1306368` `829:2612736` `832:2612736` `834:870912` `836:653184` `837:1306368` `839:217728` `840:2612736` `847:2612736` `851:1306368` `852:1306368`|
|1|3|`856:1306368` `857:653184` `859:653184` `865:12096` `866:4032` `868:108864` `874:72576` `876:24192` `877:217728` `882:217728` `890:217728` `893:108864`|
|1|3|`896:217728` `899:145152` `901:72576` `902:72576` `904:145152` `905:870912` `907:435456` `908:435456` `913:435456` `914:435456` `916:435456` `918:36288`|
|1|3|`919:108864` `920:108864` `921:870912` `924:435456` `928:31104` `929:5184` `933:62208` `936:15552` `938:124416` `940:41472`|
|2|2|`14:1008` `18:4032` `19:10080` `24:36288` `25:16128` `26:4032` `29:16128` `31:48384` `32:16128` `36:25920` `37:10368` `39:2592`|
|2|2|`76:48384` `80:24192` `83:72576` `86:30240` `87:30240` `93:32256` `97:8064` `98:40320` `101:32256` `102:20160` `103:80640` `104:8064`|
|2|2|`105:120960` `106:10080` `107:30240` `108:10080` `109:10080` `146:193536` `151:241920` `154:193536` `155:241920` `156:483840` `167:290304` `171:362880`|
|2|2|`173:290304` `179:580608` `180:362880` `181:725760` `184:30240` `189:193536` `191:967680` `193:483840` `195:241920` `197:967680` `198:483840` `199:193536`|
|2|2|`201:483840` `202:1935360` `203:483840` `204:1935360` `205:362880` `206:483840` `207:120960` `209:161280` `211:645120` `212:161280` `213:483840` `214:322560`|
|2|2|`215:161280` `216:645120` `217:322560` `218:20160` `219:20160` `245:483840` `248:120960` `249:604800` `271:870912` `280:2903040` `282:1451520` `285:217728`|
|2|2|`289:1088640` `293:1451520` `296:1451520` `297:1451520` `300:1451520` `301:7257600` `302:725760` `303:3628800` `307:483840` `311:967680` `312:1935360` `313:483840`|
|2|2|`314:1935360` `315:967680` `316:1451520` `317:1451520` `318:4838400` `319:1451520` `321:435456` `323:2177280` `331:1451520` `332:725760` `333:1451520` `334:725760`|
|2|2|`336:2177280` `337:435456` `341:1088640` `342:5443200` `344:1451520` `345:2903040` `346:5806080` `347:2903040` `348:4354560` `349:4354560` `350:3628800` `351:7257600`|
|2|2|`352:2177280` `353:2903040` `354:1088640` `356:1451520` `357:7257600` `358:2177280` `359:241920` `360:1935360` `361:483840` `362:604800` `363:120960` `364:2903040`|
|2|2|`365:1451520` `366:2903040` `367:5806080` `368:967680` `369:967680` `370:483840` `371:161280` `372:322560` `374:161280` `410:1451520` `417:1935360` `419:362880`|
|2|2|`423:362880` `424:2177280` `427:241920` `428:967680` `429:1935360` `430:967680` `431:5806080` `433:322560` `435:645120` `436:645120` `437:1290240` `439:290304`|
|2|2|`444:2903040` `451:967680` `452:1935360` `455:290304` `457:1451520` `461:2903040` `466:1451520` `467:8709120` `468:3870720` `469:1935360` `470:2903040` `471:5806080`|
|2|2|`477:483840` `479:645120` `481:2903040` `483:2903040` `485:5806080` `486:1451520` `490:1935360` `491:3870720` `493:1451520` `495:11612160` `496:2903040` `497:4354560`|
|2|2|`498:5806080` `501:2903040` `502:8709120` `503:5806080` `504:5806080` `505:4354560` `508:5806080` `509:5806080` `510:34836480` `511:5806080` `512:1935360` `513:3870720`|
|2|2|`514:1935360` `515:11612160` `516:7741440` `517:5806080` `518:15482880` `519:2903040` `520:8709120` `522:1451520` `523:1935360` `524:5806080` `525:7741440` `526:3870720`|
|2|2|`527:322560` `528:1612800` `529:322560` `530:1935360` `531:645120` `534:5806080` `535:1814400` `536:2177280` `537:4354560` `538:4354560` `539:4354560` `540:2903040`|
|2|2|`541:17418240` `542:11612160` `543:2177280` `544:6531840` `545:2903040` `546:11612160` `547:4838400` `548:967680` `549:2903040` `550:4354560` `551:2903040` `552:11612160`|
|2|2|`553:725760` `554:2903040` `555:5806080` `556:2903040` `557:725760` `558:322560` `559:322560` `581:725760` `589:1451520` `591:2903040` `595:362880` `599:181440`|
|2|2|`600:1270080` `604:241920` `606:725760` `608:2903040` `609:1451520` `613:1451520` `614:5080320` `615:1451520` `617:1451520` `618:1451520` `619:10160640` `620:483840`|
|2|2|`629:967680` `630:483840` `633:1451520` `634:1088640` `635:1451520` `638:725760` `639:4354560` `640:725760` `641:4354560` `642:3628800` `643:2177280` `644:4354560`|
|2|2|`645:7620480` `646:1451520` `647:7257600` `648:483840` `650:1088640` `651:1814400` `653:1451520` `654:362880` `659:2903040` `660:2177280` `661:2903040` `665:5806080`|
|2|2|`666:2177280` `669:967680` `674:2903040` `675:10160640` `676:1814400` `678:5806080` `679:2177280` `680:4354560` `681:8709120` `682:4354560` `683:8709120` `684:17418240`|
|2|2|`685:14515200` `686:15240960` `687:2903040` `688:14515200` `689:8709120` `690:4354560` `691:2903040` `692:8709120` `694:1088640` `695:7620480` `696:1451520` `697:7257600`|
|2|2|`698:2177280` `699:15240960` `700:4354560` `701:14515200` `703:1451520` `704:967680` `706:4354560` `707:3628800` `708:8709120` `709:4354560` `710:7257600` `711:13063680`|
|2|2|`712:13063680` `713:21772800` `714:7257600` `715:11612160` `716:3628800` `717:2177280` `718:13063680` `719:21772800` `720:2177280` `721:8709120` `722:5806080` `723:10886400`|
|2|2|`724:5806080` `725:5806080` `726:5806080` `727:5443200` `728:1088640` `729:1088640` `730:3265920` `731:8709120` `732:2903040` `733:2177280` `734:1451520` `739:290304`|
|2|2|`744:435456` `748:12096` `750:580608` `754:193536` `758:290304` `759:1161216` `760:435456` `761:870912` `762:1741824` `770:193536` `772:72576` `773:145152`|
|2|2|`775:870912` `776:580608` `779:1741824` `782:1741824` `783:1741824` `785:3483648` `786:1741824` `787:3483648` `791:1741824` `792:870912` `793:6967296` `794:580608`|
|2|2|`795:3483648` `797:1161216` `798:387072` `799:1741824` `800:1741824` `801:5225472` `802:5225472` `803:10450944` `804:870912` `805:290304` `807:870912` `808:1741824`|
|2|2|`809:6967296` `811:870912` `812:1741824` `813:580608` `815:193536` `816:1741824` `817:2612736` `818:5225472` `819:870912` `820:5225472` `821:5225472` `823:580608`|
|2|2|`824:290304` `825:2177280` `826:1741824` `827:2177280` `829:3483648` `830:1741824` `831:870912` `832:3483648` `833:10450944` `834:2903040` `835:870912` `837:1741824`|
|2|2|`838:5225472` `839:580608` `840:8709120` `841:4354560` `842:1741824` `843:10450944` `844:6967296` `845:2322432` `846:3483648` `847:3483648` `848:5225472` `849:20901888`|
|2|2|`850:2903040` `851:4354560` `852:4354560` `853:10450944` `854:17418240` `855:4354560` `856:4354560` `857:1741824` `858:4354560` `859:3919104` `860:4354560` `861:870912`|
|2|2|`862:217728` `867:36288` `868:145152` `870:145152` `875:96768` `876:32256` `878:290304` `879:290304` `880:1741824` `881:326592` `882:290304` `883:1016064`|
|2|2|`886:72576` `887:217728` `888:24192` `889:290304` `890:290304` `891:580608` `892:870912` `894:145152` `895:1306368` `896:290304` `897:870912` `898:2032128`|
|2|2|`900:48384` `902:96768` `903:96768` `904:193536` `905:1161216` `906:870912` `907:1451520` `908:580608` `909:1161216` `910:3483648` `911:4064256` `912:3483648`|
|2|2|`913:580608` `914:1451520` `915:2032128` `916:1451520` `917:3483648` `919:362880` `920:362880` `921:2903040` `922:1741824` `923:3483648` `924:2612736` `925:580608`|
|2|2|`926:1161216` `927:5184` `928:41472` `930:20736` `931:248832` `932:25920` `933:82944` `934:124416` `935:331776` `936:51840` `937:497664` `938:414720`|
|2|2|`939:290304` `940:248832` `941:10368`|
|2|3|`14:2016` `18:8064` `19:20160` `24:72576` `25:32256` `26:8064` `29:32256` `31:96768` `32:32256` `36:51840` `37:20736` `39:5184`|
|2|3|`104:12096` `109:15120` `199:145152` `206:362880` `207:90720` `210:16128` `214:241920` `217:241920` `218:60480` `317:725760` `319:725760` `337:217728`|
|2|3|`349:2177280` `352:1088640` `358:1088640` `363:362880` `365:725760` `367:2903040` `368:1451520` `369:483840` `372:483840` `373:40320` `437:483840` `498:2177280`|
|2|3|`511:2177280` `516:2903040` `518:5806080` `522:544320` `525:2903040` `526:1451520` `529:483840` `530:1935360` `542:4354560` `543:2177280` `546:4354560` `548:1451520`|
|2|3|`549:2903040` `552:4354560` `555:5806080` `557:725760` `558:483840` `559:483840` `642:1088640` `647:2177280` `651:544320` `676:544320` `685:4354560` `688:4354560`|
|2|3|`697:2177280` `701:4354560` `713:6531840` `715:8709120` `717:2177280` `719:6531840` `722:4354560` `723:3265920` `724:4354560` `726:4354560` `728:1088640` `731:6531840`|
|2|3|`732:3265920` `734:1088640` `795:870912` `803:2612736` `812:435456` `820:1306368` `821:1306368` `827:1306368` `838:1306368` `839:435456` `841:2612736` `844:1741824`|
|2|3|`849:5225472` `850:1741824` `854:10450944` `857:1306368` `858:2612736` `860:3919104` `861:435456` `862:217728` `883:217728` `898:435456` `911:870912` `915:435456`|
|2|3|`917:1741824` `920:217728` `923:1741824` `926:870912` `935:62208` `939:124416`|
|3|3|`4:36` `5:64` `20:3150` `26:30240` `32:60480` `33:6720` `37:25920` `110:630` `219:60480` `370:725760` `373:80640` `374:241920`|
|3|3|`531:645120` `553:725760` `556:2903040` `559:967680` `560:26880` `725:4354560` `733:1632960` `734:2177280` `855:2612736` `861:870912` `925:290304`|

### Block 16301 — cardinality-cap

Order 1, 2 nonzero entries. Comparison result **PASS**. The coefficient SHA-256 is `e8b85f9ac3c67c39d83fa1885f3610a1387e1586af152a88284e0fc39f1fb328`.

Actual label

```json
{"family":"cardinality-cap","bound":186}
```

| Row $i$ | Column $j$ | Variable indices and coefficients |
|---:|---:|---|
|0|0|`0:186` `1:-1`|

### Block 16307 — cardinality-scalar

Order 1, 4 nonzero entries. Comparison result **PASS**. The coefficient SHA-256 is `c472b74ed9ab5ec9263fc8c225f3c9b2dad99f0697e0d7717c7fa298c1c834a7`.

Actual label

```json
{"family":"cardinality-scalar","bound":186,"anchor":[0,6,6,0,6,0,0,0]}
```

| Row $i$ | Column $j$ | Variable indices and coefficients |
|---:|---:|---|
|0|0|`6:183` `41:-1200` `42:-10800` `43:-3375`|

### Block 16363 — cardinality-rooted

Order 1, 120 nonzero entries. Comparison result **PASS**. The coefficient SHA-256 is `34b273ab4d6d8d6d3168a3ccb7c54501b6dc823c145791956a8fbe55c8676469`.

Actual label

```json
{"family":"cardinality-rooted","bound":186,"source":{"family":"D1","anchors":[63],"k":3,"row_arity":1,"cells":[[6,7,8,9,10,11,12,13,14,15,16,17],[0,1,2,3,4,5]],"shapes":[[6,6],[6]],"columns":[[[[0,0,0,0,0,0],[1,1,1,1,1,1]],[[0,0,0,0,0,0]]]],"empty_column":false,"row_cardinalities":[2],"positive_scale_removed":64},"kept_rows":[0],"positive_scale_removed":1}
```

| Row $i$ | Column $j$ | Variable indices and coefficients |
|---:|---:|---|
|0|0|`2:184` `6:182` `7:-1170` `8:2295` `9:-540` `11:-4460` `12:-3600` `15:-3600` `16:-3375` `21:-3600` `27:-400` `41:-1160`|
|0|0|`42:-10800` `43:-3375` `44:1080` `45:3600` `46:1080` `47:960` `48:3600` `49:3120` `50:3600` `51:9000` `52:1800` `53:3000`|
|0|0|`54:1200` `55:2250` `56:-540` `57:-900` `58:-4320` `59:-2160` `60:-2700` `63:1080` `66:120` `67:-2880` `68:-3600` `69:2160`|
|0|0|`72:-600` `73:-720` `74:-3600` `84:-600` `85:-225` `88:20` `89:360` `90:900` `91:400` `111:960` `112:7200` `113:3120`|
|0|0|`114:9000` `115:3600` `116:9000` `117:12000` `118:2400` `119:9000` `120:3000` `121:-900` `122:-5760` `123:-3600` `125:1080` `126:-6480`|
|0|0|`127:6480` `128:-17280` `129:-21600` `130:6480` `133:-1800` `134:-10800` `142:-2880` `143:-14400` `144:-9600` `149:-9600` `150:-3600` `157:-600`|
|0|0|`158:1620` `159:9720` `160:8100` `163:6480` `164:16200` `165:7200` `185:360` `186:2700` `187:3600` `220:120` `221:-4320` `222:2160`|
|0|0|`223:-7200` `225:-600` `226:-6480` `227:-8100` `230:-21600` `231:-14400` `235:-7200` `242:-3600` `243:-9600` `246:-2400` `250:1620` `251:3240`|
|0|0|`252:19440` `253:16200` `256:16200` `257:7200` `267:9720` `268:24300` `269:10800` `276:16200` `277:21600` `308:900` `309:3600` `375:20`|
|0|0|`376:720` `377:1800` `379:800` `381:3240` `382:2700` `385:16200` `386:7200` `392:7200` `407:8100` `408:10800` `415:7200` `434:400`|
