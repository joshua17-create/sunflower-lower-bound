# A proof by hand that μ₃ ≥ (20 + 2√66)^{1/8}

Let f(n) be the largest size of a family of subsets of an n-element set that contains no three distinct members A, B, C with A∩B = A∩C = B∩C. Such a triple is a *3-sunflower*. The Erdős–Szemerédi sunflower-free capacity is μ₃ = lim f(n)^{1/n}. The limit exists by a standard supermultiplicativity argument; see for example the page for constant 49 in Tao's optimization-constants repository.

This note proves μ₃ ≥ (20 + 2√66)^{1/8} = 1.56642866… without a computer. Every step can be checked with pencil and paper. `proof_checks.py` confirms each individual claim by brute force, but the argument does not rely on it.

The stronger bound 1.5666969… in README.md also uses the 121-set family on Z₁₁. That family is checked only by computer (`verify.py`).

## Setting

Let V = {1,…,5} and let E be the 10 edges of the complete graph K₅ on V. These edges are the 10 points. A "graph" below means a set of edges of K₅, always taken on the full vertex set V, so some vertices may be isolated. Two edges are *disjoint*, written e ⟂ f, if they share no vertex.

A **chair** is a tree on V with degree sequence 3,2,1,1,1: a centre with three legs, one of which is extended by one more edge. A **triangle-plus-edge** is a triangle together with the edge joining the two remaining vertices. Let 𝒳 be the set of all 60 chairs and all 10 triangle-plus-edges.

The family F consists of 75 members:

* the 5 **stars** S_v = {all edges at v}, each with 4 edges;
* the 70 six-edge graphs X whose complement E ∖ X lies in 𝒳.

The complement of a triangle-plus-edge is a K₂,₃. F is invariant under all permutations of V.

### Facts about K₅

* **(D1)** Every edge is disjoint from exactly three edges: those on the other three vertices.
* **(D2)** No three edges are pairwise disjoint, since that would need six vertices.
  * There are also no four distinct edges with e₁ ⟂ e₂ ⟂ e₃ ⟂ e₄ ⟂ e₁.
  * Proof: e₁ and e₂ cover four vertices; call the fifth one v. Since e₃ ⟂ e₂, the edge e₃ lies inside e₁ ∪ {v}, and as e₃ ≠ e₁ it contains v. Likewise e₄ ⟂ e₁ and e₄ ≠ e₂ force e₄ to contain v. Then e₃ and e₄ meet, a contradiction.
* **(D3)** Two edges that meet, uv and uw, have exactly one common disjoint edge: the edge on the remaining two vertices.
* **(E)** Up to relabelling, the 4-edge graphs on V are of the six types in the table below.
  * Their numbers are 5, 5·3 = 15, 5·12 = 60, 5!/2 = 60, 5·4·3 = 60 and 10. These add up to 210 = C(10,4), so the list is complete.
  * The six types are pairwise non-isomorphic. P₅ and triangle-plus-edge share the degree sequence 2,2,2,1,1, but only the second contains a triangle.
  * The table also lists each type's disjoint pairs of edges; a, b, c, … denote distinct vertices.

| type | count | disjoint pairs of edges |
|---|---|---|
| star K₁,₄ | 5 | none |
| C₄ + K₁ (cycle abcd) | 15 | {ab, cd}, {bc, da} |
| paw + K₁ (triangle abc, edge cd) | 60 | {ab, cd} |
| path P₅ (abcde) | 60 | {ab, cd}, {ab, de}, {bc, de} |
| chair (centre c, legs ca–ab, cl, cl′) | 60 | {ab, cl}, {ab, cl′} |
| triangle-plus-edge (triangle abc, edge de) | 10 | {de, ab}, {de, bc}, {de, ca} |

* In particular, every member of 𝒳 has no isolated vertex and maximum degree at most 3.

**Lemma 1.** A 4-edge graph lies in 𝒳 exactly when it has at least two disjoint pairs of edges and they all contain one common edge.

*Proof.* Read off the table. ∎

In terms of the *disjointness graph*, whose vertices are edges of K₅ joined when they are disjoint, Lemma 1 says this: a 4-edge graph lies in 𝒳 exactly when the disjointness graph on its four edges is a star with two or three edges, possibly plus an isolated vertex.

## F is Sperner

Two 6-edge members cannot be nested, and neither can two stars. Suppose a star S_v were contained in a 6-edge member X. Then v would be isolated in E ∖ X, which is impossible for a member of 𝒳.

## F has no 3-sunflower

### Three 6-edge members

Let X, Y, Z be distinct 6-edge members with complements x, y, z ∈ 𝒳. Since X ∩ Y = E ∖ (x ∪ y), the condition X∩Y = X∩Z = Y∩Z says x ∪ y = x ∪ z = y ∪ z. Equivalently, **no edge lies in exactly one of x, y, z**.

Let t edges lie in all three sets and d edges in exactly two. Then 3t + 2d = 4 + 4 + 4 = 12, so t is even. Also t < 4, because the sets are distinct. This leaves two cases.

* **(i) t = 2, d = 3.** Then U = x ∪ y ∪ z has 5 edges. Each of x, y, z consists of the 2 common edges plus 2 of the 3 edges lying in exactly two sets. So x, y, z are U minus three different edges of U.
* **(ii) t = 0, d = 6.** Put A = z ∩ x, B = x ∩ y and C = y ∩ z. These are pairwise disjoint 2-edge sets, and x = A∪B, y = B∪C, z = C∪A.

**Lemma 2.** No 5-edge graph U has three edges e with U ∖ e ∈ 𝒳. This rules out case (i).

*Proof.* Let D_U be the disjointness graph on the five edges of U. The disjointness graph of U ∖ e is D_U − e. By Lemma 1, U ∖ e ∈ 𝒳 exactly when D_U − e is a star with at least two edges, possibly plus an isolated vertex.

By (D2), D_U has no 3-cycle or 4-cycle. So any cycle in D_U is a 5-cycle through all five vertices, and it has no chord, since a chord would create a shorter cycle. Hence D_U is either a forest or a 5-cycle. By (D1), D_U has maximum degree at most 3.

* **D_U is a 5-cycle.** Every D_U − e is a path with three edges, so no e works.
* **D_U has a vertex f of degree 3.** Let f be adjacent to g₁, g₂, g₃.
  * By (D1) these are all the edges of K₅ disjoint from f, so the fifth edge h of U meets f.
  * By (D2), no gᵢ is disjoint from another gⱼ, because f, gᵢ, gⱼ would then be pairwise disjoint.
  * By (D3), f and h have exactly one common disjoint edge. It is disjoint from f, so it is some gᵢ, say g₁.
  * An adjacency between h and gⱼ would make gⱼ a common disjoint edge of f and h. So h is adjacent to g₁ only, and the edges of D_U are fg₁, fg₂, fg₃, g₁h.
  * Then D_U − h and D_U − g₁ are stars. D_U − g₂ and D_U − g₃ are paths with three edges, and D_U − f has one edge. Exactly two choices of e work.
* **Otherwise D_U is a disjoint union of paths.** Then D_U − e must be a 3-vertex path plus an isolated vertex.
  * If D_U contains no path on at least three vertices, nothing works.
  * Otherwise D_U is P₅, P₄+K₁, P₃+P₂ or P₃+2K₁. The choices of e that work are, respectively: the 2nd and 4th vertex of P₅; the two ends of P₄; the two ends of P₂; the two isolated vertices.

In every case at most two choices of e work. ∎

**Lemma 3.** There are no pairwise disjoint 2-edge sets A, B, C with A∪B, B∪C, C∪A ∈ 𝒳. This rules out case (ii).

*Proof.* Let W = A∪B∪C, which has 6 edges, and let L = E ∖ W, which has 4. Each edge of W lies in exactly two of A∪B, B∪C, C∪A. So for each vertex v,

deg_W(v) = ½ (deg_{A∪B}(v) + deg_{B∪C}(v) + deg_{C∪A}(v)) ≥ 3/2,

because members of 𝒳 have no isolated vertex. Hence deg_W(v) ≥ 2 and deg_L(v) ≤ 2. By (E), L is C₄+K₁, P₅ or a triangle-plus-edge.

Since W ∖ A = B∪C and so on, we would need to split W into three pairs P, each with W ∖ P ∈ 𝒳. Call such a pair *admissible*. We show the split is impossible in each case.

* **L is a triangle-plus-edge** (triangle abc, edge de).
  * Then W = K₂,₃ with parts {d,e} and {a,b,c}. W contains no triangle, so each W ∖ P is a chair.
  * The centre of such a chair has three neighbours, so it is d or e (a, b, c have only two neighbours in W). Say it is d. Then the chair is {da, db, dc, ex}, and P = W ∖ chair consists of two edges at e.
  * So every admissible pair consists of two edges at d or two edges at e. The three edges at d cannot be split into such pairs.
* **L is a path P₅** (v₁v₂v₃v₄v₅).
  * Then W ∖ {v₁v₅} = {v₁v₃, v₃v₅, v₅v₂, v₂v₄, v₄v₁} is a 5-cycle.
  * For the pair P containing v₁v₅, W ∖ P is a 5-cycle minus one edge, which is a path P₅. A P₅ is not in 𝒳.
* **L is C₄+K₁** (cycle v₁v₂v₃v₄, isolated vertex u).
  * W consists of the diagonals v₁v₃ and v₂v₄ and the four edges uvᵢ. Each vᵢ has degree 2 in W.
  * W ∖ P has no isolated vertex, so an admissible P never contains both W-edges at some vᵢ, namely its diagonal and uvᵢ.
  * u has degree 4 in W but at most 3 in W ∖ P, so P contains an edge at u.
  * If P = {uvᵢ, uvⱼ}, then W ∖ P is a triangle-plus-edge when vᵢ, vⱼ are opposite corners, and a path P₅ when they are adjacent.
  * So **every admissible pair has one of two forms**: two u-edges to opposite corners, or a u-edge together with the diagonal not through its corner.
  * Two diagonals never form an admissible pair. So in a split into three admissible pairs, each diagonal is paired with a u-edge to a corner of the *other* diagonal.
  * The remaining two u-edges then go to one corner of each diagonal. These are adjacent corners, so the last pair is not admissible. ∎

### Triples containing stars

* **Three stars.** S_u ∩ S_v = {uv}, which differs from S_u ∩ S_w = {uw}.
* **Two stars S_u, S_v and a 6-edge member X.** All three intersections would equal {uv}. So X contains uv but no other edge at u or v. Then X lies inside {uv} together with the triangle on the other three vertices, which forces |X| ≤ 4. Impossible.
* **One star S_u and two 6-edge members X ≠ Y.** Suppose S_u∩X = S_u∩Y = X∩Y = K.
  * X ∖ K and Y ∖ K are disjoint and contain no edge at u. So they lie among the six edges of the K₄ on V ∖ {u}, which gives 12 − 2|K| ≤ 6, i.e. |K| ≥ 3.
  * If |K| = 4, then S_u ⊆ X and u is isolated in E ∖ X, which is impossible. So |K| = 3, and Q′ = X ∖ K and Q = Y ∖ K split that K₄ into two 3-edge sets.
  * Let e = uw be the edge at u outside K. Then E ∖ X = {e} ∪ Q and E ∖ Y = {e} ∪ Q′. Both lie in 𝒳 and have u as a leaf attached at w.
  * There are three ways this can happen. In a triangle-plus-edge, Q is the triangle missing w. In a chair, either u ends the extended leg, and Q is a K₁,₃ with w as a leaf; or u is a short leg at the centre w, and Q is a 3-edge path with w as an inner vertex.
  * So both Q and Q′ must be one of: (1) the triangle missing w, (2) a K₁,₃ with w as a leaf, (3) a 3-edge path with w as an inner vertex.
  * But Q′ is the complement of Q in K₄. The complement of (1) is the star centred at w. The complement of (2) is a triangle through w. The complement of (3) is a 3-edge path with w at an end. None of these is of the form (1)–(3). ∎

**So F is Sperner and contains no 3-sunflower.**

## The port counts

**Single edges.** Every edge e lies in 2 stars. Each member of 𝒳 has 4 edges, and the permutations of V act transitively on edges, so e lies in 70·4/10 = 28 members of 𝒳. Hence 70 − 28 = 42 six-edge members contain e, and e lies in **44** members of F.

**Pairs of edges.** Count pairs of edges inside members of 𝒳.

* *Meeting pairs.* A chair contains 6 − 2 = 4 of them and a triangle-plus-edge 6 − 3 = 3, so 60·4 + 10·3 = 270 in total. K₅ has 30 meeting pairs and the permutations act transitively on them, so each lies in 9 members of 𝒳.
* *Disjoint pairs.* There are 60·2 + 10·3 = 150 in total, spread over the 15 disjoint pairs of K₅, so each lies in 10 members of 𝒳.
* *Combining.* By inclusion–exclusion, the six-edge members containing both e and f number 70 − 28 − 28 + 9 = 23 if e and f meet, and 70 − 28 − 28 + 10 = 24 if they are disjoint. A meeting pair also lies in the one star through its common vertex. So in either case **24** members of F contain both e and f.

Now fix two distinct edges p (the *in-port*) and q (the *out-port*). Let N_ab be the number of D ∈ F with [p ∈ D] = a and [q ∈ D] = b. Then

N₁₁ = 24, N₁₀ = N₀₁ = 44 − 24 = 20, N₀₀ = 75 − 44 − 44 + 24 = 11.

## The bound

The gluing lemma of Deuber–Erdős–Gunderson–Kostochka–Meyer (1997, Lemma 4.4; a proof is in README.md) says the following. Let 𝒜 and ℬ be Sperner, 3-sunflower-free families on disjoint ground sets, and let a and b be points of them. Then

{(C∖a) ∪ D : a ∈ C ∈ 𝒜, b ∉ D ∈ ℬ} ∪ {C ∪ (D∖b) : a ∉ C ∈ 𝒜, b ∈ D ∈ ℬ}

is again Sperner and 3-sunflower-free, on the union of the two ground sets with a and b *deleted*.

Build the chain as follows.

* Let G₁ = F, with out-port q₁ = q.
* Given G_k, let G_{k+1} be the result of gluing G_k at its out-port q_k to a fresh copy of F at that copy's in-port p. The new out-port q_{k+1} is the copy of q, which survives because p ≠ q.

Each G_k is Sperner and 3-sunflower-free, on 10 + 8(k−1) = 8k + 2 points.

Let u_k and w_k count the members of G_k that avoid and contain q_k. The members of G_{k+1} arise in two ways:

* A member containing q_k is combined with a copy-member avoiding p.
* A member avoiding q_k is combined with a copy-member containing p.

Sorting the copy-members by whether they contain q gives

u_{k+1} = u_k·N₁₀ + w_k·N₀₀,  w_{k+1} = u_k·N₁₁ + w_k·N₀₁,

that is, (u_{k+1}, w_{k+1}) = (u_k, w_k)·M with

M = [[N₁₀, N₁₁], [N₀₀, N₀₁]] = [[20, 24], [11, 20]].

The characteristic polynomial of M is (20 − λ)² − 264, so its largest eigenvalue is λ = 20 + √264 = 20 + 2√66. A positive left eigenvector is (√11, √24).

The starting vector is (u₁, w₁) = (31, 44). With c = 44/√24, we have (31, 44) ≥ c·(√11, √24) entrywise. M has nonnegative entries, so this inequality is preserved under multiplication by M, and

|G_k| = u_k + w_k ≥ c(√11 + √24)·λ^{k−1}.

Every n satisfies n ≥ 8k + 2 for k = ⌊(n−2)/8⌋. Also f is nondecreasing, because a family on n points is also a family on n + 1 points. Hence

f(n)^{1/n} ≥ (c(√11+√24)·λ^{k−1})^{1/n} → λ^{1/8}   as n → ∞,

and therefore μ₃ ≥ λ^{1/8} = (20 + 2√66)^{1/8} = 1.5664286622…  ∎
