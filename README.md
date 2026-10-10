# Small base families for the 3-sunflower-free capacity

**Theorem.** Let f(n) be the largest size of a family of subsets of an n-element set containing no three distinct members A, B, C with A∩B = A∩C = B∩C, and let μ₃ = lim f(n)^{1/n} (the Erdős–Szemerédi 3-sunflower-free capacity). Then

  μ₃ ≥ (1042 + 2√261207)^{1/17} = 1.5666969787…

and, using a single base family,

  μ₃ ≥ (20 + 2√66)^{1/8} = 1.5664286622…

Run `python3 verify.py` (Python 3, standard library only, about 15 seconds).

**A proof by hand.** [PROOF.md](PROOF.md) proves the single-family bound with pencil and paper. Everything in it can be checked without a computer: the 75-set family is Sperner and 3-sunflower-free, its port counts are as stated, and the bound follows. `proof_checks.py` confirms each step of that proof by brute force.

## Previous bounds

* μ₃ > 1.551: Deuber, Erdős, Gunderson, Kostochka and Meyer (DEGKM), *Intersection statements for systems of sets*, J. Combin. Theory Ser. A 79 (1997) 118–132. This is the published record.
* μ₃ ≥ 1.554: an unpublished result of Naslund, cited in Naslund–Sawin, *Upper bounds for sunflower-free sets*, Forum Math. Sigma 5 (2017).
* μ₃ ≥ 1.5621220560: IRSI (Aspect5/sunflower-free-capacity, v1.1, 2 Oct 2026; pull request #200 to the optimization-constants repository, not yet merged).
  * They chained DEGKM's gluing lemma and counted the chain exactly with a 2×2 transfer matrix.
  * Their base family has 6,150 sets on 20 points.
  * **The exact transfer-matrix count used below is theirs.** We found it independently a week later.
* The best upper bound is μ₃ ≤ 3/2^{2/3} ≈ 1.8899 (Naslund–Sawin).

What is new here are two much smaller base families, and alternating them in the chain.

## The two base families

**F (75 sets on 10 points).** Take the points to be the 10 edges of the complete graph K₅. The members of F are the edge sets of:

* the 5 stars K₁,₄;
* the 10 complete bipartite graphs K₂,₃;
* the 60 graphs whose complement in K₅ is a *chair*: a vertex of degree 3 with one of its three legs extended by one more edge (degree sequence 3,2,1,1,1).

F is invariant under S₅ and is Sperner. A star at v lies in a 6-edge graph only if v is isolated in that graph's complement, and neither a chair nor a triangle-plus-edge (the complement of K₂,₃) has an isolated vertex. A brute-force check of all 67,525 triples confirms that F has no 3-sunflower.

Every edge lies in 44 members, and every pair of edges, adjacent or not, lies in 24.

**H (121 sets on 11 points).** Take the points to be Z₁₁. The members of H are the images of the three sets {0,1,2,3,4,5,6}, {0,1,2,3,4,6,7} and {0,1,2,6,9} under the 55 maps x ↦ ax + b, where a is a nonzero square mod 11. The first two sets give 55 images each and the third gives 11.

H is Sperner. A brute-force check of all 287,980 triples confirms that H has no 3-sunflower. The link of any point in H, meaning the members containing it with that point removed, is a copy of F.

Both families consist, apart from a few small sets, of the complements of exactly one third of all 4-element subsets.

## Gluing lemma (DEGKM 1997, Lemma 4.4)

**Lemma.** Let A and B be Sperner, 3-sunflower-free families on disjoint ground sets X and Y, and let a ∈ X and b ∈ Y. Then

G = {(C∖a) ∪ D : a ∈ C ∈ A, b ∉ D ∈ B} ∪ {C ∪ (D∖b) : a ∉ C ∈ A, b ∈ D ∈ B}

is Sperner and 3-sunflower-free on (X ∪ Y) ∖ {a, b}, and |G| = |A_a|·|B_{b̄}| + |A_{ā}|·|B_b|. Here A_a denotes the members of A containing a, and A_{ā} those not containing it.

*Proof.* Write each member of G as P ∪ Q with P ⊆ X∖a and Q ⊆ Y∖b.

Set A₁ = {C∖a : a ∈ C ∈ A}, A₀ = {C ∈ A : a ∉ C}, and define B₁ and B₀ the same way. Then G = {P ∪ Q : P ∈ A₁, Q ∈ B₀} ∪ {P ∪ Q : P ∈ A₀, Q ∈ B₁}.

A₁ and A₀ are disjoint, since P and P ∪ {a} cannot both lie in the Sperner family A. So the two parts of G are disjoint, and |G| is as stated.

*G is Sperner.* Suppose P∪Q ⊆ P′∪Q′.
* If both sets lie in the same part, Sperner-ness of A and of B gives P = P′ and Q = Q′, so the two members are equal.
* If P∪Q ∈ A₁⊗B₀ and P′∪Q′ ∈ A₀⊗B₁, then Q ⊆ Q′ ⊊ Q′∪{b} with both Q and Q′∪{b} in B, which is impossible.
* The remaining case is symmetric, using A.

*G has no 3-sunflower.* Let three distinct members Pᵢ∪Qᵢ of G have equal pairwise intersections. Then the triples (Pᵢ) and (Qᵢ) each have equal pairwise intersections.
* **All three in A₁⊗B₀.** The Qᵢ lie in B. If two coincide, say Q₁ = Q₂, then Q₁ ⊆ Q₃, so Q₁ = Q₃ by Sperner-ness. Otherwise the Qᵢ form a sunflower in B. Hence Q₁ = Q₂ = Q₃, the Pᵢ are distinct, and the sets Pᵢ∪{a} form a sunflower in A. The case of all three in A₀⊗B₁ is symmetric.
* **Two in A₁⊗B₀, one in A₀⊗B₁.** The sets Q₁, Q₂ and Q₃∪{b} lie in B, and since b ∉ Q₁, Q₂ they again have equal pairwise intersections. If Q₁ ≠ Q₂ they form a sunflower in B. If Q₁ = Q₂, then Q₁ ⊆ Q₃ ⊊ Q₃∪{b}, which contradicts Sperner-ness. The case of two in A₀⊗B₁ is symmetric, using A. ∎

## The chain and its transfer matrix

Give each base family an in-port x and an out-port z, which are two distinct points. Start with R₁ = a copy of a base family. At each step glue a fresh copy, so that R_{k+1} = G(R_k, z_k, copy, x_copy), and let z_{k+1} be the copy's out-port. Every R_k is Sperner and 3-sunflower-free by the lemma.

Let N_{ab} = #{D : [x ∈ D] = a, [z ∈ D] = b}, and let u_k and w_k count the members of R_k that avoid and contain z_k. Then

(u_{k+1}, w_{k+1}) = (u_k, w_k) · T,  T = [[N₁₀, N₁₁], [N₀₀, N₀₁]].

The reason is that sets containing z_k are combined with copy members avoiding x, and vice versa.

* For F: N = [[11,20],[20,24]] for every choice of ports, so T_F = [[20,24],[11,20]].
* For H: N = [[15,31],[31,44]], so T_H = [[31,44],[15,31]].

**One family.** Chaining copies of F adds 8 points per copy. Hence f(8k+2) ≥ |R_k| ≥ c·λ^{k−1} for a constant c > 0, with λ = 20 + 2√66 the Perron root of T_F (see PROOF.md). Since f is nondecreasing, μ₃ ≥ λ^{1/8} = 1.56642866…

**Alternating F and H.** This adds 8 + 9 = 17 points per period. The product matrix is

T_F T_H = [[980, 1624], [641, 1104]],

whose Perron root is 1042 + √1044828 = 1042 + 2√261207 ≈ 2064.168. Hence μ₃ ≥ (1042 + 2√261207)^{1/17} = 1.56669697…

The alternation gains because the two matrices have different Perron eigenvectors. No other periodic sequence of length at most 4, drawn from these and several other families (including IRSI's), does better.

`verify.py` also builds the first glued families explicitly and compares their sizes with the matrix recursion:

| Glued families | Sets | Points |
|---|---|---|
| F+F | 2,728 | 18 |
| F+H | 4,349 | 19 |
| F+F+F | 98,920 | 26 |
| F+H+F | 155,892 | 27 |

It checks the first two completely for Sperner-ness and 3-sunflowers.

## How these were found

Both families came out of computer searches. The searches used local search over all subsets, and exact CP-SAT optimisation over unions of orbits of a permutation group: S₅ acting on the edges of K₅, and the affine maps of Z₁₁. The families were then recognised by their symmetry.

* **Optimality.** Local search never found a better family on 10 or 11 points (75 and 121), and CP-SAT proves these optimal among S₅-invariant and Z₁₁-invariant families respectively.
* **Larger ground sets.** Extending H by one point gives a 187-set family on 12 points, with H as the link of every point. Its chain rate is lower (1.5588), however.
* **Other searches.** Searches on 12 to 20 points have not found a base family with a better chain rate so far. These covered many symmetry groups, uniform layers and local search.

## AI disclosure

The searches, checks, proofs, write-up and this repository were produced by Claude (Anthropic, Claude Opus 5.5).

A second Claude agent, given only the statement, re-derived everything without access to the searching code:
* the families from their verbal descriptions;
* a proof of the lemma, together with an exhaustive test of the lemma on all families on at most 5 points;
* the explicit 18- and 26-point families, set by set;
* the eigenvalues.

That agent also located the IRSI prior work.

A third Claude agent refereed the hand proof in PROOF.md line by line, with its own code, and accepted it; its wording suggestions are incorporated.

## Status

* Not peer-reviewed. **No human has checked the proof yet.** The checks so far are the two scripts here, which anyone can run, and the independent AI referees described above.
* This repository is published by a non-specialist on Claude's behalf. The publisher has not verified the mathematics personally.
* Corrections and questions are welcome as GitHub issues.
