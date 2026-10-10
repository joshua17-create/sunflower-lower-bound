#!/usr/bin/env python3
"""
Verification of two lower bounds for the Erdos-Szemeredi 3-sunflower-free capacity

    mu_3 = lim_{n -> oo} f(n)^(1/n),

f(n) = largest size of a family of subsets of an n-set with no three distinct members A, B, C
satisfying A&B = A&C = B&C.

    (1)  mu_3 >= (20 + 2*sqrt(66))^(1/8)        = 1.566428662...   (one 75-set family on 10 points)
    (2)  mu_3 >= (1042 + 2*sqrt(261207))^(1/17)  = 1.566696978...   (alternating it with a 121-set family on 11 points)

Pure Python 3, standard library only.  Runs in well under a minute.

What is checked here, and what is proved in README.md:
  * checked: both base families are Sperner (no member contains another) and 3-sunflower-free, by brute force
    over all triples; their 2x2 port matrices; the eigenvalues and the numerical bounds (with exact integer
    arithmetic); and, as a sanity check of the gluing lemma and of the counting, explicit glued families
    (18, 19 and 26/27 points) whose sizes must match the matrix recursion, the two smaller ones fully checked.
  * proved in README.md: the gluing lemma (Deuber-Erdos-Gunderson-Kostochka-Meyer 1997, Lemma 4.4) and the
    transfer-matrix count of the glued chain, which give the bounds.
"""
import itertools, sys, time
from decimal import Decimal, getcontext

# ---------------------------------------------------------------- the two base families
def k5_family():
    """75 edge sets of K5 (points = the 10 edges): 5 stars K_{1,4}, 10 K_{2,3}'s, 60 complements of chairs."""
    V = range(5)
    E = list(itertools.combinations(V, 2))
    idx = {e: i for i, e in enumerate(E)}
    full = (1 << 10) - 1
    mask = lambda es: sum(1 << idx[tuple(sorted(e))] for e in es)
    stars = {mask([(v, w) for w in V if w != v]) for v in V}
    k23 = {mask([(u, w) for u in P for w in V if w not in P]) for P in itertools.combinations(V, 2)}
    chairs = set()
    for c in V:                                   # centre of the chair (degree 3)
        for a in V:
            for b in V:
                if len({c, a, b}) < 3: continue   # leg c-a-b; the other two vertices are leaves of c
                leaves = [v for v in V if v not in (c, a, b)]
                chairs.add(mask([(c, a), (a, b)] + [(c, l) for l in leaves]))
    assert (len(stars), len(k23), len(chairs)) == (5, 10, 60)
    fam = sorted(stars | k23 | {full ^ t for t in chairs})
    assert len(fam) == 75
    return fam, 10

def z11_family():
    """121 subsets of Z_11: the images of {0..6}, {0,1,2,3,4,6,7}, {0,1,2,6,9} under x -> a*x + b, a a nonzero square."""
    p, QR = 11, [1, 3, 4, 5, 9]
    reps = [[0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 6, 7], [0, 1, 2, 6, 9]]
    fam = set()
    for r in reps:
        for a in QR:
            for b in range(p):
                fam.add(sum(1 << ((a * x + b) % p) for x in r))
    fam = sorted(fam)
    assert len(fam) == 121
    return fam, 11

# ---------------------------------------------------------------- checks
def brute_force(fam):
    for A in fam:
        for B in fam:
            if A != B and A & B == A:
                return "NOT SPERNER"
    for A, B, C in itertools.combinations(fam, 3):
        if A & B == A & C == B & C:
            return "SUNFLOWER %s %s %s" % (A, B, C)
    return "OK"

def fast_check(fam):
    """Sperner + 3-sunflower-free check for larger families.
    For each S, group the other members by their intersection K with S; a sunflower through S exists
    iff some group contains two members that are disjoint outside K."""
    S_all = sorted(set(fam))
    if len(S_all) != len(fam): return "DUPLICATES"
    for S in S_all:
        groups = {}
        for T in S_all:
            if T == S: continue
            K = S & T
            if K == S or K == T: return "NOT SPERNER"
            groups.setdefault(K, []).append(T & ~K)
        for K, rest in groups.items():
            for i in range(len(rest)):
                ri = rest[i]
                for j in range(i + 1, len(rest)):
                    if ri & rest[j] == 0:
                        return "SUNFLOWER"
    return "OK"

def port_matrix(fam, x, z):
    N = [[0, 0], [0, 0]]
    for D in fam:
        N[(D >> x) & 1][(D >> z) & 1] += 1
    # T[i][k]: previous out-bit i forces in-bit 1-i; new out-bit k
    return [[N[1][0], N[1][1]], [N[0][0], N[0][1]]], N

def matmul(P, Q):
    return [[P[0][0]*Q[0][0] + P[0][1]*Q[1][0], P[0][0]*Q[0][1] + P[0][1]*Q[1][1]],
            [P[1][0]*Q[0][0] + P[1][1]*Q[1][0], P[1][0]*Q[0][1] + P[1][1]*Q[1][1]]]

def perron(T):
    """Perron root of a positive 2x2 integer matrix as (h, D): lambda = h + sqrt(D) with h = tr/2, D = (tr/2)^2 - det
    (tr is even for the matrices used here, so h and D are integers)."""
    tr = T[0][0] + T[1][1]; det = T[0][0]*T[1][1] - T[0][1]*T[1][0]
    assert tr % 2 == 0
    return tr // 2, (tr // 2) ** 2 - det

def lam_gt_power(h, D, num, den, e):
    """Exact test of  h + sqrt(D) > (num/den)^e."""
    r = num ** e; s = den ** e          # compare  h + sqrt(D)  with  r/s
    lhs = r - h * s                     # sqrt(D) > (r - h s)/s
    if lhs < 0: return True
    return D * s * s > lhs * lhs

# ---------------------------------------------------------------- chains (explicit, for sanity checks)
def glue(cur, npos, outport, fam, m, x, z):
    """G(cur, outport, fam-copy, x): DEGKM Lemma 4.4 with a fresh copy of fam whose point x is deleted."""
    mp, p = {}, npos
    for v in range(m):
        if v != x: mp[v] = p; p += 1
    D1, D0 = [], []
    for s in fam:
        t = 0
        for v in range(m):
            if v != x and (s >> v) & 1: t |= 1 << mp[v]
        (D1 if (s >> x) & 1 else D0).append(t)
    zb, new = 1 << outport, []
    for C in cur:
        if C & zb: new.extend((C ^ zb) | D for D in D0)
        else: new.extend(C | D for D in D1)
    return new, p, mp[z]

def chain(word, fams):
    f, m, x, z = fams[word[0]]
    cur, npos, outp = list(f), m, z
    for c in word[1:]:
        f, m, x, z = fams[c]
        cur, npos, outp = glue(cur, npos, outp, f, m, x, z)
    used = 0
    for s in cur: used |= s
    return cur, bin(used).count('1')

def predicted(word, fams, mats):
    f, m, x, z = fams[word[0]]
    v = [sum(1 for D in f if not (D >> z) & 1), sum(1 for D in f if (D >> z) & 1)]
    for c in word[1:]:
        T = mats[c]
        v = [v[0]*T[0][0] + v[1]*T[1][0], v[0]*T[0][1] + v[1]*T[1][1]]
    return v[0] + v[1]

# ---------------------------------------------------------------- main
def main():
    t0 = time.time()
    getcontext().prec = 60
    A, mA = k5_family(); B, mB = z11_family()
    print("K5 family : %d sets on %d points, sizes %s" % (len(A), mA, sorted({bin(s).count('1') for s in A})))
    print("            brute force over %d triples: %s" % (len(A)*(len(A)-1)*(len(A)-2)//6, brute_force(A)))
    print("Z11 family: %d sets on %d points, sizes %s" % (len(B), mB, sorted({bin(s).count('1') for s in B})))
    print("            brute force over %d triples: %s" % (len(B)*(len(B)-1)*(len(B)-2)//6, brute_force(B)))

    TA, NA = port_matrix(A, 0, 1)        # K5: ports = edges {0,1} and {0,2}; every choice gives the same counts
    TB, NB = port_matrix(B, 0, 1)        # Z11: ports 0 and 1; every choice gives the same counts
    allA = {tuple(map(tuple, port_matrix(A, x, z)[1])) for x in range(mA) for z in range(mA) if x != z}
    allB = {tuple(map(tuple, port_matrix(B, x, z)[1])) for x in range(mB) for z in range(mB) if x != z}
    print("port counts N[a][b] = #{D : [x in D]=a, [z in D]=b}:  K5 %s (all %d port pairs alike: %s),  Z11 %s (all alike: %s)"
          % (NA, mA*(mA-1), len(allA) == 1, NB, len(allB) == 1))
    print("transfer matrices: T_K5 = %s, T_Z11 = %s" % (TA, TB))

    # bound (1)
    h, D = perron(TA)
    lam = Decimal(h) + Decimal(D).sqrt()
    print("\n(1) lambda(T_K5) = %d + sqrt(%d) = %s" % (h, D, lam))
    print("    mu_3 >= lambda^(1/8) = %s" % (lam ** (Decimal(1) / 8)))
    print("    exact check  lambda > 1.566428^8 : %s" % lam_gt_power(h, D, 1566428, 10**6, 8))
    # bound (2)
    P = matmul(TA, TB)
    h2, D2 = perron(P)
    lam2 = Decimal(h2) + Decimal(D2).sqrt()
    print("\n(2) T_K5 * T_Z11 = %s, lambda = %d + sqrt(%d) = %s" % (P, h2, D2, lam2))
    print("    mu_3 >= lambda^(1/17) = %s" % (lam2 ** (Decimal(1) / 17)))
    print("    exact check  lambda > 1.566696^17 : %s" % lam_gt_power(h2, D2, 1566696, 10**6, 17))

    # explicit chains: sizes must follow the matrix recursion; the smaller ones are checked completely
    fams = {'A': (A, mA, 0, 1), 'B': (B, mB, 0, 1)}
    mats = {'A': TA, 'B': TB}
    print("\nexplicit glued families (Lemma 4.4 applied repeatedly):")
    for word, full_check in (("AA", True), ("AB", True), ("BA", False), ("AAA", False), ("ABA", False)):
        fam, npts = chain(word, fams)
        res = fast_check(fam) if full_check else "size only"
        print("  %-4s %8d sets on %2d points, predicted %8d   %s" % (word, len(fam), npts, predicted(word, fams, mats), res))
    print("\ntotal time %.1fs" % (time.time() - t0))

if __name__ == "__main__":
    main()
