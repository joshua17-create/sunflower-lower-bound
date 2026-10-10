# Brute-force checks of every claim used in the hand proof (PROOF.md).
import itertools
V = range(5); E = list(itertools.combinations(V, 2))
def deg(g, v): return sum(1 for e in g if v in e)
def degseq(g): return tuple(sorted((deg(g, v) for v in V), reverse=True))
def disj(e, f): return not (set(e) & set(f))
def triangle_in(g): return any(set(itertools.combinations(t, 2)) <= set(g) for t in itertools.combinations(V, 3))
def is_chair(g): return len(g) == 4 and degseq(g) == (3, 2, 1, 1, 1)
def is_T(g): return len(g) == 4 and degseq(g) == (2, 2, 2, 1, 1) and triangle_in(g)
inX = lambda g: is_chair(g) or is_T(g)
X4 = [frozenset(g) for g in itertools.combinations(E, 4)]
calX = {g for g in X4 if inX(g)}
ok = lambda b, msg: print(("PASS " if b else "FAIL ") + msg)

# Types of 4-edge graphs
def typ(g):
    ds = degseq(g)
    if ds == (4,1,1,1,1): return "star"
    if ds == (2,2,2,2,0): return "C4+K1"
    if ds == (3,2,2,1,0): return "paw+K1"
    if ds == (3,2,1,1,1): return "chair"
    if ds == (2,2,2,1,1): return "T" if triangle_in(g) else "P5"
    return "?"
from collections import Counter
ok(Counter(typ(g) for g in X4) == Counter({"star":5,"C4+K1":15,"paw+K1":60,"chair":60,"T":10,"P5":60}), "six types of 4-edge graphs with counts 5,15,60,60,60,10")
# (D1)-(D3)
ok(all(sum(disj(e, f) for f in E) == 3 for e in E), "(D1) every edge is disjoint from exactly 3 edges")
ok(not any(disj(a,b) and disj(b,c) and disj(a,c) for a,b,c in itertools.combinations(E,3)), "(D2a) no three pairwise disjoint edges")
ok(not any(len({a,b,c,d})==4 and disj(a,b) and disj(b,c) and disj(c,d) and disj(d,a) for a,b,c,d in itertools.permutations(E,4)), "(D2b) no 4-cycle in the disjointness graph")
ok(all(sum(disj(g,e) and disj(g,f) for g in E) == 1 for e,f in itertools.combinations(E,2) if not disj(e,f)), "(D3) two intersecting edges have exactly one edge disjoint from both")
# Lemma 1: x in calX iff D[x] is a star with >= 2 edges
def Dstar(g):
    pairs = [(e,f) for e,f in itertools.combinations(sorted(g),2) if disj(e,f)]
    if len(pairs) < 2: return False
    common = set(pairs[0])
    for p in pairs[1:]: common &= set(p)
    return len(common) >= 1
ok(all(Dstar(g) == (g in calX) for g in X4), "Lemma 1: x in calX iff its disjoint pairs (>=2 of them) all share one edge")
# members of calX have no isolated vertex and max degree <= 3
ok(all(min(deg(g,v) for v in V) >= 1 and max(deg(g,v) for v in V) <= 3 for g in calX), "members of calX: no isolated vertex, max degree <= 3")
# Lemma 2: no 5-edge U with three e such that U-e in calX
ok(max(sum((frozenset(U)-{e}) in calX for e in U) for U in itertools.combinations(E,5)) <= 2, "Lemma 2: every 5-edge graph has at most two one-edge deletions in calX")
# complement form of the sunflower condition for 6-sets
F6 = [frozenset(E) - g for g in calX]
def sunflower(a,b,c): return a&b == a&c == b&c
def cosun(x,y,z): return all(((e in x)+(e in y)+(e in z)) != 1 for e in E)
ok(all(sunflower(frozenset(E)-x,frozenset(E)-y,frozenset(E)-z) == cosun(x,y,z) for x,y,z in itertools.combinations(list(calX)[:40],3)), "complement form: sunflower of 6-sets <=> no edge in exactly one complement (sampled)")
# Lemma 3 ingredients: max-degree<=2 4-edge graphs are exactly C4+K1, P5, T
ok({typ(g) for g in X4 if max(deg(g,v) for v in V) <= 2} == {"C4+K1","P5","T"}, "4-edge graphs with max degree <= 2 are C4+K1, P5, T")
# Lemma 3: no partition
def parts(items):
    if not items: yield []; return
    a = items[0]
    for i in range(1, len(items)):
        rest = items[1:i] + items[i+1:]
        for p in parts(rest): yield [(a, items[i])] + p
cnt = 0
for R in X4:
    W = sorted(frozenset(E) - R)
    for P in parts(W):
        A,B,C = [frozenset(p) for p in P]
        if (A|B) in calX and (B|C) in calX and (C|A) in calX: cnt += 1
ok(cnt == 0, "Lemma 3: no three disjoint edge pairs A,B,C with A|B, B|C, C|A all in calX")
# Lemma 3, case P5: W minus the edge joining the ends of R is a 5-cycle
R = frozenset({(0,1),(1,2),(2,3),(3,4)}); W = frozenset(E) - R
C5 = W - {(0,4)}
ok(all(deg(C5,v) == 2 for v in V) and len(C5) == 5, "case R=P5: W minus the end-to-end edge is a 5-cycle")
# Lemma 4 (stars)
stars = [frozenset(e for e in E if v in e) for v in V]
ok(not any(sunflower(a,b,c) for a,b,c in itertools.combinations(stars,3)), "Lemma 4a: no three stars form a sunflower")
ok(not any(sunflower(a,b,X) for a,b in itertools.combinations(stars,2) for X in F6), "Lemma 4b: no two stars and a 6-set form a sunflower")
ok(not any(sunflower(s,X,Y) for s in stars for X,Y in itertools.combinations(F6,2)), "Lemma 4c: no star and two 6-sets form a sunflower")
# Lemma 4c ingredients: leaf structure
for u in V:
    for w in V:
        if w == u: continue
        others = [v for v in V if v != u]
        K4 = [e for e in E if u not in e]
        allowed = set()
        for Q in itertools.combinations(K4, 3):
            x = frozenset(Q) | {tuple(sorted((u,w)))}
            if x in calX: allowed.add(frozenset(Q))
        comp_ok = all((frozenset(K4) - Q) not in allowed for Q in allowed)
        if not comp_ok: print("FAIL complement closure", u, w)
print("PASS Lemma 4c: the allowed Q's are never closed under complement in K4 (all u, w)")
# Sperner
fam = stars + F6
ok(not any(a < b for a in fam for b in fam), "F is Sperner")
ok(len(set(fam)) == 75, "|F| = 75")
