#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TAB Graph: Adjacency structure of 90 isomorphism classes
========================================================
Two classes are "adjacent" if they have representatives differing
by one mutation (single delta-component change OR cross-rule change).

Run:  python tab_graph.py
"""

import sys
import io

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

from itertools import product, permutations
from collections import defaultdict, Counter

S = ["P", "A", "B"]
M = [(R, C) for R in S for C in S]

def nu(x): return S[(S.index(x) + 1) % 3]
def complement(a, b): return [x for x in S if x != a and x != b][0]

CROSS_RULES = ["g1", "g2", "g3", "g4", "g5", "g6"]

def sigma_equivariant_diagonals():
    diags = []
    for dp, da, db in product(S, S, S):
        d = {}
        d[("P","P")]=dp; d[("P","A")]=da; d[("P","B")]=db
        d[("A","A")]=nu(dp); d[("A","B")]=nu(da); d[("A","P")]=nu(db)
        d[("B","B")]=nu(nu(dp)); d[("B","P")]=nu(nu(da)); d[("B","A")]=nu(nu(db))
        diags.append(((dp,da,db), d))
    return diags

def make_meet(g, d):
    def meet(x, y):
        R1,C1 = x; R2,C2 = y
        if R1==R2:
            return (R1, d[(R1,C1)]) if C1==C2 else (R1, complement(C1,C2))
        if g=="g1": col=R2
        elif g=="g2": col=C2
        elif g=="g3": col=complement(C1,C2) if C1!=C2 else C1
        elif g=="g4": col=C1
        elif g=="g5": col=R1
        elif g=="g6": col=complement(R2,C2) if R2!=C2 else R2
        return (R1, col)
    return meet

def build_cayley(meet_fn):
    return {(x,y): meet_fn(x,y) for x in M for y in M}

def canonical_cayley(table):
    s3_perms = list(permutations(S))
    canonical = None
    for perm in s3_perms:
        pm = dict(zip(S, perm))
        form = tuple(
            (x, y, (pm[table[(x,y)][0]], pm[table[(x,y)][1]]))
            for x in [(pm[e[0]], pm[e[1]]) for e in M]
            for y in [(pm[e[0]], pm[e[1]]) for e in M]
        )
        if canonical is None or form < canonical:
            canonical = form
    return canonical

def quick_sig(table):
    assoc = sum(1 for a,b,c in product(M,M,M)
                if table[(table[(a,b)],c)] == table[(a,table[(b,c)])])
    ralt = sum(1 for x,y in product(M,M)
               if table[(table[(y,x)],x)] == table[(y,table[(x,x)])])
    idemp = sum(1 for x in M if table[(x,x)] == x)
    return (assoc, ralt, idemp)


def main():
    print()
    print("TAB Graph: Adjacency of 90 Classes")
    print("=" * 70)
    print()

    diags = sigma_equivariant_diagonals()
    delta_list = [(dp,da,db) for (dp,da,db), _ in diags]

    # Build all 162 magmas and assign class IDs
    magma_index = {}  # (rule, delta) -> class_id
    class_members = defaultdict(list)  # class_id -> [(rule, delta)]

    # First pass: compute canonical forms for class assignment
    print("Computing canonical forms for 162 magmas...")
    canon_to_class = {}
    class_counter = 0

    all_tables = {}
    for g in CROSS_RULES:
        for (dp,da,db), dd in diags:
            table = build_cayley(make_meet(g, dd))
            all_tables[(g, (dp,da,db))] = table
            cf = canonical_cayley(table)
            if cf not in canon_to_class:
                canon_to_class[cf] = class_counter
                class_counter += 1
            cid = canon_to_class[cf]
            magma_index[(g, (dp,da,db))] = cid
            class_members[cid].append((g, (dp,da,db)))

    n_classes = len(canon_to_class)
    print("Classes: %d\n" % n_classes)

    # Compute adjacency: two classes are adjacent if any of their
    # members differ by one mutation
    print("Computing adjacency...")

    edges = set()

    for g in CROSS_RULES:
        for i, (dp,da,db) in enumerate(delta_list):
            cid1 = magma_index[(g, (dp,da,db))]

            # Delta mutations: change one of dp, da, db
            for new_dp in S:
                if new_dp != dp:
                    new_delta = (new_dp, da, db)
                    if new_delta in [(d[0],d[1],d[2]) for d in delta_list]:
                        cid2 = magma_index[(g, new_delta)]
                        if cid1 != cid2:
                            edges.add((min(cid1,cid2), max(cid1,cid2)))
            for new_da in S:
                if new_da != da:
                    new_delta = (dp, new_da, db)
                    if new_delta in [(d[0],d[1],d[2]) for d in delta_list]:
                        cid2 = magma_index[(g, new_delta)]
                        if cid1 != cid2:
                            edges.add((min(cid1,cid2), max(cid1,cid2)))
            for new_db in S:
                if new_db != db:
                    new_delta = (dp, da, new_db)
                    if new_delta in [(d[0],d[1],d[2]) for d in delta_list]:
                        cid2 = magma_index[(g, new_delta)]
                        if cid1 != cid2:
                            edges.add((min(cid1,cid2), max(cid1,cid2)))

            # Cross-rule mutations: change g
            for g2 in CROSS_RULES:
                if g2 != g:
                    cid2 = magma_index[(g2, (dp,da,db))]
                    if cid1 != cid2:
                        edges.add((min(cid1,cid2), max(cid1,cid2)))

    print("Edges: %d\n" % len(edges))

    # Degree distribution
    degree = Counter()
    for c1, c2 in edges:
        degree[c1] += 1
        degree[c2] += 1

    # Include isolated vertices
    for cid in range(n_classes):
        if cid not in degree:
            degree[cid] = 0

    deg_dist = Counter(degree.values())

    print("=" * 70)
    print("1. DEGREE DISTRIBUTION")
    print("=" * 70)
    print()
    print("  Degree  Count")
    for d in sorted(deg_dist.keys()):
        print("  %6d  %5d" % (d, deg_dist[d]))
    print()
    print("  Min degree: %d" % min(degree.values()))
    print("  Max degree: %d" % max(degree.values()))
    print("  Mean degree: %.1f" % (sum(degree.values()) / n_classes))
    print("  Isolated (deg=0): %d" % deg_dist.get(0, 0))

    # ================================================================
    print()
    print("=" * 70)
    print("2. INTRA-RULE vs INTER-RULE EDGES")
    print("=" * 70)
    print()

    # Classify edges
    intra = 0
    inter = 0
    for c1, c2 in edges:
        rules1 = set(g for g, d in class_members[c1])
        rules2 = set(g for g, d in class_members[c2])
        if rules1 & rules2:
            intra += 1
        else:
            inter += 1

    # Some edges might connect classes that share a rule AND differ in rule
    both = len(edges) - intra - inter
    print("  Intra-rule edges (same g, different delta): %d" % intra)
    print("  Inter-rule edges (different g, same delta): %d" % inter)
    print("  Total edges: %d" % len(edges))

    # ================================================================
    print()
    print("=" * 70)
    print("3. CONNECTED COMPONENTS")
    print("=" * 70)
    print()

    # Simple BFS for connected components
    adj = defaultdict(set)
    for c1, c2 in edges:
        adj[c1].add(c2)
        adj[c2].add(c1)

    visited = set()
    components = []
    for start in range(n_classes):
        if start in visited:
            continue
        comp = set()
        queue = [start]
        while queue:
            node = queue.pop(0)
            if node in visited:
                continue
            visited.add(node)
            comp.add(node)
            for nb in adj[node]:
                if nb not in visited:
                    queue.append(nb)
        components.append(comp)

    print("  Connected components: %d" % len(components))
    comp_sizes = sorted([len(c) for c in components], reverse=True)
    print("  Component sizes: %s" % comp_sizes)

    if len(components) == 1:
        print("  The graph is CONNECTED: every class can be reached from every other.")
    else:
        print("  Isolated components:")
        for i, comp in enumerate(components):
            if len(comp) <= 5:
                rules = set()
                for cid in comp:
                    for g, d in class_members[cid]:
                        rules.add(g)
                print("    Component %d (size %d): rules=%s" %
                      (i, len(comp), sorted(rules)))

    # ================================================================
    print()
    print("=" * 70)
    print("4. HIGH-DEGREE AND LOW-DEGREE CLASSES")
    print("=" * 70)
    print()

    # Get assoc for each class
    class_assoc = {}
    for cid in range(n_classes):
        rep_g, rep_d = class_members[cid][0]
        table = all_tables[(rep_g, rep_d)]
        a = sum(1 for x,y,z in product(M,M,M)
                if table[(table[(x,y)],z)] == table[(x,table[(y,z)])])
        class_assoc[cid] = a

    print("  Top 10 by degree:")
    for cid in sorted(degree, key=degree.get, reverse=True)[:10]:
        rep_g, rep_d = class_members[cid][0]
        print("    Class %2d: deg=%2d, assoc=%d, rule=%s, size=%d" %
              (cid, degree[cid], class_assoc[cid], rep_g, len(class_members[cid])))

    print()
    print("  Bottom 10 by degree:")
    for cid in sorted(degree, key=degree.get)[:10]:
        rep_g, rep_d = class_members[cid][0]
        print("    Class %2d: deg=%2d, assoc=%d, rule=%s, size=%d" %
              (cid, degree[cid], class_assoc[cid], rep_g, len(class_members[cid])))

    # ================================================================
    print()
    print("=" * 70)
    print("5. DEGREE vs ASSOCIATIVITY CORRELATION")
    print("=" * 70)
    print()

    # Group by assoc, show mean degree
    assoc_degrees = defaultdict(list)
    for cid in range(n_classes):
        assoc_degrees[class_assoc[cid]].append(degree[cid])

    print("  Assoc  Classes  MeanDeg  MinDeg  MaxDeg")
    for a in sorted(assoc_degrees.keys()):
        degs = assoc_degrees[a]
        print("  %5d  %5d    %5.1f    %5d    %5d" %
              (a, len(degs), sum(degs)/len(degs), min(degs), max(degs)))

    print()
    print("=" * 70)
    print("GRAPH ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
