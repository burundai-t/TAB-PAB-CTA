#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TAB Classification — Full 90-class resolution
==============================================
The algebraic signature (15 invariants) gives 56 classes.
Burnside predicts 90. This script adds Cayley table hashing
to achieve exact classification.

Run:  python tab_classify.py
"""

import sys
import io

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

from itertools import product, permutations
from collections import defaultdict

S = ["P", "A", "B"]
M = [(R, C) for R in S for C in S]


def nu(x):
    return S[(S.index(x) + 1) % 3]


def complement(a, b):
    return [x for x in S if x != a and x != b][0]


CROSS_RULES = ["g1", "g2", "g3", "g4", "g5", "g6"]
CROSS_NAMES = {
    "g1": "column-blind", "g2": "transparent",
    "g3": "su(2)-transp", "g4": "echo",
    "g5": "diagonal", "g6": "anti-compl"
}


def sigma_equivariant_diagonals():
    diags = []
    for dp, da, db in product(S, S, S):
        d = {}
        d[("P", "P")] = dp;  d[("P", "A")] = da;  d[("P", "B")] = db
        d[("A", "A")] = nu(dp);  d[("A", "B")] = nu(da);  d[("A", "P")] = nu(db)
        d[("B", "B")] = nu(nu(dp));  d[("B", "P")] = nu(nu(da));  d[("B", "A")] = nu(nu(db))
        diags.append(((dp, da, db), d))
    return diags


def make_meet(g_rule, d_dict):
    def meet(x, y):
        R1, C1 = x
        R2, C2 = y
        if R1 == R2:
            if C1 == C2:
                return (R1, d_dict[(R1, C1)])
            else:
                return (R1, complement(C1, C2))
        if g_rule == "g1":   col = R2
        elif g_rule == "g2": col = C2
        elif g_rule == "g3": col = complement(C1, C2) if C1 != C2 else C1
        elif g_rule == "g4": col = C1
        elif g_rule == "g5": col = R1
        elif g_rule == "g6": col = complement(R2, C2) if R2 != C2 else R2
        return (R1, col)
    return meet


def build_cayley(meet_fn):
    return {(x, y): meet_fn(x, y) for x in M for y in M}


def canonical_cayley(table):
    """Compute the canonical form of a Cayley table under S3 action.

    For each of the 6 permutations of {P,A,B}, relabel the table
    and return the LEXICOGRAPHICALLY SMALLEST result. Two magmas
    are isomorphic iff they have the same canonical form.
    """
    s3_perms = list(permutations(S))

    canonical = None
    for perm in s3_perms:
        pm = dict(zip(S, perm))
        # Relabel: (R,C) -> (pm[R], pm[C])
        relabeled = {}
        for (x, y), z in table.items():
            new_x = (pm[x[0]], pm[x[1]])
            new_y = (pm[y[0]], pm[y[1]])
            new_z = (pm[z[0]], pm[z[1]])
            relabeled[(new_x, new_y)] = new_z

        # Convert to a hashable sorted tuple
        form = tuple(
            (x, y, relabeled[(x, y)])
            for x in M for y in M
        )

        if canonical is None or form < canonical:
            canonical = form

    return canonical


def quick_invariants(table):
    """Quick invariant tuple for pre-grouping (avoids expensive canonical form
    for obviously different magmas)."""
    assoc = sum(1 for a, b, c in product(M, M, M)
                if table[(table[(a, b)], c)] == table[(a, table[(b, c)])])
    comm = sum(1 for x, y in product(M, M) if table[(x, y)] == table[(y, x)])
    ralt = sum(1 for x, y in product(M, M)
               if table[(table[(y, x)], x)] == table[(y, table[(x, x)])])
    lalt = sum(1 for x, y in product(M, M)
               if table[(x, table[(x, y)])] == table[(table[(x, x)], y)])
    flex = sum(1 for x, y in product(M, M)
               if table[(table[(x, y)], x)] == table[(x, table[(y, x)])])
    idemp = sum(1 for x in M if table[(x, x)] == x)
    return (assoc, comm, ralt, lalt, flex, idemp)


def main():
    print()
    print("TAB Classification - Exact 90-class resolution")
    print("=" * 70)
    print()

    diags = sigma_equivariant_diagonals()

    # Build all 162 magmas
    all_magmas = []
    for g in CROSS_RULES:
        for (dp, da, db), delta_dict in diags:
            meet = make_meet(g, delta_dict)
            table = build_cayley(meet)
            all_magmas.append({
                "rule": g,
                "delta": (dp, da, db),
                "table": table,
            })

    print("Total magmas: %d" % len(all_magmas))

    # Step 1: Group by quick invariants (pre-filter)
    print("\nStep 1: Pre-grouping by quick invariants...")
    inv_groups = defaultdict(list)
    for i, m in enumerate(all_magmas):
        inv = quick_invariants(m["table"])
        m["inv"] = inv
        inv_groups[inv].append(i)

    print("  Pre-groups: %d" % len(inv_groups))
    ambiguous = sum(1 for g in inv_groups.values() if len(g) > 1)
    print("  Groups with >1 member: %d (need canonical form)" % ambiguous)

    # Step 2: Within each pre-group, compute canonical forms
    print("\nStep 2: Computing canonical forms within ambiguous groups...")
    canon_classes = defaultdict(list)
    computed = 0
    for inv, indices in inv_groups.items():
        if len(indices) == 1:
            # Singleton: canonical form not needed, use inv as key
            idx = indices[0]
            canon_classes[("singleton", inv)].append(idx)
        else:
            # Compute canonical form for each member
            for idx in indices:
                cf = canonical_cayley(all_magmas[idx]["table"])
                canon_classes[("canon", cf)].append(idx)
                computed += 1

    print("  Canonical forms computed: %d" % computed)
    print("  Total isomorphism classes: %d" % len(canon_classes))
    print()

    # Step 3: Report
    print("=" * 70)
    print("CLASSIFICATION RESULT")
    print("=" * 70)
    print()

    # Sort classes by (assoc, size)
    class_info = []
    for key, indices in canon_classes.items():
        rep = all_magmas[indices[0]]
        assoc = rep["inv"][0]
        rules = sorted(set(all_magmas[i]["rule"] for i in indices))
        deltas = [all_magmas[i]["delta"] for i in indices]
        class_info.append({
            "size": len(indices),
            "assoc": assoc,
            "inv": rep["inv"],
            "rules": rules,
            "deltas": deltas,
            "indices": indices
        })

    class_info.sort(key=lambda c: (c["assoc"], c["size"]))

    print("  %4s  %5s  %4s  %4s  %4s  %4s  %3s  %4s  %s" %
          ("Cls", "Assoc", "Comm", "Ralt", "Lalt", "Flex", "Id", "Size", "Rules"))
    for i, c in enumerate(class_info):
        a, cm, r, l, f, idm = c["inv"]
        rules_str = "+".join(c["rules"])
        print("  %4d  %5d  %4d  %4d  %4d  %4d  %3d  %4d  [%s]" %
              (i+1, a, cm, r, l, f, idm, c["size"], rules_str))

    print()
    print("  Total classes: %d" % len(class_info))
    print()

    # Orbit size distribution
    sizes = [c["size"] for c in class_info]
    size_dist = Counter(sizes)
    print("  Orbit size distribution:")
    for s in sorted(size_dist.keys()):
        print("    Size %d: %d classes" % (s, size_dist[s]))
    total_check = sum(s * count for s, count in size_dist.items())
    print("    Total magmas: %d (check: %s)" %
          (total_check, "OK" if total_check == 162 else "FAIL"))
    print()

    # Singletons
    singletons = [c for c in class_info if c["size"] == 1]
    print("  Singleton classes (unique magmas):")
    for c in singletons:
        rep_idx = c["indices"][0]
        rep = all_magmas[rep_idx]
        print("    %s, delta=(%s,%s,%s), assoc=%d" %
              (CROSS_NAMES[rep["rule"]], rep["delta"][0], rep["delta"][1],
               rep["delta"][2], c["assoc"]))
    print()

    # CTA-compatible classes
    print("  CTA-compatible classes:")
    for i, c in enumerate(class_info):
        # Check if any member is CTA-compatible
        cta_members = []
        for idx in c["indices"]:
            m = all_magmas[idx]
            dp, da, db = m["delta"]
            if dp == da == db:
                t = m["table"]
                idemp_count = sum(1 for x in M if t[(x,x)] == x)
                i1 = all(t[(t[(x,y)], t[(x,y)])] == t[(x,x)]
                         for x, y in product(M, M))
                if idemp_count == 3 and i1:
                    cta_members.append(idx)
        if cta_members:
            rep = all_magmas[cta_members[0]]
            print("    Class %d: %s, assoc=%d, size=%d (%d CTA-compat)" %
                  (i+1, "+".join(c["rules"]), c["assoc"], c["size"], len(cta_members)))

    print()
    print("=" * 70)
    print("CLASSIFICATION COMPLETE")
    print("=" * 70)


from collections import Counter

if __name__ == "__main__":
    main()
