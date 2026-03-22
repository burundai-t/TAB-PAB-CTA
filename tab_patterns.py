#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TAB Patterns — Deep analysis of the 162-magma / 90-class space
==============================================================
Searches for structural patterns, formulas, and new constants.

Run:  python tab_patterns.py
"""

import sys
import io

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

from itertools import product, permutations
from collections import defaultdict, Counter
from math import gcd, log

S = ["P", "A", "B"]
M = [(R, C) for R in S for C in S]


def nu(x):
    return S[(S.index(x) + 1) % 3]


def complement(a, b):
    return [x for x in S if x != a and x != b][0]


CROSS_RULES = ["g1", "g2", "g3", "g4", "g5", "g6"]
CROSS_NAMES = {
    "g1": "col-blind", "g2": "transp",
    "g3": "su2-tr", "g4": "echo",
    "g5": "diag", "g6": "anti-c"
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


def compute_full(table):
    fp = {}
    fp["assoc"] = sum(1 for a, b, c in product(M, M, M)
                      if table[(table[(a, b)], c)] == table[(a, table[(b, c)])])
    fp["ralt"] = sum(1 for x, y in product(M, M)
                     if table[(table[(y, x)], x)] == table[(y, table[(x, x)])])
    fp["lalt"] = sum(1 for x, y in product(M, M)
                     if table[(x, table[(x, y)])] == table[(table[(x, x)], y)])
    fp["flex"] = sum(1 for x, y in product(M, M)
                     if table[(table[(x, y)], x)] == table[(x, table[(y, x)])])
    fp["idemp"] = sum(1 for x in M if table[(x, x)] == x)

    # Row-pattern breakdown
    pats = {"RRR": 0, "RRS": 0, "RSR": 0, "SRR": 0, "DIST": 0}
    for a, b, c in product(M, M, M):
        ra, rb, rc = a[0], b[0], c[0]
        if ra == rb == rc:      pat = "RRR"
        elif ra == rb != rc:    pat = "RRS"
        elif ra != rb == rc:    pat = "SRR"
        elif ra == rc != rb:    pat = "RSR"
        else:                   pat = "DIST"
        if table[(table[(a, b)], c)] == table[(a, table[(b, c)])]:
            pats[pat] += 1
    fp.update({"a_" + k: v for k, v in pats.items()})

    # tau_AB fixed points
    tau = {"P": "P", "A": "B", "B": "A"}
    def tau_el(e):
        return (tau[e[0]], tau[e[1]])
    fp["fix_tau"] = sum(1 for x, y in product(M, M)
                        if tau_el(table[(x, y)]) == table[(tau_el(x), tau_el(y))])

    # Diagonal fixed-point count
    fp["diag_fix"] = sum(1 for R in S for C in S
                         if table[((R, C), (R, C))][1] == C)
    return fp


def main():
    print()
    print("TAB Patterns - Deep Analysis")
    print("=" * 70)
    print()

    diags = sigma_equivariant_diagonals()
    all_magmas = []
    for g in CROSS_RULES:
        for (dp, da, db), delta_dict in diags:
            meet = make_meet(g, delta_dict)
            table = build_cayley(meet)
            fp = compute_full(table)
            all_magmas.append({
                "rule": g, "delta": (dp, da, db),
                "table": table, "fp": fp
            })

    print("Total: %d magmas\n" % len(all_magmas))

    # =================================================================
    print("=" * 70)
    print("A. ASSOCIATIVITY FORMULA SEARCH")
    print("=" * 70)
    print()
    print("  Is Assoc = f(RRR, RRS, RSR, SRR, DIST)?")
    print()

    # For each magma, verify: assoc = sum of pattern parts
    for m in all_magmas:
        fp = m["fp"]
        s = fp["a_RRR"] + fp["a_RRS"] + fp["a_RSR"] + fp["a_SRR"] + fp["a_DIST"]
        if s != fp["assoc"]:
            print("  ERROR: sum mismatch for %s %s" % (m["rule"], m["delta"]))
            break
    else:
        print("  OK: Assoc = RRR + RRS + RSR + SRR + DIST (verified for all 162)")

    print()
    # Pattern contributions by cross-rule
    print("  Pattern contributions by cross-rule (for constant delta = (P,P,P)):")
    print("  %-8s %5s %5s %5s %5s %5s  %5s" %
          ("Rule", "RRR", "RRS", "RSR", "SRR", "DIST", "Total"))
    for g in CROSS_RULES:
        for m in all_magmas:
            if m["rule"] == g and m["delta"] == ("P", "P", "P"):
                fp = m["fp"]
                print("  %-8s %5d %5d %5d %5d %5d  %5d" %
                      (CROSS_NAMES[g],
                       fp["a_RRR"], fp["a_RRS"], fp["a_RSR"],
                       fp["a_SRR"], fp["a_DIST"], fp["assoc"]))
                break

    print()
    print("  Pattern contributions for ALL 27 delta values (g1 = column-blind):")
    print("  %-12s %5s %5s %5s %5s %5s  %5s" %
          ("delta", "RRR", "RRS", "RSR", "SRR", "DIST", "Total"))
    g1_magmas = sorted([m for m in all_magmas if m["rule"] == "g1"],
                       key=lambda m: m["fp"]["assoc"])
    for m in g1_magmas:
        fp = m["fp"]
        d = m["delta"]
        print("  (%s,%s,%s)    %5d %5d %5d %5d %5d  %5d" %
              (d[0], d[1], d[2],
               fp["a_RRR"], fp["a_RRS"], fp["a_RSR"],
               fp["a_SRR"], fp["a_DIST"], fp["assoc"]))

    # =================================================================
    print()
    print("=" * 70)
    print("B. CROSS-RULE PAIRS AND SYMMETRIES")
    print("=" * 70)
    print()

    # For each cross-rule pair, check if there's a consistent relationship
    print("  Associativity relationship between rule pairs (delta = (P,P,P)):")
    print("  (comparing at same delta)")
    print()
    ref = {}
    for m in all_magmas:
        if m["delta"] == ("P", "P", "P"):
            ref[m["rule"]] = m["fp"]["assoc"]

    print("  %-10s" % "" + "  ".join("%-8s" % CROSS_NAMES[g] for g in CROSS_RULES))
    for g1 in CROSS_RULES:
        row = "  %-10s" % CROSS_NAMES[g1]
        for g2 in CROSS_RULES:
            if g1 == g2:
                row += "  %-8s" % "---"
            else:
                row += "  %-8d" % (ref[g2] - ref[g1])
        print(row)

    # =================================================================
    print()
    print("=" * 70)
    print("C. RRR VALUES AND FIBER ASSOCIATIVITY")
    print("=" * 70)
    print()

    # RRR depends only on the diagonal (fiber), not on the cross-rule
    print("  RRR values (should depend only on delta, not on cross-rule):")
    rrr_by_delta = defaultdict(set)
    for m in all_magmas:
        rrr_by_delta[m["delta"]].add(m["fp"]["a_RRR"])

    rrr_independent = all(len(v) == 1 for v in rrr_by_delta.values())
    print("  RRR independent of cross-rule: %s" % rrr_independent)
    print()

    if rrr_independent:
        print("  RRR by delta (P-row params):")
        print("  %-12s  %5s  %s" % ("delta", "RRR", "RRR/3 (= fiber assoc)"))
        for delta in sorted(rrr_by_delta.keys()):
            rrr = list(rrr_by_delta[delta])[0]
            print("  (%s,%s,%s)     %5d  %5d/27" % (delta[0], delta[1], delta[2], rrr, rrr//3))

    # =================================================================
    print()
    print("=" * 70)
    print("D. SRR PATTERN — THE BUFFER THEOREM")
    print("=" * 70)
    print()

    # SRR = (R_a != R_b = R_c) — the buffer pattern
    # In PAB (g1): SRR is always associative (162/162)
    # What about other rules?
    print("  SRR values by cross-rule:")
    for g in CROSS_RULES:
        vals = sorted(set(m["fp"]["a_SRR"] for m in all_magmas if m["rule"] == g))
        print("  %-10s: %s" % (CROSS_NAMES[g], vals))

    print()
    print("  Which rules have SRR = 162 (always associative buffer)?")
    for g in CROSS_RULES:
        all_162 = all(m["fp"]["a_SRR"] == 162 for m in all_magmas if m["rule"] == g)
        print("  %-10s: %s" % (CROSS_NAMES[g], "YES (SRR=162 always)" if all_162 else "no"))

    # =================================================================
    print()
    print("=" * 70)
    print("E. DIST PATTERN — ALL-DISTINCT TYPES")
    print("=" * 70)
    print()

    print("  DIST values by cross-rule:")
    for g in CROSS_RULES:
        vals = sorted(set(m["fp"]["a_DIST"] for m in all_magmas if m["rule"] == g))
        print("  %-10s: %s" % (CROSS_NAMES[g], vals))

    print()
    print("  DIST = 162 (always assoc) vs DIST = 0 (never assoc):")
    for g in CROSS_RULES:
        d162 = sum(1 for m in all_magmas if m["rule"] == g and m["fp"]["a_DIST"] == 162)
        d0 = sum(1 for m in all_magmas if m["rule"] == g and m["fp"]["a_DIST"] == 0)
        d54 = sum(1 for m in all_magmas if m["rule"] == g and m["fp"]["a_DIST"] == 54)
        print("  %-10s: DIST=162: %2d  DIST=54: %2d  DIST=0: %2d" %
              (CROSS_NAMES[g], d162, d54, d0))

    # =================================================================
    print()
    print("=" * 70)
    print("F. THE DIAGONAL FIXED-POINT COUNT")
    print("=" * 70)
    print()

    # diag_fix = number of (R,C) where delta(R,C) = C (= fixed points of delta)
    print("  diag_fix = |{(R,C) : delta(R,C) = C}| (diagonal idempotent count)")
    print()
    df_vals = sorted(set(m["fp"]["diag_fix"] for m in all_magmas))
    print("  Distinct values: %s" % df_vals)
    print("  (These are 0, 3, 6, 9 = multiples of n = 3)")
    print()

    # Relationship between diag_fix and other invariants
    print("  Relationship: diag_fix vs idemp, RRR, assoc (for g1):")
    print("  %-5s  %-5s  %-5s  %-6s" % ("dfix", "idemp", "RRR", "assoc"))
    seen = set()
    for m in sorted(all_magmas, key=lambda m: m["fp"]["diag_fix"]):
        if m["rule"] != "g1":
            continue
        fp = m["fp"]
        key = (fp["diag_fix"], fp["idemp"], fp["a_RRR"])
        if key in seen:
            continue
        seen.add(key)
        print("  %-5d  %-5d  %-5d  %-6d" %
              (fp["diag_fix"], fp["idemp"], fp["a_RRR"], fp["assoc"]))

    # =================================================================
    print()
    print("=" * 70)
    print("G. ASSOCIATIVITY AS FUNCTION OF DELTA (within each rule)")
    print("=" * 70)
    print()

    # For g1 (column-blind), assoc depends on delta through RRR + SRR + RSR + RRS + DIST
    # SRR is always 162 for g1. DIST is always 0 for g1.
    # So: assoc(g1, delta) = RRR(delta) + RRS(delta) + RSR(delta) + 162 + 0
    #                       = RRR + RRS + RSR + 162
    print("  For g1 (column-blind): DIST=0, SRR=162 always.")
    print("  So: Assoc = RRR + RRS + RSR + 162")
    print()

    print("  Verify:")
    for m in all_magmas:
        if m["rule"] != "g1":
            continue
        fp = m["fp"]
        pred = fp["a_RRR"] + fp["a_RRS"] + fp["a_RSR"] + 162
        if pred != fp["assoc"]:
            print("  FAIL: delta=%s, pred=%d, actual=%d" % (m["delta"], pred, fp["assoc"]))
            break
    else:
        print("  OK: Assoc(g1) = RRR + RRS + RSR + 162 for all 27 deltas")

    # Same for other rules
    print()
    for g in CROSS_RULES:
        members = [m for m in all_magmas if m["rule"] == g]
        # Check if SRR and DIST are constant
        srr_vals = set(m["fp"]["a_SRR"] for m in members)
        dist_vals = set(m["fp"]["a_DIST"] for m in members)
        srr_const = len(srr_vals) == 1
        dist_const = len(dist_vals) == 1
        if srr_const and dist_const:
            srr_v = list(srr_vals)[0]
            dist_v = list(dist_vals)[0]
            print("  %s: SRR=%d (const), DIST=%d (const)" % (CROSS_NAMES[g], srr_v, dist_v))
            print("    => Assoc = RRR + RRS + RSR + %d" % (srr_v + dist_v))
        else:
            print("  %s: SRR varies (%s), DIST varies (%s)" %
                  (CROSS_NAMES[g], sorted(srr_vals), sorted(dist_vals)))

    # =================================================================
    print()
    print("=" * 70)
    print("H. KEY NUMBERS IN THE TAB SPECTRUM")
    print("=" * 70)
    print()

    assoc_spectrum = sorted(set(m["fp"]["assoc"] for m in all_magmas))
    print("  All 28 assoc values as fractions of 729:")
    for a in assoc_spectrum:
        g = gcd(a, 729)
        print("  %3d/729 = %3d/%3d = %.6f" % (a, a//g, 729//g, a/729))

    print()
    print("  Differences between consecutive values:")
    diffs = [assoc_spectrum[i+1] - assoc_spectrum[i] for i in range(len(assoc_spectrum)-1)]
    print("  %s" % diffs)
    print("  Distinct diffs: %s" % sorted(set(diffs)))

    print()
    print("  GCD of all assoc values: %d" % gcd(*assoc_spectrum) if len(assoc_spectrum) > 1 else "N/A")
    print("  All divisible by 3: %s" % all(a % 3 == 0 for a in assoc_spectrum))
    print("  All divisible by 6: %s" % all(a % 6 == 0 for a in assoc_spectrum))

    # =================================================================
    print()
    print("=" * 70)
    print("I. ECHO (g4) STRUCTURE — THE MAXIMUM")
    print("=" * 70)
    print()

    echo_magmas = [m for m in all_magmas if m["rule"] == "g4"]
    print("  Echo rule: 27 magmas")
    print("  %-12s %5s %5s %5s %5s %5s  %5s  %3s" %
          ("delta", "RRR", "RRS", "RSR", "SRR", "DIST", "Total", "dfix"))
    for m in sorted(echo_magmas, key=lambda m: m["fp"]["assoc"]):
        fp = m["fp"]
        d = m["delta"]
        print("  (%s,%s,%s)    %5d %5d %5d %5d %5d  %5d  %3d" %
              (d[0], d[1], d[2],
               fp["a_RRR"], fp["a_RRS"], fp["a_RSR"],
               fp["a_SRR"], fp["a_DIST"], fp["assoc"], fp["diag_fix"]))

    # =================================================================
    print()
    print("=" * 70)
    print("PATTERN ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
