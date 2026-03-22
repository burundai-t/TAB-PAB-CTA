#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TAB Spectrum: Analysis of the 28 associativity values
=====================================================
Run:  python tab_spectrum.py
"""

import sys
import io

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

from itertools import product
from collections import defaultdict, Counter
from math import gcd

S = ["P", "A", "B"]
M = [(R, C) for R in S for C in S]

def nu(x): return S[(S.index(x) + 1) % 3]
def complement(a, b): return [x for x in S if x != a and x != b][0]

CROSS_RULES = ["g1", "g2", "g3", "g4", "g5", "g6"]
CROSS_NAMES = {"g1": "col-blind", "g2": "transp", "g3": "su2-tr",
               "g4": "echo", "g5": "diag", "g6": "anti-c"}

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

def count_assoc(table):
    return sum(1 for a,b,c in product(M,M,M)
               if table[(table[(a,b)],c)] == table[(a,table[(b,c)])])

def main():
    print()
    print("TAB Spectrum Analysis")
    print("=" * 70)
    print()

    diags = sigma_equivariant_diagonals()
    all_magmas = []
    for g in CROSS_RULES:
        for (dp,da,db), dd in diags:
            table = build_cayley(make_meet(g, dd))
            a = count_assoc(table)
            all_magmas.append({"rule": g, "delta": (dp,da,db), "assoc": a})

    spectrum = sorted(set(m["assoc"] for m in all_magmas))
    print("28 values: %s\n" % spectrum)

    # ================================================================
    print("=" * 70)
    print("1. FACTORIZATION OF SPECTRUM VALUES")
    print("=" * 70)
    print()
    print("  %5s  %10s  %6s  %s" % ("Value", "Fraction", "Mod 9", "Factorization"))
    for v in spectrum:
        g = gcd(v, 729)
        # Factor v
        n = v
        factors = []
        for p in [2,3,5,7,11,13,17,19,23]:
            while n % p == 0:
                factors.append(p)
                n //= p
            if n == 1:
                break
        if n > 1:
            factors.append(n)
        print("  %5d  %4d/%-4d  %6d  %s" %
              (v, v//g, 729//g, v % 9, " x ".join(str(f) for f in factors)))

    # ================================================================
    print()
    print("=" * 70)
    print("2. DIFFERENCES AND GAPS")
    print("=" * 70)
    print()

    diffs = [spectrum[i+1] - spectrum[i] for i in range(len(spectrum)-1)]
    print("  Consecutive differences: %s" % diffs)
    print("  Distinct diffs: %s" % sorted(set(diffs)))
    print("  All diffs divisible by 3: %s" % all(d % 3 == 0 for d in diffs))
    print("  Diffs / 3: %s" % [d//3 for d in diffs])
    print()

    # Factor the diffs
    print("  Diffs factored:")
    for d in sorted(set(diffs)):
        n = d
        factors = []
        for p in [2,3,5,7,11,13]:
            while n % p == 0:
                factors.append(p)
                n //= p
        if n > 1:
            factors.append(n)
        count = diffs.count(d)
        print("    %3d = %-15s (appears %dx)" %
              (d, " x ".join(str(f) for f in factors), count))

    # ================================================================
    print()
    print("=" * 70)
    print("3. SPECTRUM BY CROSS-RULE")
    print("=" * 70)
    print()

    for g in CROSS_RULES:
        vals = sorted(set(m["assoc"] for m in all_magmas if m["rule"] == g))
        rng = vals[-1] - vals[0]
        print("  %s: %d values, range [%d, %d], span=%d" %
              (CROSS_NAMES[g], len(vals), vals[0], vals[-1], rng))
        print("    values: %s" % vals)
        if len(vals) > 1:
            d = [vals[i+1]-vals[i] for i in range(len(vals)-1)]
            print("    diffs:  %s" % d)
        print()

    # ================================================================
    print()
    print("=" * 70)
    print("4. OVERLAPS BETWEEN CROSS-RULES")
    print("=" * 70)
    print()

    rule_sets = {}
    for g in CROSS_RULES:
        rule_sets[g] = set(m["assoc"] for m in all_magmas if m["rule"] == g)

    print("  Pairwise overlaps (shared assoc values):")
    for i, g1 in enumerate(CROSS_RULES):
        for g2 in CROSS_RULES[i+1:]:
            overlap = rule_sets[g1] & rule_sets[g2]
            if overlap:
                print("    %s & %s: %s" %
                      (CROSS_NAMES[g1], CROSS_NAMES[g2], sorted(overlap)))

    print()
    # Which values appear in multiple rules?
    val_rules = defaultdict(set)
    for m in all_magmas:
        val_rules[m["assoc"]].add(m["rule"])

    multi = {v: rules for v, rules in val_rules.items() if len(rules) > 1}
    print("  Values shared by multiple rules:")
    for v in sorted(multi.keys()):
        print("    %d: %s" % (v, ", ".join(CROSS_NAMES[r] for r in sorted(multi[v]))))

    single = {v: rules for v, rules in val_rules.items() if len(rules) == 1}
    print()
    print("  Values unique to one rule: %d out of %d" % (len(single), len(spectrum)))

    # ================================================================
    print()
    print("=" * 70)
    print("5. ARITHMETIC STRUCTURE")
    print("=" * 70)
    print()

    # Check: are values of the form 3k where k has a pattern?
    k_vals = [v // 3 for v in spectrum]
    print("  Values / 3: %s" % k_vals)
    print()

    # Check modular structure
    for mod in [6, 9, 18, 27]:
        residues = sorted(set(v % mod for v in spectrum))
        print("  Values mod %d: %s (%d classes)" % (mod, residues, len(residues)))

    # ================================================================
    print()
    print("=" * 70)
    print("6. COMPLEMENT PAIRS: v + v' = 729?")
    print("=" * 70)
    print()

    spec_set = set(spectrum)
    print("  Checking if v and 729-v are both in the spectrum:")
    complement_pairs = []
    for v in spectrum:
        cp = 729 - v
        if cp in spec_set and cp >= v:
            complement_pairs.append((v, cp))
            print("    %d + %d = 729  (%s)" %
                  (v, cp, "self" if v == cp else "pair"))

    print()
    if not complement_pairs:
        print("  No complement pairs found.")
    print("  Missing complements:")
    for v in spectrum:
        if 729 - v not in spec_set:
            print("    %d -> 729-%d = %d (NOT in spectrum)" % (v, v, 729-v))

    # ================================================================
    print()
    print("=" * 70)
    print("7. KNOWN CTA NUMBERS IN SPECTRUM")
    print("=" * 70)
    print()

    # Express each value in terms of CTA structural numbers
    n = 3
    cta = {
        "n": 3, "n2": 9, "n2-1": 8, "n2-n": 6, "n2-n-1": 5,
        "n2-n+1": 7, "n2+n+1": 13, "n3": 27, "n4": 81,
        "n6": 729, "(n-1)/n2 * n6": 162, "n4*(n2-n+1)/n2": 567
    }

    print("  Spectrum values as CTA expressions:")
    for v in spectrum:
        exprs = []
        # Try simple multiples
        if v % 81 == 0: exprs.append("%d * 81" % (v//81))
        if v % 27 == 0: exprs.append("%d * 27" % (v//27))
        if v % 9 == 0: exprs.append("%d * 9" % (v//9))
        # Known values
        if v == 162: exprs.append("2/9 * 729 = A_min * n^6")
        if v == 219: exprs.append("73/243 * 729 = PAB")
        if v == 243: exprs.append("1/3 * 729 = n^5")
        if v == 567: exprs.append("7/9 * 729 = A_max")
        if v == 189: exprs.append("7 * 27 = D3 * n^3")
        if v == 297: exprs.append("11 * 27")
        if v == 351: exprs.append("13 * 27 = |PG| * n^3")
        if v == 513: exprs.append("19 * 27")
        if v == 489: exprs.append("163/243 * 729")

        print("  %5d: %s" % (v, "; ".join(exprs) if exprs else "---"))

    # ================================================================
    print()
    print("=" * 70)
    print("8. THE SIX RULE-RANGES AS BANDS")
    print("=" * 70)
    print()

    print("  Rule spectrum ranges (non-overlapping?):")
    bands = []
    for g in CROSS_RULES:
        vals = sorted(set(m["assoc"] for m in all_magmas if m["rule"] == g))
        bands.append((vals[0], vals[-1], g))
    bands.sort()
    for lo, hi, g in bands:
        bar = " " * ((lo - 162) // 8) + "#" * max(1, (hi - lo) // 8)
        print("  %s: [%3d, %3d]  %s" % (CROSS_NAMES[g], lo, hi, bar))

    # Check overlaps between ranges
    print()
    for i in range(len(bands)):
        for j in range(i+1, len(bands)):
            lo1, hi1, g1 = bands[i]
            lo2, hi2, g2 = bands[j]
            if lo1 <= hi2 and lo2 <= hi1:
                overlap = min(hi1,hi2) - max(lo1,lo2)
                print("  Overlap: %s [%d,%d] & %s [%d,%d] = %d" %
                      (CROSS_NAMES[g1], lo1, hi1, CROSS_NAMES[g2], lo2, hi2, overlap))

    print()
    print("=" * 70)
    print("SPECTRUM ANALYSIS COMPLETE")
    print("=" * 70)

if __name__ == "__main__":
    main()
