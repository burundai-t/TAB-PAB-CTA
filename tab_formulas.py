#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TAB Formulas: Resolve remaining gaps
=====================================
1. RRR as function of delta
2. SRR for g2, g6
3. The g1/g5 split at (fix=3, st=3)
4. Full analytic formula attempt

Run:  python tab_formulas.py
"""

import sys
import io

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

from itertools import product
from collections import defaultdict

S = ["P", "A", "B"]
M = [(R, C) for R in S for C in S]


def nu(x):
    return S[(S.index(x) + 1) % 3]


def complement(a, b):
    return [x for x in S if x != a and x != b][0]


CROSS_RULES = ["g1", "g2", "g3", "g4", "g5", "g6"]
CROSS_NAMES = {
    "g1": "col-blind", "g2": "transp", "g3": "su2-tr",
    "g4": "echo", "g5": "diag", "g6": "anti-c"
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


def compute_all_patterns(table):
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
    return pats


def delta_full_props(d_dict):
    """Extended delta properties."""
    p = {}
    # fix: delta(R,C) = C
    p["fix"] = sum(1 for R in S for C in S if d_dict[(R, C)] == C)
    # self_type: delta(R,C) = R
    p["self_type"] = sum(1 for R in S for C in S if d_dict[(R, C)] == R)
    # diag_self: delta(R,R) = R
    p["diag_self"] = sum(1 for R in S if d_dict[(R, R)] == R)
    # constant per row
    p["const"] = all(len(set(d_dict[(R, C)] for C in S)) == 1 for R in S)
    # P-row distinct count
    p["p_dist"] = len(set(d_dict[("P", C)] for C in S))
    # specific P-row values
    p["dPP"] = d_dict[("P", "P")]
    p["dPA"] = d_dict[("P", "A")]
    p["dPB"] = d_dict[("P", "B")]
    return p


def main():
    print()
    print("TAB Formulas: Resolving Gaps")
    print("=" * 70)

    diags = sigma_equivariant_diagonals()

    # Build all data
    all_data = []
    for g in CROSS_RULES:
        for (dp, da, db), d_dict in diags:
            meet = make_meet(g, d_dict)
            table = build_cayley(meet)
            pats = compute_all_patterns(table)
            props = delta_full_props(d_dict)
            all_data.append({
                "rule": g, "delta": (dp, da, db),
                "pats": pats, "props": props
            })

    # ================================================================
    print()
    print("=" * 70)
    print("1. RRR AS FUNCTION OF DELTA")
    print("=" * 70)
    print()

    # RRR depends only on delta (proven). What property of delta?
    # Candidates: fix, self_type, diag_self, const, p_dist

    # Group by delta (cross-rule independent)
    rrr_by_delta = {}
    for d in all_data:
        key = d["delta"]
        if key not in rrr_by_delta:
            rrr_by_delta[key] = d["pats"]["RRR"]

    # Test: RRR = f(fix)?
    fix_rrr = defaultdict(set)
    for d in all_data:
        if d["rule"] == "g1":
            fix_rrr[d["props"]["fix"]].add(d["pats"]["RRR"])

    print("  RRR vs fix:")
    for f in sorted(fix_rrr.keys()):
        print("    fix=%d: RRR=%s" % (f, sorted(fix_rrr[f])))
    print("  RRR = f(fix) alone? %s" %
          all(len(v) == 1 for v in fix_rrr.values()))

    # Test: RRR = f(self_type)?
    st_rrr = defaultdict(set)
    for d in all_data:
        if d["rule"] == "g1":
            st_rrr[d["props"]["self_type"]].add(d["pats"]["RRR"])

    print()
    print("  RRR vs self_type:")
    for st in sorted(st_rrr.keys()):
        print("    self_type=%d: RRR=%s" % (st, sorted(st_rrr[st])))
    print("  RRR = f(self_type)? %s" %
          all(len(v) == 1 for v in st_rrr.values()))

    # Test: RRR = f(diag_self)?
    ds_rrr = defaultdict(set)
    for d in all_data:
        if d["rule"] == "g1":
            ds_rrr[d["props"]["diag_self"]].add(d["pats"]["RRR"])

    print()
    print("  RRR vs diag_self:")
    for ds in sorted(ds_rrr.keys()):
        print("    diag_self=%d: RRR=%s" % (ds, sorted(ds_rrr[ds])))

    # Test: RRR = f(const)?
    c_rrr = defaultdict(set)
    for d in all_data:
        if d["rule"] == "g1":
            c_rrr[d["props"]["const"]].add(d["pats"]["RRR"])

    print()
    print("  RRR vs const:")
    for c in sorted(c_rrr.keys()):
        print("    const=%s: RRR=%s" % (c, sorted(c_rrr[c])))

    # Test: RRR = f(p_dist)?
    pd_rrr = defaultdict(set)
    for d in all_data:
        if d["rule"] == "g1":
            pd_rrr[d["props"]["p_dist"]].add(d["pats"]["RRR"])

    print()
    print("  RRR vs p_dist (P-row distinct values):")
    for pd in sorted(pd_rrr.keys()):
        print("    p_dist=%d: RRR=%s" % (pd, sorted(pd_rrr[pd])))

    # Combined: RRR = f(self_type, diag_self)?
    print()
    combo_rrr = defaultdict(set)
    for d in all_data:
        if d["rule"] == "g1":
            key = (d["props"]["self_type"], d["props"]["diag_self"])
            combo_rrr[key].add(d["pats"]["RRR"])

    print("  RRR vs (self_type, diag_self):")
    all_unique = True
    for key in sorted(combo_rrr.keys()):
        vals = sorted(combo_rrr[key])
        u = len(vals) == 1
        if not u:
            all_unique = False
        print("    (st=%d, ds=%d): RRR=%s %s" %
              (key[0], key[1], vals, "" if u else "<-- AMBIGUOUS"))
    print("  RRR = f(self_type, diag_self)? %s" % all_unique)

    # ================================================================
    print()
    print("=" * 70)
    print("2. SRR FOR g2 AND g6")
    print("=" * 70)
    print()

    for g in ["g2", "g6"]:
        print("  --- %s ---" % CROSS_NAMES[g])
        fix_srr = defaultdict(set)
        st_srr = defaultdict(set)
        ds_srr = defaultdict(set)
        for d in all_data:
            if d["rule"] != g:
                continue
            fix_srr[d["props"]["fix"]].add(d["pats"]["SRR"])
            st_srr[d["props"]["self_type"]].add(d["pats"]["SRR"])
            ds_srr[d["props"]["diag_self"]].add(d["pats"]["SRR"])

        print("  SRR vs fix:")
        for f in sorted(fix_srr.keys()):
            print("    fix=%d: SRR=%s" % (f, sorted(fix_srr[f])))
        srr_fix = all(len(v) == 1 for v in fix_srr.values())
        print("  SRR = f(fix)? %s" % srr_fix)

        if srr_fix:
            print("  Formula: SRR = ", end="")
            vals = [(f, list(fix_srr[f])[0]) for f in sorted(fix_srr.keys())]
            # Try linear fit
            if len(vals) >= 2:
                x0, y0 = vals[0]
                x1, y1 = vals[-1]
                if x1 != x0:
                    slope = (y1 - y0) / (x1 - x0)
                    intercept = y0 - slope * x0
                    # Verify
                    ok = all(abs(intercept + slope * x - y) < 0.01 for x, y in vals)
                    if ok:
                        if intercept == int(intercept) and slope == int(slope):
                            print("%d + %d * fix" % (int(intercept), int(slope)))
                        else:
                            print("%.1f + %.1f * fix" % (intercept, slope))
                    else:
                        print("(not linear)")
        print()

    # ================================================================
    print()
    print("=" * 70)
    print("3. THE g1/g5 SPLIT: diag_self AS THIRD PARAMETER")
    print("=" * 70)
    print()

    for g in ["g1", "g5"]:
        print("  --- %s ---" % CROSS_NAMES[g])
        combo3 = defaultdict(lambda: {"rrs": set(), "rsr": set(), "srr": set()})
        for d in all_data:
            if d["rule"] != g:
                continue
            key = (d["props"]["fix"], d["props"]["self_type"], d["props"]["diag_self"])
            combo3[key]["rrs"].add(d["pats"]["RRS"])
            combo3[key]["rsr"].add(d["pats"]["RSR"])
            combo3[key]["srr"].add(d["pats"]["SRR"])

        all_resolved = True
        print("  (fix, st, ds) -> (RRS, RSR, SRR)")
        for key in sorted(combo3.keys()):
            c = combo3[key]
            rrs_u = len(c["rrs"]) == 1
            rsr_u = len(c["rsr"]) == 1
            srr_u = len(c["srr"]) == 1
            if not (rrs_u and rsr_u and srr_u):
                all_resolved = False
            rrs_s = list(c["rrs"])[0] if rrs_u else sorted(c["rrs"])
            rsr_s = list(c["rsr"])[0] if rsr_u else sorted(c["rsr"])
            srr_s = list(c["srr"])[0] if srr_u else sorted(c["srr"])
            flag = "" if (rrs_u and rsr_u and srr_u) else " <-- AMBIGUOUS"
            print("    (%d,%d,%d) -> RRS=%-6s RSR=%-6s SRR=%-6s%s" %
                  (key[0], key[1], key[2], rrs_s, rsr_s, srr_s, flag))

        print("  Fully resolved with (fix, st, ds)? %s" % all_resolved)
        print()

    # ================================================================
    print()
    print("=" * 70)
    print("4. FULL ANALYTIC FORMULAS")
    print("=" * 70)
    print()

    # For each rule, attempt: Assoc = f(RRR, fix, self_type, diag_self)
    for g in CROSS_RULES:
        print("  --- %s ---" % CROSS_NAMES[g])
        combo_full = defaultdict(lambda: {"assoc": set(), "rrr": set(),
                                          "rrs": set(), "rsr": set(),
                                          "srr": set(), "dist": set()})
        for d in all_data:
            if d["rule"] != g:
                continue
            key = (d["props"]["fix"], d["props"]["diag_self"])
            for pat in ["RRR", "RRS", "RSR", "SRR", "DIST"]:
                combo_full[key][pat.lower()].add(d["pats"][pat])
            combo_full[key]["assoc"].add(d["pats"]["RRR"] + d["pats"]["RRS"] +
                                         d["pats"]["RSR"] + d["pats"]["SRR"] +
                                         d["pats"]["DIST"])

        resolved = True
        for key in sorted(combo_full.keys()):
            c = combo_full[key]
            for pat in ["rrr", "rrs", "rsr", "srr", "dist"]:
                if len(c[pat]) > 1:
                    resolved = False
                    break

        if resolved:
            print("  RESOLVED with (fix, diag_self)!")
            print("  (fix, ds) -> RRR  RRS  RSR  SRR  DIST  Total")
            for key in sorted(combo_full.keys()):
                c = combo_full[key]
                vals = {pat: list(c[pat])[0] for pat in ["rrr", "rrs", "rsr", "srr", "dist"]}
                total = sum(vals.values())
                print("    (%d,%d) -> %4d %4d %4d %4d %4d  %5d" %
                      (key[0], key[1], vals["rrr"], vals["rrs"], vals["rsr"],
                       vals["srr"], vals["dist"], total))
        else:
            # Try (fix, self_type, diag_self)
            combo_full3 = defaultdict(lambda: {"assoc": set()})
            for d in all_data:
                if d["rule"] != g:
                    continue
                key3 = (d["props"]["fix"], d["props"]["self_type"], d["props"]["diag_self"])
                total = sum(d["pats"][p] for p in ["RRR", "RRS", "RSR", "SRR", "DIST"])
                combo_full3[key3]["assoc"].add(total)

            resolved3 = all(len(v["assoc"]) == 1 for v in combo_full3.values())
            print("  NOT resolved with (fix, ds). Try (fix, st, ds): %s" % resolved3)
            if resolved3:
                print("  (fix, st, ds) -> Total")
                for key3 in sorted(combo_full3.keys()):
                    print("    (%d,%d,%d) -> %d" %
                          (key3[0], key3[1], key3[2], list(combo_full3[key3]["assoc"])[0]))
        print()

    # ================================================================
    print()
    print("=" * 70)
    print("5. VERIFICATION: PREDICT ALL 162 ASSOCIATIVITIES")
    print("=" * 70)
    print()

    # Build prediction tables for rules that resolved with (fix, ds)
    for g in CROSS_RULES:
        lookup = {}
        all_resolved_2 = True
        for d in all_data:
            if d["rule"] != g:
                continue
            key = (d["props"]["fix"], d["props"]["diag_self"])
            total = sum(d["pats"][p] for p in ["RRR", "RRS", "RSR", "SRR", "DIST"])
            if key in lookup:
                if lookup[key] != total:
                    all_resolved_2 = False
            else:
                lookup[key] = total

        if all_resolved_2:
            # Verify all 27 magmas
            errors = 0
            for d in all_data:
                if d["rule"] != g:
                    continue
                key = (d["props"]["fix"], d["props"]["diag_self"])
                predicted = lookup[key]
                actual = sum(d["pats"][p] for p in ["RRR", "RRS", "RSR", "SRR", "DIST"])
                if predicted != actual:
                    errors += 1
            print("  %s: Assoc = f(fix, diag_self) -- %d errors / 27" %
                  (CROSS_NAMES[g], errors))
        else:
            # Try 3-param
            lookup3 = {}
            for d in all_data:
                if d["rule"] != g:
                    continue
                key3 = (d["props"]["fix"], d["props"]["self_type"], d["props"]["diag_self"])
                total = sum(d["pats"][p] for p in ["RRR", "RRS", "RSR", "SRR", "DIST"])
                lookup3[key3] = total

            errors3 = 0
            for d in all_data:
                if d["rule"] != g:
                    continue
                key3 = (d["props"]["fix"], d["props"]["self_type"], d["props"]["diag_self"])
                if lookup3[key3] != sum(d["pats"][p] for p in ["RRR", "RRS", "RSR", "SRR", "DIST"]):
                    errors3 += 1
            print("  %s: Assoc = f(fix, st, ds) -- %d errors / 27" %
                  (CROSS_NAMES[g], errors3))

    total_magmas = len(all_data)
    print()
    print("  Total: all 162 magmas covered by at most 3 delta-parameters")

    print()
    print("=" * 70)
    print("FORMULAS COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
