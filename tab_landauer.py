#!/usr/bin/env python3
"""
tab_landauer.py — Information-theoretic analysis of TAB cross-rules
====================================================================
Computes:
  1. H_storage = I(C_out; C1, C2)          — unconditional MI
  2. H_access  = I(C_out; C1, C2 | R1, R2) — conditional MI
  3. Correlation with associativity across all 162 magmas
  4. Summary: three classifications of cross-rules

Run: python tab_landauer.py
"""

import sys, io, math
from itertools import product
from collections import Counter

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")

S = ["P", "A", "B"]
M = [(R, C) for R in S for C in S]


def nu(x):
    return S[(S.index(x) + 1) % 3]


def compl(a, b):
    return [x for x in S if x != a and x != b][0]


def cross_out(g, R1, C1, R2, C2):
    """Output column for cross-row operation (R1 != R2)."""
    if g == "g1": return R2
    if g == "g2": return C2
    if g == "g3": return compl(C1, C2) if C1 != C2 else C1
    if g == "g4": return C1
    if g == "g5": return R1
    if g == "g6": return compl(R2, C2) if R2 != C2 else R2


def H(counts, total):
    """Shannon entropy in bits."""
    return -sum(c/total * math.log2(c/total)
                for c in counts.values() if c > 0)


# ---- Constants ----
RULES = ["g1", "g2", "g3", "g4", "g5", "g6"]
NAMES = {"g1": "column-blind", "g2": "transparent", "g3": "su(2)-transp",
         "g4": "echo", "g5": "diagonal", "g6": "anti-compl"}
cross_rows = [(R1, R2) for R1 in S for R2 in S if R1 != R2]  # 6
col_pairs  = [(C1, C2) for C1 in S for C2 in S]              # 9
N_cross = 54  # 6 x 9

print("=" * 72)
print("  LANDAUER ANALYSIS: Information Content of TAB Cross-Rules")
print("=" * 72)

# ====================================================================
# 1. Unconditional MI: I(C_out; C1, C2)
# ====================================================================
print("\n--- 1. Unconditional MI: I(C_out; C1, C2) ---")
print("    How much column info does the OUTPUT carry?\n")

mi_uncond = {}
for g in RULES:
    out_cnt = Counter()
    joint_cnt = Counter()  # (C1, C2, C_out)
    for R1, R2 in cross_rows:
        for C1, C2 in col_pairs:
            co = cross_out(g, R1, C1, R2, C2)
            out_cnt[co] += 1
            joint_cnt[(C1, C2, co)] += 1

    H_out = H(out_cnt, N_cross)
    H_cols = math.log2(9)        # uniform over 9 column pairs
    H_joint = H(joint_cnt, N_cross)
    mi = H_out + H_cols - H_joint
    mi_uncond[g] = mi

    print(f"  {g} ({NAMES[g]:14s}):  I = {mi:.4f} bits   "
          f"H(out)={H_out:.3f}  H(C1,C2)={H_cols:.3f}  H(joint)={H_joint:.3f}")

print(f"\n  Reference: log2(3) = {math.log2(3):.4f} bits")

# ====================================================================
# 2. Conditional MI: I(C_out; C1, C2 | R1, R2)
# ====================================================================
print("\n--- 2. Conditional MI: I(C_out; C1, C2 | R1, R2) ---")
print("    Given rows, how much column info is USED in computation?\n")

mi_cond = {}
for g in RULES:
    # Average H(C_out | R1, R2) over all row-pairs
    sum_H_given_rows = 0.0
    for R1, R2 in cross_rows:
        out_cnt = Counter()
        for C1, C2 in col_pairs:
            out_cnt[cross_out(g, R1, C1, R2, C2)] += 1
        sum_H_given_rows += H(out_cnt, 9)
    H_cond_rows = sum_H_given_rows / 6

    # H(C_out | R1,R2,C1,C2) = 0 (deterministic)
    cmi = H_cond_rows
    mi_cond[g] = cmi

    print(f"  {g} ({NAMES[g]:14s}):  I_cond = {cmi:.4f} bits   "
          f"H(out|rows) = {H_cond_rows:.3f}")

# ====================================================================
# 3. Detail: output table for each rule
# ====================================================================
print("\n--- 3. Cross-row output tables ---\n")
for g in RULES:
    print(f"  {g} ({NAMES[g]}): (R1,C1)*(R2,C2) -> column = ...")
    for R1, R2 in cross_rows[:2]:  # show 2 row-pairs as examples
        row = []
        for C1, C2 in col_pairs:
            row.append(cross_out(g, R1, C1, R2, C2))
        cols_str = " ".join(f"{c}" for c in row)
        print(f"    ({R1},{R2}): {cols_str}")
    print()

# ====================================================================
# 4. Associativity correlation
# ====================================================================
print("--- 4. Associativity by rule at PAB delta (P,P,P) ---\n")


def sigma_diags():
    diags = []
    for dp, da, db in product(S, S, S):
        d = {}
        d[("P","P")] = dp;       d[("P","A")] = da;       d[("P","B")] = db
        d[("A","A")] = nu(dp);   d[("A","B")] = nu(da);   d[("A","P")] = nu(db)
        d[("B","B")] = nu(nu(dp)); d[("B","P")] = nu(nu(da)); d[("B","A")] = nu(nu(db))
        diags.append(((dp, da, db), d))
    return diags


def make_meet(g, dd):
    def meet(x, y):
        R1, C1 = x; R2, C2 = y
        if R1 == R2:
            return (R1, dd[(R1, C1)]) if C1 == C2 else (R1, compl(C1, C2))
        return (R1, cross_out(g, R1, C1, R2, C2))
    return meet


def assoc_count(mf):
    n = 0
    for a in M:
        for b in M:
            ab = mf(a, b)
            for c in M:
                if mf(ab, c) == mf(a, mf(b, c)):
                    n += 1
    return n


# PAB delta: d(R,C) = R for all
pab_d = {(R, C): R for R in S for C in S}

for g in RULES:
    a = assoc_count(make_meet(g, pab_d))
    print(f"  {g} ({NAMES[g]:14s}): {a}/729 = {a/729:.4f}  ({100*a/729:.1f}%)")

# ====================================================================
# 5. Full 162-magma analysis: Assoc ranges by H_access group
# ====================================================================
print("\n--- 5. Assoc ranges by H_access group (all 162 magmas) ---\n")

all_diags = sigma_diags()
by_rule = {g: [] for g in RULES}

for g in RULES:
    for (dp, da, db), dd in all_diags:
        a = assoc_count(make_meet(g, dd))
        by_rule[g].append(a)

# Group by H_access
group0 = []  # H_access = 0: g1, g5
group1 = []  # H_access = log2(3): g2, g3, g4, g6
for g in RULES:
    if mi_cond[g] < 0.01:
        group0.extend(by_rule[g])
    else:
        group1.extend(by_rule[g])

print(f"  H_access = 0     (g1, g5):       "
      f"n={len(group0):3d}  range=[{min(group0)}, {max(group0)}]  "
      f"mean={sum(group0)/len(group0):.1f}")
print(f"  H_access = log2(3) (g2,g3,g4,g6): "
      f"n={len(group1):3d}  range=[{min(group1)}, {max(group1)}]  "
      f"mean={sum(group1)/len(group1):.1f}")

print("\n  Per-rule breakdown:")
for g in RULES:
    v = by_rule[g]
    print(f"    {g} ({NAMES[g]:14s}): [{min(v):3d}, {max(v):3d}]  "
          f"mean={sum(v)/len(v):.1f}  H_acc={'0':>7s}" if mi_cond[g] < 0.01
          else f"    {g} ({NAMES[g]:14s}): [{min(v):3d}, {max(v):3d}]  "
               f"mean={sum(v)/len(v):.1f}  H_acc={'log2(3)':>7s}")

# ====================================================================
# 6. Summary table
# ====================================================================
print("\n" + "=" * 72)
print("  SUMMARY: Three Classifications of Cross-Rules")
print("=" * 72)
print()
print("  Rule            H_storage  H_access  C2-indep  SRR=162  DIST=0")
print("  --------------- ---------  --------  --------  -------  ------")

for g in RULES:
    hs = "0" if mi_uncond[g] < 0.01 else "log2(3)"
    ha = "0" if mi_cond[g] < 0.01 else "log2(3)"
    c2 = "yes" if g in ("g1", "g4", "g5") else "no"
    srr = "yes" if g in ("g1", "g4", "g5") else "no"
    dist0 = "yes" if g == "g1" else "no"
    print(f"  {g} {NAMES[g]:14s} {hs:>7s}    {ha:>7s}   {c2:>5s}     "
          f"{srr:>4s}     {dist0:>3s}")

print()
print("  Classification 1 — H_storage = 0:    {g1, g5, g6}")
print("     Output carries no column info (erased or mixed away)")
print()
print("  Classification 2 — H_access = 0:     {g1, g5}")
print("     Computation USES no column info (pure row-determined)")
print()
print("  Classification 3 — C2-independent:   {g1, g4, g5}")
print("     Buffer theorem: SRR = 162")
print()
print("  LANDAUER SELECTION:")
print("    Step 1: Minimum H_access -> {g1, g5}")
print("    Step 2: Causal direction (output = R2, not R1) -> g1")
print("    Result: PAB uniquely selected.")
print()
print("  NOTE: H_access does NOT correlate with associativity.")
print("  g5 (H_access=0) has HIGH assoc; g4 (H_access=log2(3)) has HIGHEST.")
print("  The Landauer argument selects g1 for INFORMATION reasons,")
print("  not associativity reasons. Both paths converge to g1.")
