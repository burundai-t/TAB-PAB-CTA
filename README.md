# TAB-PAB-CTA
σ-EQUIVARIANT MAGMAS ON AG(2, 3): CLASSIFICATION AND INFORMATION-THEORETIC SELECTION

## Repository contents

- `pab verify.py` — Cayley table from AX0–AX3; Assoc = 219/729; |Aut| = 6; 3 idempotents; identities I1–I3; STS fiber; 3 Hamiltonian paths (out of 32); Tβ-count = 2 pair

- `tab classify.py` — All 162 magmas; Burnside = 90 iso classes; 18 canonical; 6 standard; uniqueness of PAB

- `tab patterns.py` — T1–T4, T6, T7, T8, T9 for all 162 magmas PASS

- `tab formulas.py` — Master formula (Theorem 4.9): five-parameter prediction matches exhaustive count for all 162 magmas

- `tab spectrum.py` — Exactly 28 distinct values; all divisible by 3; anchor points; gap structure

- `tab graph.py` — Mutation graph: 162 vertices, 891 edges, 11-regular, connected (BFS from every vertex)

- `tab landauer.py` — Hstorage, Haccess for all 6 rules over all 54 inputs; binary classification (Theorem 5.3);g6: Hstorage = 0, Haccess = log2 3

- `tab landauer step2.py` — Intersection theorem (Theorem 5.6): min(Haccess) ∩ min(λ) = {g1}; Fβ(g1) < Fβ(g) for all g ̸= g1 and β ∈ {0.01, 0.1, 1, 5, 100} (Corollary 5.7)

## Requirements

- Python 3
- No external dependencies

## Run

`python pab verify.py`

`python tab classify.py`

`python tab patterns.py`

`python tab formulas.py`

`python tab spectrum.py`

`python tab graph.py`

`python tab landauer.py`

`python tab landauer step2.py`

## Expected result

A successful run ends with:
ALL VERIFICATIONS PASSED / PASS

## Scope
This code is intended as a computational verification artifact for a finite algebraic construction. It is not a general-purpose package and does not require installation.

## Associated paper
The accompanying paper is published separately. This repository contains the code artifact used to verify the corresponding Appendix C claims.
https://doi.org/10.5281/zenodo.19163946
