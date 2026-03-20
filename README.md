# TAB-PAB-CTA
A Nine-Element Non-Associative Magma from Type Decoherence of su(2): Algebra, Classification, and Constrained Hamiltonian Structure

# PAB Magma — Computational Verification Suite

This repository contains a single Python script, `pab_verify.py`, which exhaustively verifies the finite computational claims stated in Appendix C of the associated paper.

The script covers:

- **C.1** Cayley table from axioms AX0–AX3
- **C.2** associativity count and related algebraic properties
- **C.3** automorphism group
- **C.4** classification of 162 sigma-invariant magmas
- **C.5** isomorphism classes via Burnside counting
- **C.6** Hamiltonian path enumeration
- **C.7** decoherence verification

## Repository contents

- `pab_verify.py` — standalone verification script

## Requirements

- Python 3
- No external dependencies

## Run

`python pab_verify.py`

## Expected result

A successful run ends with:
ALL VERIFICATIONS PASSED

The script prints intermediate results for each verification block, including the Cayley table, associativity statistics, automorphism count, classification summaries, Burnside counts, Hamiltonian paths, and decoherence checks.

## Scope
This code is intended as a computational verification artifact for a finite algebraic construction. It is not a general-purpose package and does not require installation.

## Associated paper
The accompanying paper is published separately. This repository contains the code artifact used to verify the corresponding Appendix C claims.
