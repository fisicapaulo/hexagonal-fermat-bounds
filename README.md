# Hexagonal Fermat Bounds — Exact Verification & Benchmarks

Companion code for the paper:

> **Paulo Vieira**, *"Bounds for Fermat quotients via hexagonal interpolation determinants"* (2026).

## Overview

Although the paper is completely self-contained and all numerical comparisons in Section 4.2 can be checked by hand, this minimal repository provides a dependency-free script (`verify_and_demo.py`) using **only Python's standard library** (`fractions.Fraction` and exact integers) to:

1. **Certify all rational inequalities in Section 4** (`certify_section_4`) using exact rational series bounds (zero floating-point arithmetic), confirming the total normalized costs `1851/1000 < 119/64` (general and quadratic coprime case, $C = 15$) and `1857/1000 < 119/64` (base-2 and quadratic case, $C = 13$).
2. **Benchmark explicit upper bounds for $v_p(2^{p-1}-1)$** (`benchmark_bounds`), comparing the elementary bound, Yamada (2010, constant $283$), Niedbala Giraudin (2026, constant $69$), and Theorem 1.1 ($C = 13$ and quadratic refinement $13/2$).
3. **Illustrate the selective weight $\mathcal{B}(p_j, p_i)$ for Odd Perfect Numbers** (`benchmark_opn_edge`) in Theorem 5.1, highlighting the automatic quadratic-residue halving ($\kappa_{i,j} = 15/2$) induced by Euler's form when $j \neq \star$ and $p_i \equiv \pm 1 \pmod 8$.

## Requirements & Usage

Requires **Python 3.9+** (no external packages needed). Run:

```bash
python3 verify_and_demo.py
```

## License

Released under the MIT License.
