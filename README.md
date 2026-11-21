# Hexagonal Fermat Bounds — Exact Rational Verification & Benchmarks

[![DOI](https://zenodo.org/badge/1466177217.svg)](https://doi.org/10.5281/zenodo.24412370)

Companion code for the paper:

> **Paulo Vieira**, *"Bounds for Fermat quotients via hexagonal interpolation determinants"* (2026).

## Overview

Although the paper is completely self-contained and all numerical inequalities in Section 4 can be checked by hand, this repository provides a single dependency-free script (`verify_and_demo.py`) written in pure Python (`fractions.Fraction` and exact integers only; **no floating-point arithmetic is used anywhere**, not even when rounding table entries).

Specifically, the script performs three tasks:

1. **Exact rational certification of Section 4 (`certify_section_4`):**
   Encloses all logarithms via the rational series bounds from Section 4.2, verifies the monotonicity hypotheses (including $\log(0.9 \cdot 15) > 1$ and $\log(13) > 1$), checks the square-root and Stirling/Taylor bounds, and confirms the normalized inequalities
   - `1851/1000 < 119/64` for the general coprime case and its quadratic refinement ($C = 15$);
   - `1857/1000 < 119/64` for the base-$2$ specialization and its quadratic refinement ($C = 13$).
   All checks use an explicit `require(...)` guard rather than `assert`, so they remain active even under `python3 -O`.

2. **Base-$2$ upper bounds with certified integer rounding (`benchmark_bounds`):**
   Finds deterministic primes $p \equiv \pm 1 \pmod{24}$ (checked via deterministic Miller–Rabin below $3.18 \times 10^{23}$ and explicit Euler-criterion tests for bases $2$ and $3$) and compares upper bounds for $v_p(2^{p-1} - 1)$ across:
   - the elementary bound $\lceil (p-1)\log 2 / \log p \rceil - 1$;
   - Yamada (2010, constant $283$);
   - Niedbala Giraudin (2026, constant $69$);
   - Theorem 1.1 ($C = 13$ and quadratic refinement $13/2$).
   Every `floor` and `ceil` value is certified via rational interval enclosure (`rounded_log_ratio`), returning an integer only when the lower and upper rational endpoints round to the exact same integer.

3. **Illustrative admissible weight $\mathcal{B}(p_j, p_i)$ for Section 5 (`benchmark_opn_edge`):**
   Computes and certifies the **exact multiplicative order** $f_{i,j} = \operatorname{ord}_{p_i}(p_j) = 500{,}019$ for the admissible non-Eulerian pair $(p_j, p_i) = (10{,}000{,}019,\; 1{,}000{,}039)$ with $p_i \equiv 7 \pmod 8$, verifies minimality of the order across all prime divisors of $p_i - 1$, and compares the resulting weight $\mathcal{B}(p_j, p_i) = \min\{\text{elementary},\, \text{quadratic hexagonal}\}$ against unrefined bounds. As noted in the output, this illustrates the local weight on an admissible pair and makes no claim about the existence of an odd perfect number.

## Requirements & Usage

Requires **Python 3.9+** (standard library only; zero external dependencies). Run:

```bash
python3 verify_and_demo.py
```

## Sample Output

```text
1. EXACT RATIONAL CHECKS FOR SECTION 4
   [OK] General and quadratic cases: 1851/1000 < 119/64
   [OK] Base 2 and quadratic cases:  1857/1000 < 119/64

2. BASE-2 UPPER BOUNDS (CERTIFIED INTEGER ROUNDING)
          Prime p |     Elementary |         Yamada |    N. Giraudin |        Hex. 13 |      Hex. 13/2
------------------------------------------------------------------------------------------------------
        1,594,559 |         77,387 |      4,354,784 |      1,061,770 |         77,389 |         38,695
       10,000,079 |        430,045 |     21,443,031 |      5,228,162 |        381,057 |        190,529
    1,000,000,007 |     33,447,777 |  1,297,161,958 |    316,269,173 |     23,051,378 |     11,525,690
1,000,000,000,177 | 25,085,832,976 | 729,653,595,729 | 177,901,406,735 | 12,966,399,449 |  6,483,199,725

3. ILLUSTRATIVE WEIGHT FOR THE ODD PERFECT NUMBER APPLICATION
   Source -> target: 10,000,019 -> 1,000,039; actual order: 500,019
   This is an admissible pair, not a claim about an existing odd perfect number.
   Elementary bound using the actual order           583,353
   Yamada (283)                                   17,277,943
   Niedbala Giraudin (69)                          4,212,646
   Hexagonal bound (15)                              878,034
   Hexagonal quadratic bound (15/2)                  439,018
   Weight B = min(elementary, quadratic)             439,018
```

## License

Released under the MIT License.
