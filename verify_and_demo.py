#!/usr/bin/env python3
"""
Companion script for:
"Bounds for Fermat quotients via hexagonal interpolation determinants"
Author: Paulo Vieira

Pure Python 3 (standard library only, exact rational arithmetic via fractions.Fraction).
Performs:
  1. Exact rational certification of all numerical inequalities in Section 4.
  2. Explicit comparison table: Elementary vs. Yamada (2010) vs. Giraudin (2026) vs. Theorem 1.1.
  3. Selective weight B(p_j, p_i) comparison for Odd Perfect Numbers (Theorem 5.1).
"""

from fractions import Fraction
import math


# ==============================================================================
# Deterministic Primality Check (Miller-Rabin exact for all n < 3.3 * 10^24)
# ==============================================================================

def is_prime(n: int) -> bool:
    if n < 2:
        return False
    small_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True


def next_prime_pm1_mod24(start: int) -> int:
    """Return the smallest prime p >= start such that p == +-1 (mod 24)."""
    p = start if start % 2 == 1 else start + 1
    while True:
        if p % 24 in (1, 23) and is_prime(p):
            return p
        p += 2


# ==============================================================================
# PART 1: Exact Rational Certification of Section 4 (Zero Floating-Point!)
# ==============================================================================

def log_bounds_series(z: Fraction, m: int = 12) -> tuple[Fraction, Fraction]:
    """
    Exact rational lower and upper bounds for log(z) (z >= 1) using Section 4.2:
    z = 2^h * r with 1 <= r < 2, and t = (w - 1)/(w + 1) <= 1/3.
    """
    assert z >= 1
    h = 0
    r = z
    while r >= 2:
        r /= 2
        h += 1

    def unit_interval_bounds(w: Fraction) -> tuple[Fraction, Fraction]:
        if w == 1:
            return Fraction(0), Fraction(0)
        t = (w - 1) / (w + 1)
        t2 = t * t
        power = t
        partial = Fraction(0)
        for j in range(m):
            partial += power / (2 * j + 1)
            power *= t2
        low = 2 * partial
        rem = (2 * power) / ((2 * m + 1) * (1 - t2))
        return low, low + rem

    low2, high2 = unit_interval_bounds(Fraction(2))
    low_r, high_r = unit_interval_bounds(r)
    return h * low2 + low_r, h * high2 + high_r


def exp_3_halves_upper() -> Fraction:
    """Exact rational upper bound for e^(3/2) / 2 from Section 4.2 (degree 12 + tail)."""
    x = Fraction(3, 2)
    fact = 1
    power = Fraction(1)
    total = Fraction(0)
    for k in range(13):
        if k > 0:
            fact *= k
            power *= x
        total += power / fact
    fact13 = fact * 13
    tail = (power * x) / (fact13 * (1 - Fraction(3, 28)))
    return (total + tail) / 2


def certify_section_4() -> None:
    print("=" * 88)
    print("1. EXACT RATIONAL CERTIFICATION OF SECTION 4 INEQUALITIES")
    print("=" * 88)

    log2_lo, log2_hi = log_bounds_series(Fraction(2), m=12)
    log3_lo, log3_hi = log_bounds_series(Fraction(3), m=12)
    assert Fraction(693147, 1000000) < log2_lo and log2_hi < Fraction(693148, 1000000)
    assert Fraction(1098612, 1000000) < log3_lo and log3_hi < Fraction(1098613, 1000000)
    assert log2_lo * log3_lo > Fraction(76, 100)
    assert 13 * log2_hi * log3_hi < Fraction(99, 10)

    e32_half_hi = exp_3_halves_upper()
    assert e32_half_hi < Fraction(2241, 1000)

    # Case C = 15 (including quadratic refinement g >= (p-1)/2 for p > 2^15 => p >= 32771)
    log_p_gen_lo, log_p_gen_hi = log_bounds_series(Fraction(32771), m=12)
    assert Fraction(1039, 100) < log_p_gen_lo and log_p_gen_hi < Fraction(104, 10)
    A_g_gen_lo = Fraction(16384) * (log2_lo * log3_lo) / (log_p_gen_hi ** 2)
    assert A_g_gen_lo > 113
    s_gen = (20 * 113) // 3
    assert s_gen >= 750

    paren1 = Fraction(1000001, 1000000) + Fraction(2982, 100000)
    assert Fraction(1) + Fraction(1, 3 * 750 * 16384) < Fraction(1000001, 1000000) ** 2
    assert Fraction(2 + Fraction(1, 16384), 3 * 750) < Fraction(2982, 100000) ** 2
    lhs_sq_gen = Fraction(36, 15) * Fraction(751, 750) * (paren1 ** 2)
    assert lhs_sq_gen < Fraction(1598, 1000) ** 2

    log_1002_lo, _ = log_bounds_series(Fraction(1002, 1000), m=12)
    assert (2 * log2_hi) / 751 < log_1002_lo
    b0_over_w_gen = Fraction(1002, 1000) * Fraction(2241, 1000) * Fraction(1598, 1000) / 4
    assert b0_over_w_gen < Fraction(9, 10)

    _, log_13_5_hi = log_bounds_series(Fraction(27, 2), m=12)
    assert log_13_5_hi < Fraction(2603, 1000)
    cheb_gen = log_13_5_hi / (15 * log2_lo)
    assert cheb_gen < Fraction(251, 1000)

    _, log_N_gen_hi = log_bounds_series(Fraction(27 * 751), m=12)
    assert log_N_gen_hi < 10
    had_gen = Fraction(10) / (750 * Fraction(1039, 100))
    assert had_gen < Fraction(2, 1000)

    total_gen = Fraction(1598, 1000) + Fraction(251, 1000) + Fraction(2, 1000)
    assert total_gen == Fraction(1851, 1000) < Fraction(119, 64)

    # Case C = 13 (base 2, p > 3^13 = 1594323 => p >= 1594325, g >= 797162)
    log_p_b2_lo, log_p_b2_hi = log_bounds_series(Fraction(1594325), m=12)
    assert Fraction(1428, 100) < log_p_b2_lo and log_p_b2_hi < Fraction(143, 10)
    A_g_b2_lo = Fraction(797161) * (log2_lo * log3_lo) / (log_p_b2_hi ** 2)
    assert A_g_b2_lo > 2950

    paren2 = Fraction(1000001, 1000000) + Fraction(627, 100000)
    assert Fraction(1) + Fraction(1, 3 * 17000 * 797161) < Fraction(1000001, 1000000) ** 2
    assert Fraction(2 + Fraction(1, 797161), 3 * 17000) < Fraction(627, 100000) ** 2
    lhs_sq_b2 = Fraction(36, 13) * Fraction(17001, 17000) * (paren2 ** 2)
    assert lhs_sq_b2 < Fraction(1676, 1000) ** 2

    log_1001_lo, _ = log_bounds_series(Fraction(1001, 1000), m=12)
    assert (2 * log2_hi) / 17001 < log_1001_lo
    b0_over_w_b2 = Fraction(1001, 1000) * Fraction(2241, 1000) * Fraction(1676, 1000) / 4
    assert b0_over_w_b2 < 1

    _, log_13_hi = log_bounds_series(Fraction(13), m=12)
    assert log_13_hi < Fraction(2565, 1000)
    cheb_b2 = log_13_hi / (13 * log3_lo)
    assert cheb_b2 < Fraction(180, 1000)

    _, log_N_b2_hi = log_bounds_series(Fraction(27 * 17001), m=12)
    assert log_N_b2_hi < 14
    had_b2 = Fraction(14) / (17000 * Fraction(1428, 100))
    assert had_b2 < Fraction(1, 1000)

    total_b2 = Fraction(1676, 1000) + Fraction(180, 1000) + Fraction(1, 1000)
    assert total_b2 == Fraction(1857, 1000) < Fraction(119, 64)

    print("[OK] All rational inequalities in Section 4 verified with exact Fractions:")
    print(f"     General / Quadratic case (C = 15): {float(total_gen):.3f} = 1851/1000 < 119/64 = {119/64:.6f}")
    print(f"     Base-2 / Quadratic case  (C = 13): {float(total_b2):.3f} = 1857/1000 < 119/64 = {119/64:.6f}")
    print()


# ==============================================================================
# PART 2: Benchmark Table — Elementary vs. Yamada vs. Giraudin vs. Theorem 1.1
# ==============================================================================

def benchmark_bounds() -> None:
    print("=" * 88)
    print("2. UPPER BOUNDS FOR v_p(2^(p-1) - 1) ON PRIMES p == +-1 (mod 24)")
    print("=" * 88)
    header = (
        f"{'Prime p':>17} | {'Elementary':>12} | {'Yamada (283)':>14} | "
        f"{'Giraudin (69)':>14} | {'Thm 1.1 (13)':>13} | {'Thm 1.1 (13/2)':>14}"
    )
    print(header)
    print("-" * 88)

    targets = [3**13 + 1, 10**7, 10**9, 10**12]
    primes_24 = [next_prime_pm1_mod24(t) for t in targets]

    for p in primes_24:
        assert is_prime(p) and p % 24 in (1, 23)
        elem = math.ceil((p - 1) * math.log(2) / math.log(p)) - 1
        # In base 2 with y = 3: log(y)*log(xy) = log(3)*log(6)
        yamada = math.floor(283 * (p - 1) * math.log(3) * math.log(6) / (math.log(p) ** 2)) + 4
        giraudin = math.floor(69 * (p - 1) * math.log(3) * math.log(6) / (math.log(p) ** 2)) + 4
        thm_13 = math.floor(13 * (p - 1) * math.log(2) * math.log(3) / (math.log(p) ** 2)) + 2
        thm_half = math.floor(6.5 * (p - 1) * math.log(2) * math.log(3) / (math.log(p) ** 2)) + 2
        print(
            f"{p:>17,d} | {elem:>12,d} | {yamada:>14,d} | "
            f"{giraudin:>14,d} | {thm_13:>13,d} | {thm_half:>14,d}"
        )
    print()


# ==============================================================================
# PART 3: Odd Perfect Number Selective Weight Demo (Theorem 5.1)
# ==============================================================================

def benchmark_opn_edge() -> None:
    print("=" * 88)
    print("3. SELECTIVE WEIGHT B(p_j, p_i) IN THEOREM 5.1 (p_i == 7 mod 8, j != *)")
    print("=" * 88)
    p_i = 1_000_039
    p_j = 10_000_019
    assert is_prime(p_i) and is_prime(p_j) and p_i % 8 == 7
    f_ij = (p_i - 1) // 2  # Odd order dividing p_i - 1 => automatic QR by Euler!

    elem_ord = math.ceil(f_ij * math.log(p_j) / math.log(p_i)) - 1
    yamada_edge = math.floor(283 * (p_i - 1) * math.log(2) * math.log(2 * p_j) / (math.log(p_i) ** 2)) + 4
    giraudin_edge = math.floor(69 * (p_i - 1) * math.log(2) * math.log(2 * p_j) / (math.log(p_i) ** 2)) + 4
    hex_15 = math.floor(15 * (p_i - 1) * math.log(2) * math.log(p_j) / (math.log(p_i) ** 2)) + 2
    hex_half = math.floor(7.5 * (p_i - 1) * math.log(2) * math.log(p_j) / (math.log(p_i) ** 2)) + 2

    print(f"Active pair (p_j -> p_i) = ({p_j:,d} -> {p_i:,d}), p_i == 7 (mod 8), f_{{i,j}} = (p_i - 1)/2:")
    print(f"  - Elementary order bound           : {elem_ord:>12,d}")
    print(f"  - Yamada (2010) [283, log(2*p_j)]  : {yamada_edge:>12,d}")
    print(f"  - Giraudin (2026) [69, log(2*p_j)] : {giraudin_edge:>12,d}")
    print(f"  - Hexagonal bound (kappa = 15)     : {hex_15:>12,d}")
    print(f"  - Euler-QR bound  (kappa = 15/2)   : {hex_half:>12,d}  <-- automatic for j != *!")
    print("=" * 88)


if __name__ == "__main__":
    certify_section_4()
    benchmark_bounds()
    benchmark_opn_edge()
