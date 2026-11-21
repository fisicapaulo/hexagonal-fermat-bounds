#!/usr/bin/env python3
"""Companion code for Paulo Vieira's hexagonal Fermat bounds.

Requires Python 3.9+; standard library only.
1. Checks the numerical inequalities used in Section 4.
2. Compares upper bounds for the base-2 Fermat numerator valuation.
3. Illustrates the weight B(p_j, p_i) from Section 5.

All numerical checks and table roundings use integers or Fraction.
"""

from fractions import Fraction
from functools import lru_cache
from math import ceil, floor

F = Fraction
MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
MR_LIMIT = 318_665_857_834_031_151_167_461


def require(condition: bool, message: str) -> None:
    """Check a condition even when Python is run with -O."""
    if not condition:
        raise AssertionError(message)


def is_prime(n: int) -> bool:
    """Deterministic Miller–Rabin for integers below MR_LIMIT."""
    if n < 2:
        return False
    if n >= MR_LIMIT:
        raise ValueError(f"Outside the deterministic range: {n}")
    for p in MR_BASES:
        if n % p == 0:
            return n == p

    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1

    for a in MR_BASES:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


def next_prime_pm1_mod24(start: int) -> int:
    """Return the smallest prime p >= start with p == +/-1 (mod 24)."""
    p = max(3, start)
    p += (p % 2 == 0)
    while True:
        if p % 24 in (1, 23) and is_prime(p):
            return p
        p += 2


def distinct_prime_factors(n: int) -> list[int]:
    factors = []
    q = 2
    while q * q <= n:
        if n % q == 0:
            factors.append(q)
            while n % q == 0:
                n //= q
        q = 3 if q == 2 else q + 2
    if n > 1:
        factors.append(n)
    return factors


def multiplicative_order(a: int, p: int) -> int:
    if not is_prime(p) or a % p == 0:
        raise ValueError("Expected a prime modulus and a nonzero residue.")
    order = p - 1
    for q in distinct_prime_factors(order):
        while order % q == 0 and pow(a, order // q, p) == 1:
            order //= q
    return order


@lru_cache(maxsize=None)
def log_bounds(z: Fraction, m: int = 12) -> tuple[Fraction, Fraction]:
    """Enclose log(z) by rational bounds from the series in Section 4."""
    r = F(z)
    if r < 1 or m < 1:
        raise ValueError("Expected z >= 1 and m >= 1.")

    h = 0
    while r >= 2:
        r /= 2
        h += 1

    def reduced_bounds(w: Fraction) -> tuple[Fraction, Fraction]:
        t = (w - 1) / (w + 1)
        t2 = t * t
        power = t
        partial = F(0)
        for j in range(m):
            partial += power / (2 * j + 1)
            power *= t2
        low = 2 * partial
        tail = 2 * power / ((2 * m + 1) * (1 - t2))
        return low, low + tail

    lo2, hi2 = reduced_bounds(F(2))
    lor, hir = reduced_bounds(r)
    return h * lo2 + lor, h * hi2 + hir


def exp_3_halves_upper() -> Fraction:
    """Upper bound for exp(3/2)/2 by a degree-12 sum and a tail."""
    term = total = F(1)
    for k in range(1, 13):
        term *= F(3, 2 * k)
        total += term
    first_omitted = term * F(3, 26)
    return (total + first_omitted / (1 - F(3, 28))) / 2


def certify_case(
    *, C, cutoff, g_min, log_p_min, log_p_cap, A_min, s_min,
    sqrt_bound, G_max, power_cap, b_ratio, y_min,
    log_y_min, log_cheb_cap, cheb_cap, log_cost_cap, cost_cap,
) -> Fraction:
    """Check one parameter regime, including its quadratic refinement."""
    require(cutoff >= 8, "Cutoff outside the monotonicity range")
    lo, hi = log_bounds(F(cutoff))
    require(log_p_min < lo and hi < log_p_cap, f"C={C}: cutoff log")
    require(F(cutoff - 1, 2) <= g_min <= cutoff // 2, f"C={C}: g lower bound")

    # The proof uses monotonicity of (t-1)/log(t)^2 for t >= 8.
    require(F(cutoff - 1, 2) * F("0.76") / log_p_cap**2 > A_min,
            f"C={C}: A_g lower bound")
    require(floor(F(4 * C * A_min, 9)) >= s_min, f"C={C}: s lower bound")

    first = F("1.000001")
    require(1 + F(1, 3 * s_min * g_min) < first**2,
            f"C={C}: first square root")
    require((2 + F(1, g_min)) / (3 * s_min) < sqrt_bound**2,
            f"C={C}: second square root")
    G_squared = F(36, C) * (1 + F(1, s_min)) * (first + sqrt_bound)**2
    require(G_squared < G_max**2 and G_max < 2, f"C={C}: bound for G")

    _, log2_hi = log_bounds(F(2))
    power_lo, _ = log_bounds(power_cap)
    require(2 * log2_hi / (s_min + 1) < power_lo, f"C={C}: power of two")
    require(power_cap * F("2.241") * G_max / 4 < b_ratio, f"C={C}: b_0/w")

    # For w >= C, log(b_ratio*w)/w is decreasing if log(b_ratio*C) > 1.
    cheb_lo, cheb_hi = log_bounds(b_ratio * C)
    require(1 < cheb_lo and cheb_hi < log_cheb_cap, f"C={C}: Chebyshev log")
    y_lo, _ = log_bounds(F(y_min))
    require(log_y_min < y_lo, f"C={C}: lower bound for log(y)")
    require(log_cheb_cap / (C * log_y_min) < cheb_cap, f"C={C}: Chebyshev cost")

    _, cost_hi = log_bounds(F(27 * (s_min + 1)))
    require(cost_hi < log_cost_cap, f"C={C}: combinatorial log")
    require(log_cost_cap / (s_min * log_p_min) < cost_cap,
            f"C={C}: combinatorial cost")

    total = G_max + cheb_cap + cost_cap
    require(total < F(119, 64), f"C={C}: final numerical criterion")
    return total


def certify_section_4() -> None:
    print("1. EXACT RATIONAL CHECKS FOR SECTION 4")
    log2_lo, log2_hi = log_bounds(F(2))
    log3_lo, log3_hi = log_bounds(F(3))
    require(F("0.693147") < log2_lo <= log2_hi < F("0.693148"), "log(2)")
    require(F("1.098612") < log3_lo <= log3_hi < F("1.098613"), "log(3)")
    require(log2_lo * log3_lo > F("0.76"), "Lower bound for log(2) log(3)")
    require(13 * log2_hi * log3_hi < F(99, 10), "Base-2 decimal coefficient")
    require(exp_3_halves_upper() < F("2.241"), "exp(3/2)/2")

    general = certify_case(
        C=15, cutoff=2**15, g_min=16_384,
        log_p_min=F("10.39"), log_p_cap=F("10.4"),
        A_min=113, s_min=750, sqrt_bound=F("0.02982"),
        G_max=F("1.598"), power_cap=F("1.002"), b_ratio=F("0.9"),
        y_min=2, log_y_min=F("0.693147"), log_cheb_cap=F("2.603"),
        cheb_cap=F("0.251"), log_cost_cap=F(10), cost_cap=F("0.002"),
    )
    base_two = certify_case(
        C=13, cutoff=3**13, g_min=797_161,
        log_p_min=F("14.28"), log_p_cap=F("14.3"),
        A_min=2950, s_min=17_000, sqrt_bound=F("0.00627"),
        G_max=F("1.676"), power_cap=F("1.001"), b_ratio=F(1),
        y_min=3, log_y_min=F("1.098612"), log_cheb_cap=F("2.565"),
        cheb_cap=F("0.180"), log_cost_cap=F(14), cost_cap=F("0.001"),
    )
    require(general == F(1851, 1000), "Unexpected general total")
    require(base_two == F(1857, 1000), "Unexpected base-2 total")
    print(f"   [OK] General and quadratic cases: {general} < 119/64")
    print(f"   [OK] Base 2 and quadratic cases:  {base_two} < 119/64")
    print()


def rounded_log_ratio(
    scale: Fraction, numerator: tuple[int, ...], p: int,
    denominator_power: int, rounding,
) -> int:
    """Round scale * product(log(numerator)) / log(p)^denominator_power.

    Return an integer only when both rational endpoints give that integer.
    """
    for m in (12, 24, 48, 96):
        plo, phi = log_bounds(F(p), m)
        low = scale / phi**denominator_power
        high = scale / plo**denominator_power
        for a in numerator:
            alo, ahi = log_bounds(F(a), m)
            low *= alo
            high *= ahi
        rounded_low, rounded_high = rounding(low), rounding(high)
        if rounded_low == rounded_high:
            return rounded_low

    raise ArithmeticError(
        f"Rounding not certified: scale={scale}, numerator={numerator}, "
        f"p={p}, denominator power={denominator_power}"
    )


def elementary_bound(x: int, p: int, exponent: int) -> int:
    return rounded_log_ratio(F(exponent), (x,), p, 1, ceil) - 1


def determinant_bound(
    coefficient: Fraction, a: int, b: int, p: int, additive: int,
) -> int:
    return rounded_log_ratio(coefficient * (p - 1), (a, b), p, 2, floor) + additive


def benchmark_bounds() -> None:
    print("2. BASE-2 UPPER BOUNDS (CERTIFIED INTEGER ROUNDING)")
    print(f"{'Prime p':>17} | {'Elementary':>14} | {'Yamada':>14} | "
          f"{'N. Giraudin':>14} | {'Hex. 13':>14} | {'Hex. 13/2':>14}")
    print("-" * 102)

    for start in (3**13 + 1, 10**7, 10**9, 10**12):
        p = next_prime_pm1_mod24(start)
        require(all(pow(a, (p - 1) // 2, p) == 1 for a in (2, 3)),
                f"Quadratic hypotheses failed at p={p}")
        values = (
            elementary_bound(2, p, p - 1),
            determinant_bound(F(283), 3, 6, p, 4),
            determinant_bound(F(69), 3, 6, p, 4),
            determinant_bound(F(13), 2, 3, p, 2),
            determinant_bound(F(13, 2), 2, 3, p, 2),
        )
        print(f"{p:>17,d} | " + " | ".join(f"{value:>14,d}" for value in values))
    print()


def benchmark_opn_edge() -> None:
    print("3. ILLUSTRATIVE WEIGHT FOR THE ODD PERFECT NUMBER APPLICATION")
    p_i = 1_000_039
    require(is_prime(p_i) and p_i % 8 == 7, "Invalid target prime")

    # Keep the original starting point; select the first prime with
    # a nontrivial odd order, as required for a non-Eulerian source.
    for p_j in range(10_000_019, 10_002_020, 2):
        if not is_prime(p_j):
            continue
        # Since p_i == 7 (mod 8), (p_i-1)/2 is odd.
        if pow(p_j, (p_i - 1) // 2, p_i) != 1:
            continue
        f_ij = multiplicative_order(p_j, p_i)
        if f_ij > 1 and f_ij % 2 == 1:
            break
    else:
        raise RuntimeError("No suitable illustrative source prime found.")

    require(pow(2, (p_i - 1) // 2, p_i) == 1, "Auxiliary base is not a square")
    require(pow(p_j, f_ij, p_i) == 1, "Invalid order")
    require(all(pow(p_j, f_ij // q, p_i) != 1
                for q in distinct_prime_factors(f_ij)), "Order is not minimal")

    elem = elementary_bound(p_j, p_i, f_ij)
    yamada = determinant_bound(F(283), 2, 2 * p_j, p_i, 4)
    giraudin = determinant_bound(F(69), 2, 2 * p_j, p_i, 4)
    hex_general = determinant_bound(F(15), 2, p_j, p_i, 2)
    hex_quadratic = determinant_bound(F(15, 2), 2, p_j, p_i, 2)
    weight = min(elem, hex_quadratic)

    print(f"   Source -> target: {p_j:,d} -> {p_i:,d}; actual order: {f_ij:,d}")
    print("   This is an admissible pair, not a claim about an existing odd perfect number.")
    for label, value in (
        ("Elementary bound using the actual order", elem),
        ("Yamada (283)", yamada),
        ("Niedbala Giraudin (69)", giraudin),
        ("Hexagonal bound (15)", hex_general),
        ("Hexagonal quadratic bound (15/2)", hex_quadratic),
        ("Weight B = min(elementary, quadratic)", weight),
    ):
        print(f"   {label:<44} {value:>12,d}")


if __name__ == "__main__":
    certify_section_4()
    benchmark_bounds()
    benchmark_opn_edge()
