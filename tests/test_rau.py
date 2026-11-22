"""Fast test suite for the Register Arithmetic Unit.

Runs in a few seconds and is intended for continuous integration; the
exhaustive experiments of the paper live in reproduce_results.py.

    python3 -m unittest discover tests
"""
import os
import random
import sys
import unittest
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from rau import (value, canonicalize, normalize, reg_mul_norm, gcd,
                 reg_add, reg_sub, reg_mul, reduce_fraction, is_admissible, is_normal,
                 frac_add, frac_sub, frac_mul, frac_div, frac_eq,
                 frac_add_lazy, frac_sub_lazy, frac_mul_lazy, frac_div_lazy,
                 frac_add_opt, frac_sub_opt, frac_mul_opt, frac_div_opt,
                 from_int, from_ratio, from_ratio_reduced, as_ratio, to_ratio)

SEED = 20260713


class CanonicalForm(unittest.TestCase):
    def test_recovers_value_and_is_canonical(self):
        for n in range(-5000, 5001):
            K = canonicalize(n)
            self.assertIn(K[1], (0, 1))
            self.assertEqual(value(K), n)
            self.assertEqual(normalize(K), K)      # idempotent

    def test_paper_examples(self):
        self.assertEqual(canonicalize(7), (2, 1))
        self.assertEqual(canonicalize(-1), (-2, 1))
        self.assertEqual(canonicalize(-3), (-3, 1))
        self.assertEqual(canonicalize(0), (0, 0))

    def test_large_magnitudes(self):
        for e in (64, 256, 1024):
            n = 2 ** e - 1
            self.assertEqual(value(canonicalize(n)), n)
            self.assertEqual(value(canonicalize(-n)), -n)

    def test_formulas_3a_3b_and_3c_agree(self):
        for n in range(-2000, 2001):
            expected = (n // 2, 0) if n % 2 == 0 else ((n - 3) // 2, 1)
            self.assertEqual(canonicalize(n), expected)

    def test_canonical_representative_is_unique(self):
        # no second pair with b in {0, 1} has the same value (Lemma 1)
        for n in range(-200, 201):
            hits = [(a, b) for b in (0, 1) for a in range(-400, 401)
                    if 2 * a + 3 * b == n]
            self.assertEqual(hits, [canonicalize(n)])

    def test_register_space_normalisation(self):
        # normalize never evaluates the pair, yet equals kappa(phi(R))
        for a in range(-40, 41):
            for b in range(-40, 41):
                self.assertEqual(normalize((a, b)), canonicalize(value((a, b))))


class IntegerHomomorphism(unittest.TestCase):
    def test_exhaustive_small(self):
        regs = [(a, b) for a in range(-8, 9) for b in range(-8, 9)]
        for X in regs:
            vx = value(X)
            for Y in regs:
                vy = value(Y)
                self.assertEqual(value(reg_add(X, Y)), vx + vy)
                self.assertEqual(value(reg_sub(X, Y)), vx - vy)
                self.assertEqual(value(reg_mul(X, Y)), vx * vy)

    def test_random_large(self):
        rng = random.Random(SEED)
        for _ in range(2000):
            X = (rng.getrandbits(128) - 2 ** 127, rng.getrandbits(128) - 2 ** 127)
            Y = (rng.getrandbits(128) - 2 ** 127, rng.getrandbits(128) - 2 ** 127)
            self.assertEqual(value(reg_mul(X, Y)), value(X) * value(Y))
            self.assertEqual(reg_mul_norm(X, Y), canonicalize(value(X) * value(Y)))

    def test_raw_product_is_not_canonical_table_3(self):
        # formula (5) returns the raw pair; Table 3 lists the normalised one
        self.assertEqual(reg_mul((1, 0), (0, 1)), (0, 2))        # 2 x 3, raw
        self.assertEqual(reg_mul_norm((1, 0), (0, 1)), (3, 0))   # 2 x 3, normalised
        for a in range(-20, 21):
            for c in range(-20, 21):
                self.assertEqual(reg_mul((a, 0), (c, 1)), (2*a*c, 2*a))
                self.assertEqual(reg_mul_norm((a, 0), (c, 0)), (2*a*c, 0))
                self.assertEqual(reg_mul_norm((a, 0), (c, 1)), (2*a*c + 3*a, 0))
                self.assertEqual(reg_mul_norm((a, 1), (c, 0)), (2*a*c + 3*c, 0))
                self.assertEqual(reg_mul_norm((a, 1), (c, 1)), (2*a*c + 3*a + 3*c + 3, 1))

    def test_repeated_squaring_of_one_grows_without_normalisation(self):
        # Section 5.3: every pair represents 1, the first register doubles in length
        acc, seen = (-1, 1), []
        for k in range(1, 11):
            acc = reg_mul(acc, acc)
            seen.append(acc)
            self.assertEqual(value(acc), 1)
            self.assertEqual(acc[0], 2 ** (2 ** k - 1))
            self.assertEqual(normalize(acc), (-1, 1))
        self.assertEqual(seen[:3], [(2, -1), (8, -5), (128, -85)])


class RationalArithmetic(unittest.TestCase):
    def test_section_9_trace(self):
        r = frac_mul(frac_add(from_ratio(1, 3), from_ratio(1, 6)), from_int_pair(2))
        self.assertEqual(to_ratio(r), (1, 1))

    def test_differential_against_fraction(self):
        rng = random.Random(SEED)
        for _ in range(3000):
            nb = rng.choice((32, 64, 128))
            p, q, r, s = (rng.getrandbits(nb) or 1 for _ in range(4))
            if rng.getrandbits(1):
                p = -p
            if rng.getrandbits(1):
                r = -r
            X, Y = from_ratio(p, q), from_ratio(r, s)
            F, G = Fraction(p, q), Fraction(r, s)
            self.assertEqual(to_ratio(frac_add(X, Y)), (F + G).as_integer_ratio())
            self.assertEqual(to_ratio(frac_sub(X, Y)), (F - G).as_integer_ratio())
            self.assertEqual(to_ratio(frac_mul(X, Y)), (F * G).as_integer_ratio())
            self.assertEqual(to_ratio(frac_div(X, Y)), (F / G).as_integer_ratio())
            self.assertEqual(frac_eq(X, Y), F == G)

    def test_lazy_variants_agree_semantically(self):
        rng = random.Random(SEED + 1)
        for _ in range(1000):
            p, q, r, s = (rng.getrandbits(64) or 1 for _ in range(4))
            X, Y = from_ratio(p, q), from_ratio(r, s)
            F, G = Fraction(p, q), Fraction(r, s)
            for lazy, expected in ((frac_add_lazy(X, Y), F + G),
                                   (frac_sub_lazy(X, Y), F - G),
                                   (frac_mul_lazy(X, Y), F * G),
                                   (frac_div_lazy(X, Y), F / G)):
                self.assertEqual(Fraction(value(lazy[0]), value(lazy[1])), expected)

    def test_optimised_variants_match_algorithm_3(self):
        rng = random.Random(SEED + 2)
        for _ in range(2000):
            p, q, r, s = (rng.getrandbits(96) or 1 for _ in range(4))
            if rng.getrandbits(1):
                p = -p
            if rng.getrandbits(1):
                r = -r
            X, Y = from_ratio(p, q), from_ratio(r, s)
            self.assertEqual(to_ratio(frac_add_opt(X, Y)), to_ratio(frac_add(X, Y)))
            self.assertEqual(to_ratio(frac_sub_opt(X, Y)), to_ratio(frac_sub(X, Y)))
            self.assertEqual(to_ratio(frac_mul_opt(X, Y)), to_ratio(frac_mul(X, Y)))
            self.assertEqual(to_ratio(frac_div_opt(X, Y)), to_ratio(frac_div(X, Y)))

    def test_denominator_sign_is_normalised(self):
        for p, q in ((1, -3), (-1, -3), (5, -10)):
            num, den = to_ratio(from_ratio(p, q))
            self.assertGreater(den, 0)
            self.assertEqual(Fraction(num, den), Fraction(p, q))

    def test_division_by_zero_valued_fraction_is_rejected(self):
        X = from_ratio(1, 2)
        for Z in (from_ratio(0, 1), ((3, -2), (-1, 1)), ((-6, 4), (1, 0))):
            for div in (frac_div, frac_div_lazy, frac_div_opt):
                with self.assertRaises(ZeroDivisionError):
                    div(X, Z)

    def test_zero_valued_denominator_is_rejected(self):
        # Q = (3, -2) is not (0, 0) but represents zero (Definition 2)
        self.assertEqual(value((3, -2)), 0)
        bad, half, one = ((-1, 1), (3, -2)), from_ratio(1, 2), from_ratio(1, 1)
        self.assertFalse(is_admissible(bad))
        with self.assertRaises(ZeroDivisionError):
            reduce_fraction((-1, 1), (3, -2))
        ops = (frac_add, frac_sub, frac_mul, frac_div, frac_eq,
               frac_add_lazy, frac_sub_lazy, frac_mul_lazy, frac_div_lazy,
               frac_add_opt, frac_sub_opt, frac_mul_opt, frac_div_opt)
        for op in ops:
            with self.assertRaises(ZeroDivisionError):
                op(bad, half)
            with self.assertRaises(ZeroDivisionError):
                op(half, bad)
        # the transitivity counter-example 0/1 ~ 0/(3,-2) ~ 1/1 cannot be formed
        with self.assertRaises(ZeroDivisionError):
            frac_eq(((0, 0), (3, -2)), one)

    def test_distinct_representatives_of_zero(self):
        zero, one, half = from_ratio(0, 1), from_ratio(1, 1), from_ratio(1, 2)
        for k in (1, -1, 5, -9):
            Z = ((3 * k, -2 * k), (2, 1))                 # 0/7
            self.assertTrue(frac_eq(Z, zero))
            self.assertFalse(frac_eq(Z, one))
            self.assertEqual(reduce_fraction(*Z), ((0, 0), (-1, 1)))
            self.assertEqual(frac_add(half, Z), half)

    def test_equal_fractions_with_different_encodings(self):
        rng = random.Random(SEED + 5)
        for _ in range(500):
            p = rng.randint(-10 ** 6, 10 ** 6); q = rng.randint(1, 10 ** 6)
            ref = from_ratio_reduced(p, q)
            self.assertTrue(is_normal(ref))
            self.assertEqual(as_ratio(ref), Fraction(p, q).as_integer_ratio())
            for m in (-5, -1, 2, 9):
                k1, k2 = rng.randint(-7, 7), rng.randint(-7, 7)
                a, b = canonicalize(p * m); c, d = canonicalize(q * m)
                E = ((a + 3 * k1, b - 2 * k1), (c + 3 * k2, d - 2 * k2))
                self.assertTrue(frac_eq(E, ref))
                self.assertEqual(reduce_fraction(*E), ref)
                self.assertEqual(frac_mul(E, E), frac_mul(ref, ref))

    def test_optimised_variants_return_normal_form(self):
        # normal-form inputs -> register tuples identical to Algorithm 3
        rng = random.Random(SEED + 6)
        for _ in range(1000):
            p, q, r, s = (rng.getrandbits(96) or 1 for _ in range(4))
            X, Y = from_ratio_reduced(-p, q), from_ratio_reduced(r, s)
            for opt, ref in ((frac_add_opt, frac_add), (frac_sub_opt, frac_sub),
                             (frac_mul_opt, frac_mul), (frac_div_opt, frac_div)):
                out = opt(X, Y)
                self.assertTrue(is_normal(out))
                self.assertEqual(out, ref(X, Y))
        # the precondition matters: unreduced inputs give an unreduced result
        self.assertEqual(as_ratio(frac_add_opt(from_ratio(2, 4), from_ratio(1, 3))), (10, 12))
        self.assertEqual(as_ratio(frac_add(from_ratio(2, 4), from_ratio(1, 3))), (5, 6))


class Helpers(unittest.TestCase):
    def test_gcd_matches_math_gcd(self):
        import math
        rng = random.Random(SEED + 3)
        for _ in range(2000):
            x, y = rng.getrandbits(64), rng.getrandbits(64)
            self.assertEqual(gcd(x, y), math.gcd(x, y))
            self.assertEqual(gcd(-x, y), math.gcd(x, y))

    def test_reduce_fraction_is_in_lowest_terms(self):
        import math
        rng = random.Random(SEED + 4)
        for _ in range(2000):
            p, q = rng.getrandbits(64) or 1, rng.getrandbits(64) or 1
            P, Q = reduce_fraction(canonicalize(p), canonicalize(q))
            self.assertEqual(math.gcd(value(P), value(Q)), 1)
            self.assertEqual(Fraction(value(P), value(Q)), Fraction(p, q))

    def test_from_int_round_trip(self):
        for n in (-2 ** 100, -7, 0, 1, 2, 3, 2 ** 100):
            self.assertEqual(value(from_int(n)), n)


def from_int_pair(n):
    """The integer n as a register fraction n/1."""
    return from_ratio(n, 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
