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

from rau import (value, canonicalize, normalize, gcd, reg_add, reg_sub, reg_mul,
                 reduce_fraction, frac_add, frac_sub, frac_mul, frac_div, frac_eq,
                 frac_add_lazy, frac_sub_lazy, frac_mul_lazy, frac_div_lazy,
                 frac_add_opt, frac_sub_opt, frac_mul_opt, frac_div_opt,
                 from_int, from_ratio, to_ratio)

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
        Z = from_ratio(0, 1)
        with self.assertRaises(AssertionError):
            frac_div(X, Z)


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
