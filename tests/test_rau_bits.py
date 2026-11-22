"""The shift-and-mask variants must return exactly the pairs rau.py returns."""
import os, random, sys, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from rau import canonicalize, normalize, value, reg_mul, reg_mul_norm
from rau_bits import (canonicalize_bits, normalize_bits, value_bits,
                      reg_mul_bits, reg_mul_canon, reg_mul_via_value)

SEED = 20260713


class BitVariants(unittest.TestCase):
    def test_canonicalize_and_value(self):
        for n in range(-20000, 20001):
            self.assertEqual(canonicalize_bits(n), canonicalize(n))
        for a in range(-40, 41):
            for b in range(-40, 41):
                self.assertEqual(value_bits((a, b)), value((a, b)))
                self.assertEqual(normalize_bits((a, b)), normalize((a, b)))

    def test_raw_product_all_small_pairs(self):
        regs = [(a, b) for a in range(-9, 10) for b in range(-9, 10)]
        for X in regs:
            for Y in regs:
                self.assertEqual(reg_mul_bits(X, Y), reg_mul(X, Y))
                self.assertEqual(reg_mul_via_value(X, Y), reg_mul_norm(X, Y))

    def test_canonical_product(self):
        for m in range(-300, 301):
            for n in range(-300, 301):
                X, Y = canonicalize(m), canonicalize(n)
                self.assertEqual(reg_mul_canon(X, Y), reg_mul_norm(X, Y))
                self.assertEqual(reg_mul_canon(X, Y), canonicalize(m * n))

    def test_random_large(self):
        rng = random.Random(SEED)
        for bits in (64, 1024, 16384):
            for _ in range(300):
                m = rng.getrandbits(bits) - (1 << (bits - 1))
                n = rng.getrandbits(bits) - (1 << (bits - 1))
                X, Y = canonicalize(m), canonicalize(n)
                self.assertEqual(reg_mul_canon(X, Y), canonicalize(m * n))
                self.assertEqual(reg_mul_bits(X, Y), reg_mul(X, Y))
                U = (rng.getrandbits(bits) - (1 << (bits - 1)), rng.getrandbits(bits) - (1 << (bits - 1)))
                V = (rng.getrandbits(bits) - (1 << (bits - 1)), rng.getrandbits(bits) - (1 << (bits - 1)))
                self.assertEqual(reg_mul_bits(U, V), reg_mul(U, V))
                self.assertEqual(reg_mul_via_value(U, V), reg_mul_norm(U, V))
                self.assertEqual(normalize_bits(U), normalize(U))

    def test_canonical_variant_needs_canonical_operands(self):
        # documented precondition: for other representatives the result is wrong
        X, Y = (-6, 7), (0, 2)                       # values 9 and 6
        self.assertEqual(value(reg_mul_via_value(X, Y)), 54)
        self.assertNotEqual(value(reg_mul_canon(X, Y)), 54)


if __name__ == "__main__":
    unittest.main()
