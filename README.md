# RegisterAlgebra

[![tests](https://github.com/rkalendar/RegisterAlgebra/actions/workflows/tests.yml/badge.svg)](https://github.com/rkalendar/RegisterAlgebra/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21630724.svg)](https://doi.org/10.5281/zenodo.24432373)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

🌐 **Online version:** <https://digitalgens.org/algebra.html>

Reference implementation of **Register Algebra** — an exact integer and rational
arithmetic in which every integer `N` is represented as a pair of integer
registers `(a, b)` under the evaluation map

```
N = 2a + 3b
```

and every rational is kept as an *unevaluated image of division*: an ordered
pair of register pairs whose denominator has non-zero **value**, and whose
quotient is never formed. The registers are arbitrary-precision integers
(Python `int`), not fixed-width machine words. All operations are closed over
the integers, and the arithmetic itself (`rau.py`) executes no floating-point
instruction of any kind — floats occur only in `reproduce_results.py`, where
wall-clock timings are measured. Exact integer division is used: halving in
`canonicalize`/`normalize`, the remainder steps of the Euclidean algorithm, and
the exact division by the gcd in `reduce_fraction`.

This code accompanies the manuscript:

> R. Kalendar. *Register Algebra: Numbers as Linear Functions of Two Bases for
> Exact Integer and Rational Arithmetic.* (under review, 2026)

<!-- Once the article is published, replace the line above with the full
     reference and add:  DOI: https://doi.org/XX.XXXX/XXXXX  -->

## Contents

| File | Description |
| --- | --- |
| `rau.py` | The Register Arithmetic Unit: transcription of Algorithms 1–3 of the paper, plus the lazy variants of formulas (7a)–(7d) and the classical rational operations of Section 11.4 |
| `reproduce_results.py` | Regenerates every number reported in Section 11 (Tables 9–12) from a fixed seed |
| `tests/test_rau.py` | Fast test suite (runs in a few seconds; used by continuous integration) |
| `results/section11_reference_log.txt` | Frozen output of the full experiment suite, as reported in the paper |

Requires Python ≥ 3.8. There are no dependencies beyond the standard library.
Continuous integration runs the worked example, the test suite, and experiments
`e1`–`e4` (including `e2t`) on CPython 3.8, 3.10, 3.12, and 3.13.

## Quick start

```python
from rau import from_ratio, frac_add, frac_mul, to_ratio

x = from_ratio(1, 3)          # image of division 1/3
y = from_ratio(1, 6)          # image of division 1/6
z = from_ratio(2, 1)

r = frac_mul(frac_add(x, y), z)
print(to_ratio(r))            # (1, 1) — exact; the quotient 1/3 is never formed
```

The integer layer is deliberately small:

```python
from rau import canonicalize, value, reg_add, reg_mul

X = canonicalize(7)           # (2, 1),  since 2*2 + 3*1 = 7
Y = canonicalize(-4)          # (-2, 0)
value(reg_add(X, Y))          # 3
value(reg_mul(X, Y))          # -28
reg_mul(X, Y)                 # (-8, -4): the raw product is NOT canonical
reg_mul_norm(X, Y)            # (-14, 0): normalised product
```

A register fraction is *admissible* only if its denominator has non-zero value.
`Q != (0, 0)` is not enough — `(3, -2)` has value 0 — and every entry point
raises `ZeroDivisionError` for such a denominator.

## The operation set

Everything `rau.py` exports, in the order the paper introduces it. A *register
integer* is a tuple `(a, b)`; a *register fraction* is a tuple of two of them.

| Layer | Functions | Notes |
| --- | --- | --- |
| Integers (Algorithm 1) | `value`, `canonicalize`, `reg_add`, `reg_sub`, `reg_mul` | Formulas (1), (3)–(5); thirteen lines in total. `reg_mul` returns the raw, generally non-canonical product |
| Normalisation (Algorithm 2) | `normalize`, `reg_mul_norm`, `gcd`, `reduce_fraction`, `is_admissible`, `is_normal` | `normalize` works in register space, formula (3d), without evaluating the pair; `reduce_fraction` evaluates the pair, strictly for cancellation, and rejects zero-valued denominators |
| Rationals (Algorithm 3) | `frac_add`, `frac_sub`, `frac_mul`, `frac_div`, `frac_eq` | Formulas (6b), (7a)–(7d); any admissible operands in, normal form out |
| Lazy variants | `frac_add_lazy`, `frac_sub_lazy`, `frac_mul_lazy`, `frac_div_lazy` | The same formulas without reduction — the growth studied by `e4`; zero-valued denominators and zero divisors are rejected here too |
| Classical variants (Section 11.4) | `frac_add_opt`, `frac_sub_opt`, `frac_mul_opt`, `frac_div_opt` | Decode the registers, run Knuth's algorithm with input-side cancellation on plain integers, re-encode. They do **not** execute formulas (7a)–(7d) and require operands in normal form |
| Constructors | `from_int`, `from_ratio`, `from_ratio_reduced`, `as_ratio`, `to_ratio` | Interface code, not part of the operation set. `from_ratio` encodes verbatim (unreduced); `from_ratio_reduced` returns normal form, like `Fraction(p, q)` |

## Reproducing the results of the paper

```bash
python3 rau.py                     # runs the worked example of Section 9
python3 reproduce_results.py       # full suite, about four minutes on one core
python3 reproduce_results.py e1 e4 # any subset of the experiments
python3 -m unittest discover tests # fast test suite
```

The suite performs the following, with the seed fixed at `20260713`:

| Experiment | What it establishes |
| --- | --- |
| `e1` | Exhaustive checks over every ordered pair drawn from the 961 register pairs with \|a\|, \|b\| ≤ 15 (923,521 combinations): homomorphism identities (4), (5), register-space normalisation (3d) and the normalised product; the canonical form on every integer in [−10⁵, 10⁵]; the closed forms of Table 3; ring properties of raw pairs — 4,803,343 checks |
| `e2` | Differential testing against Python's exact `fractions.Fraction` on 125,000 random operand sets of 32–256 bits, zero numerators included: Algorithm 3 with normal-form output, the equality predicate (6b), all four lazy operations, re-encoded operands, and the classical variants — 2,250,000 checks |
| `e2t` | Targeted cases: zero-valued denominators such as (3, −2), zero divisors, distinct representatives of zero, equal fractions with different encodings — 29,160 checks |
| `e3` | Literal integer-operation counts of the register operations in the unit-cost model (Table 5) |
| `e4` | Register growth without reduction: Θ(n log n) bits against Θ(n) for the reduced form |
| `e5` | Timings: isolated products by parity class, raw and normalised (Table 10); chains of products under three normalisation policies (Table 11); a chain of rational operations under matched invariants (Table 12). Every measurement is the fastest of three timed loops, repeated in seven independent runs; medians and ranges are reported |

Together, 7,082,503 individual checks, with zero discrepancies. Absolute timings
and, to a lesser degree, ratios depend on the machine; the processor and
interpreter are printed at the head of `e5`. The reference run in `results/`
records the machine on which it was produced.

## What this is, and what it is not

Register Algebra is an alternative encoding of ordinary exact arithmetic, not a
faster arithmetic and not a new algebraic structure: in canonical form the pair
(a, b) carries the same information as (⌊N/2⌋, N mod 2). The measurements in
this repository make the trade-offs explicit:

- Register multiplication costs 7 integer multiplications and 2 additions per
  product in the unit-cost model, of which four are products of two registers.
  These are operation counts, not running times: the registers are
  arbitrary-precision integers and every operation costs at least time linear
  in their length. For canonical operands (`b ∈ {0,1}`) only one full-size
  product remains, so at 16,384-bit operands canonical register multiplication
  cannot be distinguished from native multiplication within run-to-run
  variation; for balanced non-canonical operands the cost is four to five times
  native.
- The choice of representative matters. `reg_mul` applied to canonical pairs
  returns a non-canonical pair, and a chain of raw products is not bounded by
  the values it represents: squaring κ(1) = (−1, 1) twenty times with `reg_mul`
  yields a pair that still represents 1 and whose first register has 1,048,576
  bits. Normalising after every product (`reg_mul_norm`) keeps the pair
  canonical at the price of a linear-time step.
- A chain of 2,000 mixed rational operations executed with formulas (7a)–(7d)
  and `reduce_fraction` is about 60× slower than `fractions.Fraction`, because
  the full cross-multiplied pair is reduced at each step. The `*_opt` functions
  come within a factor of about 1.7 of `Fraction`, but they do so by decoding
  the registers and running the classical algorithm on plain integers; they
  measure the cost of the encoding as a storage format, not the speed of the
  register operations.

The value of the formalism is didactic: a small, fully executable example of a
non-injective encoding, its kernel, its canonical representatives and the
field-of-fractions construction — not performance.

## Citing

If you use this code, please cite both the software and the article. GitHub's
"Cite this repository" button renders the entry from `CITATION.cff`.

Software (all versions, resolves to the latest):

> Kalendar, R. (2026). *RegisterAlgebra: reference implementation of exact
> register arithmetic* [Computer software]. Zenodo.
> https://doi.org/10.5281/zenodo.21630724

The DOI above is the concept DOI; version 1.0.0 specifically is archived at
https://doi.org/10.5281/zenodo.21630725. Releases are also preserved in
[Software Heritage](https://archive.softwareheritage.org/).

## License

MIT — see `LICENSE`.
