# RegisterAlgebra

[![tests](https://github.com/rkalendar/RegisterAlgebra/actions/workflows/tests.yml/badge.svg)](https://github.com/rkalendar/RegisterAlgebra/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24434511.svg)](https://doi.org/10.5281/zenodo.24434511)
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
| `rau_bits.py`, `bench_bits.py`, `tests/test_rau_bits.py` | Shift-and-mask variants of the register operations, their benchmark and their tests — not part of the paper, see [below](#beyond-the-paper-shifts-masks-and-machine-words) |
| `rau_word.c` | The same operations in C with every register held in one 64-bit machine word — not part of the paper |
| `results/bench_bits_log.txt`, `results/rau_word_log.txt` | Output of the two benchmarks above |

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
from rau import canonicalize, value, reg_add, reg_mul, reg_mul_norm

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

## Beyond the paper: shifts, masks and machine words

The files in this section were added after release 1.1.0. They are not part of
the manuscript and not in the archived version: every number in the paper comes
from `rau.py` and `reproduce_results.py`, and `rau.py` is unchanged. They answer
a follow-up question: how much do bit operations help, and can a register
operation be brought down to a single machine instruction?

**Python: `rau_bits.py`.** Each function returns exactly the pair its
counterpart in `rau.py` returns. Scaling by 2 is a left shift and scaling by 3 a
shift and an addition; halving and parity are `N >> 1` and `N & 1`; and for
canonical operands, whose second registers are 0 or 1, the products `a*d`, `b*c`
and `b*d` become the masks `a & -d`, `c & -b` and `b & d`, so that a single
multiplication remains. Time of the normalised product of two canonical pairs,
relative to native integer multiplication (median of seven runs):

| Operation | 64-bit | 1,024-bit | 16,384-bit |
| --- | --- | --- | --- |
| `reg_mul_norm` (reference, `rau.py`) | 5.40 | 1.56 | 1.06 |
| `reg_mul_canon` (shifts and masks, one multiplication) | 3.85 | 1.37 | 1.05 |
| `reg_mul_via_value` (decode, multiply, encode) | 4.86 | 1.45 | 1.04 |

Masks save about 30% at 64 bits and about 12% at 1,024 bits; at 16,384 bits the
one long multiplication dominates and the variants cannot be told apart.
Replacing only the scalings by 2 and 3 with shifts changes nothing measurable.
`canonicalize_bits` is up to four times faster than `canonicalize`, but that
step is a small fraction of a product.

**C: `rau_word.c`.** Every register is one signed 64-bit word and operand values
stay below 2^30, so that all products fit. Cost per operation in core cycles
(fastest of 15 runs; the three runs made agree to within 0.1 cycle):

| Operation | Throughput, cycles | Latency, cycles |
| --- | --- | --- |
| Native `x * y` | 1.1 | 3.0 |
| Formula (5), then normalisation (3d) | 5.6 | 10.8 |
| Shifts and masks, one multiplication | 3.4 | 7.2 |
| Decode, multiply, encode | 2.4 | 8.1 |
| Native `x + y` | 1.1 | 1.1 |
| Register sum (4), then normalisation (3d) | 2.8 | 2.7 |

Throughput is measured on independent operations, latency on a chain in which
each result is an operand of the next operation. No register operation runs in
one cycle. The best variants remain 2.2 to 2.7 times slower than the native
instruction: the register product contains the same machine multiplication,
preceded and followed by the instructions that encode and decode the pair.
Word-size registers also bound the range, and a chain of raw products overflows
at the sixth squaring of κ(1) = (−1, 1), so only normalised chains are timed.

```bash
python3 bench_bits.py                                        # about four minutes
gcc -O2 -fno-tree-vectorize -o rau_word rau_word.c && ./rau_word   # about one minute, x86-64
```

Both programs verify their variants against the reference before timing
anything: `rau_word` checks 13.8 million cases against the native result, and
`./rau_word --dump` prints vectors that can be compared with `rau.py`. The
figures above were measured on the machine named at the head of each log in
`results/`. The time-stamp counter of that machine ticks at its nominal rate, so
cycles were obtained by calibration, taking the latency of a native 64-bit
multiplication as 3 cycles; the ratios to native do not depend on that
assumption.

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

Software (version 1.1.0, GitHub release `v.1.1.0`):

> Kalendar, R. (2026). *RegisterAlgebra: reference implementation of exact
> register arithmetic* (Version 1.1.0) [Computer software]. Zenodo.
> https://doi.org/10.5281/zenodo.24434511

The DOI above identifies version 1.1.0, the version described in the revised
manuscript. Version 1.0.0 is archived at
https://doi.org/10.5281/zenodo.21630725. Releases are also preserved in
[Software Heritage](https://archive.softwareheritage.org/).

## License

MIT — see `LICENSE`.
