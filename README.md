# RegisterAlgebra

[![tests](https://github.com/rkalendar/RegisterAlgebra/actions/workflows/tests.yml/badge.svg)](https://github.com/rkalendar/RegisterAlgebra/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21630724.svg)](https://doi.org/10.5281/zenodo.21630724)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

🌐 **Online version:** <[https://digitalgens.org/algebra.html](https://digitalgens.org/algebra.html)>

Reference implementation of **Register Algebra** — an exact integer and rational
arithmetic in which every integer `N` is represented as a pair of integer
registers `(a, b)` under the evaluation map

```
N = 2a + 3b
```

and every rational is kept as an *unevaluated image of division*: an ordered
pair of register pairs, whose quotient is never formed. All operations are
closed over the integers, and no floating-point instruction is executed at any
point in this repository.

This code accompanies the manuscript:

> R. Kalendar. *Register Algebra: Numbers as Linear Functions of Two Bases for
> Exact Integer and Rational Arithmetic.* (under review, 2026)

<!-- Once the article is published, replace the line above with the full
     reference and add:  DOI: https://doi.org/XX.XXXX/XXXXX  -->

## Contents

| File | Description |
| --- | --- |
| `rau.py` | The Register Arithmetic Unit: verbatim transcription of Algorithms 1–3 of the paper, plus the lazy variants of formulas (7a)–(7d) and the optimised rational operations of Section 11.4 |
| `reproduce_results.py` | Regenerates every number reported in Section 11 (Tables 9–11) from a fixed seed |
| `tests/test_rau.py` | Fast test suite (runs in a few seconds; used by continuous integration) |
| `results/section11_reference_log.txt` | Frozen output of the full experiment suite, as reported in the paper |

Requires Python ≥ 3.8. There are no dependencies beyond the standard library.

## Quick start

```python
from rau import from_ratio, frac_add, frac_mul, to_ratio

x = from_ratio(1, 3)          # image of division 1/3
y = from_ratio(1, 6)          # image of division 1/6
z = from_ratio(2, 1)

r = frac_mul(frac_add(x, y), z)
print(to_ratio(r))            # (1, 1) — exact, computed without division
```

The integer layer is deliberately small:

```python
from rau import canonicalize, value, reg_add, reg_mul

X = canonicalize(7)           # (2, 1),  since 2*2 + 3*1 = 7
Y = canonicalize(-4)          # (-2, 0)
value(reg_add(X, Y))          # 3
value(reg_mul(X, Y))          # -28
```

## Reproducing the results of the paper

```bash
python3 rau.py                     # runs the worked example of Section 9
python3 reproduce_results.py       # full suite, ~45 s on one core
python3 reproduce_results.py e1 e4 # any subset of the five experiments
python3 -m unittest discover tests # fast test suite
```

The suite performs the following, with the seed fixed at `20260713`:

| Experiment | What it establishes |
| --- | --- |
| `e1` | Exhaustive verification of the homomorphism identities (4), (5) over all 923,521 register pairs with \|a\|, \|b\| ≤ 15, and of the canonical form on every integer in [−10⁵, 10⁵] — 3,370,566 checks |
| `e2` | Differential testing of operations (6), (7a)–(7d) against Python's exact `fractions.Fraction` on 125,000 random operand sets of 32–256 bits — 875,000 checks |
| `e3` | Literal integer-operation counts of the register operations (Table 5) |
| `e4` | Register growth without reduction: Θ(n log n) bits against Θ(n) for the reduced form |
| `e5` | Microbenchmarks of register multiplication and of rational chains (Tables 10–11) |

Together, more than 4.2 million individual checks, with zero discrepancies.
Absolute timings depend on the machine; the ratios reported in the paper are
stable. The reference run in `results/` was produced on CPython 3.12.3, single
core.

## What this is, and what it is not

Register Algebra is a *structural* result, not a faster arithmetic. The
measurements in this repository make the trade-off explicit:

- Register multiplication costs 7 integer multiplications and 2 additions per
  product, of which four are full-precision. For canonical operands (`b ∈ {0,1}`)
  only one full-precision product remains, so at 16,384-bit operands register
  multiplication runs within 3 % of native multiplication; for balanced
  non-canonical operands the cost settles at the predicted factor of four.
- A chain of 2,000 mixed rational operations, executed with Algorithm 3 exactly
  as specified, is 64.6× slower than `fractions.Fraction`, because the algorithm
  reduces the full cross-multiplied pair at each step. Since register fractions
  are ordinary fractions under the isomorphism of Theorem 5, the classical
  optimisations of rational arithmetic (Knuth, *TAOCP* vol. 2, §4.5.1) transfer
  unchanged; the optimised variants included here close the gap to 2.0×.

The value of the formalism is exactness, closure over the integers, and a
transparent, directly implementable route to the field-of-fractions
construction — not performance.

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
