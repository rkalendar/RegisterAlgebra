"""rau.py — Reference implementation of the Register Arithmetic Unit (RAU).

Transcription of Algorithms 1-3 of the paper

    "Register Algebra: Numbers as Linear Functions of Two Bases for
     Exact Integer and Rational Arithmetic"  (R. Kalendar, 2026)

A register integer is a pair (a, b) of Python integers with semantic
value  2*a + 3*b   (formula (1)).  A register fraction is a pair (P, Q)
of register pairs with value(Q) != 0, whose semantic value is
value(P)/value(Q); the quotient of numerator by denominator is never
formed (Section 6, "image of division").

The registers are arbitrary-precision integers (Python int), not
fixed-width machine words.  Only integer arithmetic is used: the module
contains no floating-point operation of any kind.  Exact integer
division does occur, in three places: division by 2 in canonicalize and
normalize, the remainder steps of the Euclidean algorithm, and the exact
division by the gcd in reduce_fraction.

Version 1.1 (revision).  Changes with respect to 1.0:
  * a register fraction is admissible only if value(Q) != 0 (previously
    only Q != (0, 0) was excluded); every entry point enforces this and
    raises ZeroDivisionError otherwise, including the lazy variants;
  * normalize works in register space (formula (3d)) and no longer
    evaluates the pair; reg_mul_norm is the normalised product of Table 3;
  * the input/output invariants of every rational operation are stated
    explicitly, and from_ratio_reduced builds a fraction in normal form.

Supplementary Material S1.  Requires Python >= 3.8; no dependencies.
"""

# --------------------------------------------------------------------------
# Algorithm 1 — fundamental operations on register integers  (Section 8.1)
# --------------------------------------------------------------------------
# The lines marked with  #core  constitute the complete operation set of
# the integer layer (cf. Table 8).

def value(R):                                            # formula (1)  #core
    return 2 * R[0] + 3 * R[1]                                          #core

def canonicalize(N):                                     # formula (3c) #core
    q, r = divmod(N, 2)                # N = 2q + r,  r in {0, 1}       #core
    return (q - r, r)                  # kappa(N) = (q - r, r)          #core

def reg_add(X, Y):                                       # formula (4)  #core
    return (X[0] + Y[0], X[1] + Y[1])                                   #core

def reg_sub(X, Y):                                                      #core
    return (X[0] - Y[0], X[1] - Y[1])                                   #core

def reg_mul(X, Y):                                       # formula (5)  #core
    """Raw product.  The result is in general NOT canonical, even when
    both operands are (Section 5.3)."""
    a1, b1 = X
    a2, b2 = Y
    return (2 * a1 * a2, 2 * (a1 * b2 + b1 * a2) + 3 * b1 * b2)         #core


# --------------------------------------------------------------------------
# Algorithm 2 — normalisation and GCD reduction  (Section 8.2)
# --------------------------------------------------------------------------

def normalize(R):                                        # formula (3d)
    """Canonical representative of the class of R, computed in register
    space:  (a, b) -> (a + 3k, b - 2k)  with  k = floor(b / 2).

    Equal to canonicalize(value(R)) (Proposition 2), but the pair is never
    evaluated: one halving of b, one scaling by 3 and one addition.
    """
    k, r = divmod(R[1], 2)
    return (R[0] + 3 * k, r)

def reg_mul_norm(X, Y):
    """Normalised product: normalize(reg_mul(X, Y)).  Canonical output for
    arbitrary operands; for canonical operands it agrees with the closed
    forms of Table 3."""
    return normalize(reg_mul(X, Y))

from math import gcd as _gcd    # C-implemented Euclidean gcd — the same
                                 # primitive Fraction uses, so benchmarks
                                 # isolate the register-layer overhead.

def gcd(x, y):
    """Euclidean algorithm (explicit loop, as spelled in Algorithm 2)."""
    x = abs(x); y = abs(y)
    while y != 0:
        x, y = y, x % y
    return x

def is_admissible(X):
    """A register fraction (P, Q) is admissible iff value(Q) != 0.

    Q != (0, 0) is NOT sufficient: every pair (3k, -2k) has value zero.
    """
    return value(X[1]) != 0

def _require_admissible(X):
    if value(X[1]) == 0:
        raise ZeroDivisionError("register fraction with a zero-valued denominator")

def reduce_fraction(P, Q):
    """Normal form of the admissible register fraction P/Q.

    Input : any register pairs P, Q with value(Q) != 0.
    Output: (P', Q') with P', Q' canonical, gcd(value(P'), value(Q')) = 1
            and value(Q') > 0.  Zero is returned as (0, 0)/(-1, 1) = 0/1.

    The integer values are evaluated here, strictly for cancellation.
    """
    num = value(P)
    den = value(Q)
    if den == 0:
        raise ZeroDivisionError("register fraction with a zero-valued denominator")
    g = _gcd(num, den)                # g >= 1 because den != 0
    if den < 0:                       # positive-denominator convention
        g = -g
    return (canonicalize(num // g), canonicalize(den // g))

def is_normal(X):
    """True iff X is in the normal form produced by reduce_fraction."""
    P, Q = X
    p, q = value(P), value(Q)
    return (P[1] in (0, 1) and Q[1] in (0, 1) and q > 0 and _gcd(p, q) == 1)


# --------------------------------------------------------------------------
# Algorithm 3 — rational arithmetic over register fractions  (Section 8.3)
# --------------------------------------------------------------------------
# A register fraction is represented as the pair (P, Q) of register pairs.
#
# Invariants of the three families of operations:
#
#   frac_*       input : any admissible fractions, any representatives
#                output: normal form (reduced, canonical, denominator > 0)
#   frac_*_lazy  input : any admissible fractions, any representatives
#                output: admissible, NOT reduced and NOT canonical in general
#   frac_*_opt   input : normal form (required)
#                output: normal form
#
# All of them raise ZeroDivisionError on a zero-valued denominator and on
# division by a zero-valued fraction.

def frac_add(X, Y):                                      # formula (7a)
    num = reg_add(reg_mul(X[0], Y[1]), reg_mul(Y[0], X[1]))
    den = reg_mul(X[1], Y[1])
    return reduce_fraction(num, den)

def frac_sub(X, Y):                                      # formula (7b)
    num = reg_sub(reg_mul(X[0], Y[1]), reg_mul(Y[0], X[1]))
    den = reg_mul(X[1], Y[1])
    return reduce_fraction(num, den)

def frac_mul(X, Y):                                      # formula (7c)
    return reduce_fraction(reg_mul(X[0], Y[0]), reg_mul(X[1], Y[1]))

def frac_div(X, Y):                                      # formula (7d)
    _require_admissible(Y)            # Y.Q moves to the numerator, so it
    if value(Y[0]) == 0:              # must be checked here explicitly
        raise ZeroDivisionError("division by a zero-valued register fraction")
    return reduce_fraction(reg_mul(X[0], Y[1]), reg_mul(X[1], Y[0]))

def frac_eq(X, Y):                                       # relation (6b)
    """Equality of values by cross-multiplication; no quotient is formed.

    Defined on admissible fractions only: on that domain the relation is
    an equivalence, whereas a zero-valued denominator would break
    transitivity (0/Q with value(Q) = 0 would equal every fraction).
    """
    _require_admissible(X); _require_admissible(Y)
    return value(reg_mul(X[0], Y[1])) == value(reg_mul(X[1], Y[0]))

# ---- lazy variants: formulas (7a)-(7d) without reduce_fraction ----------

def frac_add_lazy(X, Y):
    _require_admissible(X); _require_admissible(Y)
    return (reg_add(reg_mul(X[0], Y[1]), reg_mul(Y[0], X[1])),
            reg_mul(X[1], Y[1]))

def frac_sub_lazy(X, Y):
    _require_admissible(X); _require_admissible(Y)
    return (reg_sub(reg_mul(X[0], Y[1]), reg_mul(Y[0], X[1])),
            reg_mul(X[1], Y[1]))

def frac_mul_lazy(X, Y):
    _require_admissible(X); _require_admissible(Y)
    return (reg_mul(X[0], Y[0]), reg_mul(X[1], Y[1]))

def frac_div_lazy(X, Y):
    _require_admissible(X); _require_admissible(Y)
    if value(Y[0]) == 0:
        raise ZeroDivisionError("division by a zero-valued register fraction")
    return (reg_mul(X[0], Y[1]), reg_mul(X[1], Y[0]))


# --------------------------------------------------------------------------
# Convenience constructors (interface code, not part of the operation set)
# --------------------------------------------------------------------------

def from_int(n):
    """Encode a Python integer as a canonical register pair."""
    return canonicalize(n)

def from_ratio(p, q):
    """Encode p/q verbatim as an admissible register fraction (image of
    division).  The registers are canonical; the fraction is NOT reduced
    and the sign of q is kept."""
    if q == 0:
        raise ZeroDivisionError("zero denominator")
    return (canonicalize(p), canonicalize(q))

def from_ratio_reduced(p, q):
    """Encode p/q in normal form (reduced, denominator > 0), i.e. with the
    same construction-time invariant as fractions.Fraction(p, q)."""
    if q == 0:
        raise ZeroDivisionError("zero denominator")
    return reduce_fraction(canonicalize(p), canonicalize(q))

def as_ratio(X):
    """(numerator, denominator) exactly as stored — no reduction."""
    return value(X[0]), value(X[1])

def to_ratio(X):
    """(numerator, denominator) of a register fraction, reduced."""
    P, Q = reduce_fraction(X[0], X[1])
    return value(P), value(Q)


# --------------------------------------------------------------------------
# Classical rational arithmetic behind a register interface (Section 11.4)
# --------------------------------------------------------------------------
# These functions do NOT execute formulas (7a)-(7d).  They decode the four
# registers with value(), run the classical integer algorithm with gcd
# cancellation on the *inputs* (Knuth, TAOCP vol. 2, sec. 4.5.1 [11]) and
# re-encode the result with canonicalize().  They are included to measure
# what the register encoding costs on top of the classical algorithm.
#
# Precondition: both operands in normal form (as returned by
# from_ratio_reduced, reduce_fraction or any frac_* / frac_*_opt call).
# For operands that are not reduced the result has the correct value but
# need not be reduced: 2/4 + 1/3 would be returned as 10/12.

def frac_add_opt(X, Y, _sub=False):
    p1, q1 = value(X[0]), value(X[1])
    p2, q2 = value(Y[0]), value(Y[1])
    if q1 == 0 or q2 == 0:
        raise ZeroDivisionError("register fraction with a zero-valued denominator")
    g = _gcd(q1, q2)
    if g == 1:
        num, den = p1 * q2 + (-p2 if _sub else p2) * q1, q1 * q2
    else:
        t = p1 * (q2 // g) + (-p2 if _sub else p2) * (q1 // g)
        g2 = _gcd(t, g)
        num, den = t // g2, (q1 // g) * (q2 // g2)
    if den < 0:
        num, den = -num, -den
    return (canonicalize(num), canonicalize(den))

def frac_sub_opt(X, Y):
    return frac_add_opt(X, Y, _sub=True)

def frac_mul_opt(X, Y):
    p1, q1 = value(X[0]), value(X[1])
    p2, q2 = value(Y[0]), value(Y[1])
    if q1 == 0 or q2 == 0:
        raise ZeroDivisionError("register fraction with a zero-valued denominator")
    g1 = _gcd(p1, q2); g2 = _gcd(p2, q1)
    num = (p1 // g1) * (p2 // g2)
    den = (q1 // g2) * (q2 // g1)
    if den < 0:
        num, den = -num, -den
    return (canonicalize(num), canonicalize(den))

def frac_div_opt(X, Y):
    _require_admissible(Y)
    if value(Y[0]) == 0:
        raise ZeroDivisionError("division by a zero-valued register fraction")
    return frac_mul_opt(X, (Y[1], Y[0]))


if __name__ == "__main__":
    # The worked example of Section 9:  (1/3 + 1/6) * 2  =  1
    X = from_ratio(1, 3)
    Y = from_ratio(1, 6)
    Z = from_ratio(2, 1)
    R = frac_mul(frac_add(X, Y), Z)
    assert R == ((-1, 1), (-1, 1)) and to_ratio(R) == (1, 1)
    print("(1/3 + 1/6) x 2 =", "%d/%d" % to_ratio(R), "- Section 9 trace OK")
