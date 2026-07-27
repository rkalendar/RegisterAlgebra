"""rau.py — Reference implementation of the Register Arithmetic Unit (RAU).

Faithful, line-for-line transcription of Algorithms 1-3 of the paper

    "Register Algebra: Numbers as Linear Functions of Two Bases for
     Exact Integer and Rational Arithmetic"  (R. Kalendar, 2026)

A register integer is a pair (a, b) of Python integers with semantic
value  2*a + 3*b   (formula (1)).  A register fraction is a pair of
register pairs ((P, Q)) whose semantic value is value(P)/value(Q); the
quotient itself is never computed (Section 6, "image of division").

Only arbitrary-precision integer arithmetic is used: the module contains
no floating-point operation of any kind.

Supplementary Material S1.  Requires Python >= 3.8; no dependencies.
"""

# --------------------------------------------------------------------------
# Algorithm 1 — fundamental operations on register integers  (Section 8.1)
# --------------------------------------------------------------------------
# The eleven lines marked with  #core  constitute the complete operation
# set of the formalism (cf. Table 8, "~16 lines of code").

def value(R):                                            # formula (1)  #core
    return 2 * R[0] + 3 * R[1]                                          #core

def canonicalize(N):                                     # formula (3)  #core
    if N % 2 == 0:                                                      #core
        return (N // 2, 0)                               # (3a)         #core
    return ((N - 3) // 2, 1)                             # (3b)         #core

def reg_add(X, Y):                                       # formula (4)  #core
    return (X[0] + Y[0], X[1] + Y[1])                                   #core

def reg_sub(X, Y):                                                      #core
    return (X[0] - Y[0], X[1] - Y[1])                                   #core

def reg_mul(X, Y):                                       # formula (5)  #core
    a1, b1 = X
    a2, b2 = Y
    return (2 * a1 * a2, 2 * (a1 * b2 + b1 * a2) + 3 * b1 * b2)         #core


# --------------------------------------------------------------------------
# Algorithm 2 — normalisation and GCD reduction  (Section 8.2)
# --------------------------------------------------------------------------

def normalize(R):
    """Reduce any register pair to canonical form (b in {0, 1})."""
    return canonicalize(value(R))

from math import gcd as _gcd    # C-implemented Euclidean gcd — the same
                                 # primitive Fraction uses, so benchmarks
                                 # isolate the register-layer overhead.

def gcd(x, y):
    """Euclidean algorithm (explicit loop, as spelled in Algorithm 2)."""
    x = abs(x); y = abs(y)
    while y != 0:
        x, y = y, x % y
    return x

def reduce_fraction(P, Q):
    """Reduce the register fraction P/Q to lowest terms.

    The single point of the formalism at which an integer value is
    explicitly evaluated (Section 6.4) — strictly for cancellation.
    """
    num = value(P)
    den = value(Q)
    g = _gcd(abs(num), abs(den))
    if den < 0:                       # ensure a positive denominator
        g = -g
    return (canonicalize(num // g), canonicalize(den // g))


# --------------------------------------------------------------------------
# Algorithm 3 — rational arithmetic over register fractions  (Section 8.3)
# --------------------------------------------------------------------------
# A register fraction is represented as the pair (P, Q) of register pairs.
# Algorithm 3 reduces after every operation; the *_lazy variants apply
# formulas (7a)-(7d) verbatim without reduction (Section 11.3 of the
# paper studies the difference).

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
    assert value(Y[0]) != 0, "division by a zero-valued fraction"
    return reduce_fraction(reg_mul(X[0], Y[1]), reg_mul(X[1], Y[0]))

def frac_eq(X, Y):                                       # formula (6)
    """Equality by cross-multiplication — no division is performed."""
    return value(reg_mul(X[0], Y[1])) == value(reg_mul(X[1], Y[0]))

# ---- lazy variants: formulas (7a)-(7d) without reduce_fraction ----------

def frac_add_lazy(X, Y):
    return (reg_add(reg_mul(X[0], Y[1]), reg_mul(Y[0], X[1])),
            reg_mul(X[1], Y[1]))

def frac_sub_lazy(X, Y):
    return (reg_sub(reg_mul(X[0], Y[1]), reg_mul(Y[0], X[1])),
            reg_mul(X[1], Y[1]))

def frac_mul_lazy(X, Y):
    return (reg_mul(X[0], Y[0]), reg_mul(X[1], Y[1]))

def frac_div_lazy(X, Y):
    return (reg_mul(X[0], Y[1]), reg_mul(X[1], Y[0]))


# --------------------------------------------------------------------------
# Convenience constructors (interface code, not part of the operation set)
# --------------------------------------------------------------------------

def from_int(n):
    """Encode a Python integer as a canonical register pair."""
    return canonicalize(n)

def from_ratio(p, q):
    """Encode the rational p/q as a register fraction (image of division)."""
    if q == 0:
        raise ZeroDivisionError("zero denominator")
    return (canonicalize(p), canonicalize(q))

def to_ratio(X):
    """Return (numerator, denominator) of a register fraction, reduced."""
    P, Q = reduce_fraction(X[0], X[1])
    return value(P), value(Q)


if __name__ == "__main__":
    # The worked example of Section 9:  (1/3 + 1/6) * 2  =  1
    X = from_ratio(1, 3)
    Y = from_ratio(1, 6)
    Z = from_ratio(2, 1)
    R = frac_mul(frac_add(X, Y), Z)
    assert to_ratio(R) == (1, 1)
    print("(1/3 + 1/6) x 2 =", "%d/%d" % to_ratio(R), "- Section 9 trace OK")


# --------------------------------------------------------------------------
# Optimised rational operations (Section 11.4 of the paper)
# --------------------------------------------------------------------------
# Algorithm 3 above follows formulas (7a)-(7d) verbatim and therefore
# reduces the full cross-multiplied pair.  Because register fractions are
# ordinary fractions under the isomorphism of Theorem 5, the classical
# optimisations of rational arithmetic (Knuth, TAOCP vol. 2, sec. 4.5.1
# [11]) transfer unchanged: cancel by gcds of the *inputs* before
# multiplying, so that the expensive full-size gcd is never formed.

def frac_add_opt(X, Y, _sub=False):
    p1, q1 = value(X[0]), value(X[1])
    p2, q2 = value(Y[0]), value(Y[1])
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
    g1 = _gcd(p1, q2); g2 = _gcd(p2, q1)
    num = (p1 // g1) * (p2 // g2)
    den = (q1 // g2) * (q2 // g1)
    if den < 0:
        num, den = -num, -den
    return (canonicalize(num), canonicalize(den))

def frac_div_opt(X, Y):
    p2 = value(Y[0])
    assert p2 != 0, "division by a zero-valued fraction"
    return frac_mul_opt(X, (Y[1], Y[0]))
