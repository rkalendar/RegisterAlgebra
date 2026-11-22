"""rau_bits.py — shift-and-mask variants of the register operations.

Experimental companion to rau.py.  Every function here returns exactly the
same register pair as its counterpart in rau.py; only the way the result is
computed differs:

  * scaling by 2 is a left shift, scaling by 3 is a shift and an addition;
  * halving and parity are a right shift and a mask (N >> 1, N & 1);
  * for canonical operands the second registers are 0 or 1, so the products
    a*d, b*c and b*d become masks (a & -d, c & -b, b & d) and a single
    multiplication, a*c, remains.

Python's >> is a floor shift and & 1 is the parity also for negative
integers, so the functions are valid on all of Z.

These are constant-factor rewrites.  The registers are arbitrary-precision
integers: a shift, a mask or an addition of an n-bit register still costs
time linear in n, and the multiplication still costs M(n).  rau.py is left
untouched so that it stays identical to the version described in the paper.
"""

def canonicalize_bits(N):                       # = rau.canonicalize, formula (3c)
    r = N & 1
    return ((N >> 1) - r, r)

def normalize_bits(R):                          # = rau.normalize, formula (3d)
    b = R[1]
    k = b >> 1
    return (R[0] + (k << 1) + k, b & 1)

def value_bits(R):                              # = rau.value, formula (1)
    b = R[1]
    return (R[0] << 1) + (b << 1) + b

def reg_mul_bits(X, Y):                         # = rau.reg_mul, formula (5), raw
    """Raw product for arbitrary representatives.  The four products of two
    registers remain; the three scalings by 2, 2 and 3 become shifts."""
    a1, b1 = X
    a2, b2 = Y
    bb = b1 * b2
    return ((a1 * a2) << 1, ((a1 * b2 + b1 * a2) << 1) + (bb << 1) + bb)

def reg_mul_canon(X, Y):                        # = rau.reg_mul_norm on canonical pairs
    """Normalised product of two CANONICAL pairs (second registers in {0, 1}),
    by the closed forms of Table 3: one multiplication, the rest shifts and
    masks.  The result is canonical.  Not valid for other representatives."""
    a, b = X
    c, d = Y
    t = (a & -d) + (c & -b) + (b & d)           # a*d + c*b + b*d, no multiplication
    return (((a * c) << 1) + (t << 1) + t, b & d)

def reg_mul_via_value(X, Y):                    # = rau.reg_mul_norm, any representatives
    """Decode both pairs, multiply the integers, re-encode: the route by which
    the canonical product reduces to one native multiplication."""
    N = value_bits(X) * value_bits(Y)
    r = N & 1
    return ((N >> 1) - r, r)
