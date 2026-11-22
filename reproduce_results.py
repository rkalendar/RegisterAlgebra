"""reproduce_results.py — reproduces every number in Section 11 of the paper.

Usage:  python3 reproduce_results.py          (full suite, a few minutes on one core)
        python3 reproduce_results.py e1 e4    (any subset, in the given order)

Experiments
  E1  exhaustive checks on small registers: homomorphism identities,
      canonical form, register-space normalisation, Table 3, ring axioms
  E2  randomized differential testing against fractions.Fraction, plus
      targeted cases (zero divisors, representatives of zero, equal
      fractions with different encodings)
  E3  integer-operation counts (empirical confirmation of Table 5)
  E4  register growth: lazy formulas (7a)-(7d) vs Algorithm 3
  E5  timings: isolated products (Table 10), product chains under three
      normalisation policies (Table 11), rational chain (Table 12)

Timings are wall-clock and machine dependent.  Every measurement is the
fastest of three timed loops; it is repeated in RUNS independent runs, and
the median and the min-max range over those runs are reported.

Supplementary Material S1 (version 1.1). Pure Python, no dependencies.
"""
import platform, random, statistics, sys, time
from fractions import Fraction

from rau import (value, canonicalize, reg_add, reg_sub, reg_mul,
                 normalize, reg_mul_norm, gcd, reduce_fraction, is_normal,
                 frac_add, frac_sub, frac_mul, frac_div, frac_eq,
                 frac_add_lazy, frac_sub_lazy, frac_mul_lazy, frac_div_lazy,
                 frac_add_opt, frac_sub_opt, frac_mul_opt, frac_div_opt,
                 from_ratio, from_ratio_reduced, as_ratio, to_ratio)

SEED = 20260713
random.seed(SEED)
RUNS = 7                       # independent timing runs per measurement
REPORT = {}
TOTAL_CHECKS = 0

def _count(n):
    global TOTAL_CHECKS
    TOTAL_CHECKS += n

def raises_zero(fn, *args):
    try:
        fn(*args)
    except ZeroDivisionError:
        return True
    return False

# ---------------------------------------------------------------- E1
def e1_exhaustive(B=15, NRANGE=100000, T=3, C=50):
    # (a) homomorphism identities (4), (5) and the normalised product
    regs = [(a, b) for a in range(-B, B + 1) for b in range(-B, B + 1)]
    hom = nrm = 0
    for X in regs:
        vx = value(X)
        N = normalize(X)
        assert N == canonicalize(vx) and N[1] in (0, 1)          # Proposition 2
        nrm += 1
        for Y in regs:
            vy = value(Y)
            assert value(reg_add(X, Y)) == vx + vy
            assert value(reg_sub(X, Y)) == vx - vy
            assert value(reg_mul(X, Y)) == vx * vy
            hom += 3
            assert reg_mul_norm(X, Y) == canonicalize(vx * vy)   # normalised product
            nrm += 1
    # (b) canonical form on a range of integers
    can = 0
    for N in range(-NRANGE, NRANGE + 1):
        K = canonicalize(N)
        assert K[1] in (0, 1) and value(K) == N and normalize(K) == K
        assert K == ((N // 2, 0) if N % 2 == 0 else ((N - 3) // 2, 1))   # (3a), (3b)
        can += 4
    # (c) closed forms of Table 3 for canonical operands
    tab = 0
    for a in range(-C, C + 1):
        for c in range(-C, C + 1):
            assert reg_mul_norm((a, 0), (c, 0)) == (2*a*c, 0)
            assert reg_mul_norm((a, 0), (c, 1)) == (2*a*c + 3*a, 0)
            assert reg_mul_norm((a, 1), (c, 0)) == (2*a*c + 3*c, 0)
            assert reg_mul_norm((a, 1), (c, 1)) == (2*a*c + 3*a + 3*c + 3, 1)
            tab += 4
    # (d) raw pairs form a commutative ring without identity; ker(phi) is an ideal
    small = [(a, b) for a in range(-T, T + 1) for b in range(-T, T + 1)]
    ring = 0
    for X in small:
        for Y in small:
            assert reg_mul(X, Y) == reg_mul(Y, X)
            ring += 1
            for Z in small:
                assert reg_mul(reg_mul(X, Y), Z) == reg_mul(X, reg_mul(Y, Z))
                assert reg_mul(X, reg_add(Y, Z)) == reg_add(reg_mul(X, Y), reg_mul(X, Z))
                ring += 2
    for k in range(-B, B + 1):
        for Y in regs:
            a, b = reg_mul((3 * k, -2 * k), Y)
            assert value((a, b)) == 0 and 2 * a == -3 * b
            ring += 1
    assert all(reg_mul(E, (1, 0)) != (1, 0) for E in regs)       # no identity element
    total = hom + nrm + can + tab + ring
    _count(total)
    REPORT['e1'] = dict(pairs=len(regs) ** 2, hom=hom, nrm=nrm, can=can, tab=tab,
                        ring=ring, total=total)
    print(f"E1  exhaustive |a|,|b|<={B}: {len(regs)**2:,} ordered pairs")
    print(f"      homomorphism identities (4),(5)            {hom:>10,} checks")
    print(f"      register-space normalisation, normalised product {nrm:>4,} checks")
    print(f"      canonical form on [-10^5, 10^5]            {can:>10,} checks")
    print(f"      closed forms of Table 3, |a|,|c|<={C}        {tab:>10,} checks")
    print(f"      ring axioms on raw pairs, kernel is an ideal {ring:>8,} checks")
    print(f"      E1 total                                   {total:>10,} checks")

# ---------------------------------------------------------------- E2
def _reencode(X, rng):
    """Another representative of the same rational: scale numerator and
    denominator by a common non-zero factor and shift every register pair
    by an element of ker(phi)."""
    p, q = as_ratio(X)
    m = rng.choice((-7, -3, -2, -1, 2, 3, 5, 11))
    out = []
    for v in (p * m, q * m):
        a, b = canonicalize(v)
        k = rng.randint(-9, 9)
        out.append((a + 3 * k, b - 2 * k))
    return (out[0], out[1])

def e2_differential(trials=125000, bits=(32, 64, 128, 256)):
    rng = random.Random(SEED + 2)
    def rnd(nb, allow_zero=False):
        if allow_zero and rng.random() < 0.02:
            return 0
        x = rng.getrandbits(nb) or 1
        return -x if rng.getrandbits(1) else x
    c_eager = c_eq = c_lazy = c_repr = c_opt = zero_divs = 0
    for i in range(trials):
        nb = bits[i % len(bits)]
        p, q, r, s = rnd(nb, True), rnd(nb), rnd(nb, True), rnd(nb)
        X, Y = from_ratio(p, q), from_ratio(r, s)
        FX, FY = Fraction(p, q), Fraction(r, s)
        # Algorithm 3: results in normal form, compared WITHOUT further reduction
        S, D, M = frac_add(X, Y), frac_sub(X, Y), frac_mul(X, Y)
        assert as_ratio(S) == (FX + FY).as_integer_ratio() and is_normal(S)
        assert as_ratio(D) == (FX - FY).as_integer_ratio() and is_normal(D)
        assert as_ratio(M) == (FX * FY).as_integer_ratio() and is_normal(M)
        if r != 0:
            Q = frac_div(X, Y)
            assert as_ratio(Q) == (FX / FY).as_integer_ratio() and is_normal(Q)
        else:
            zero_divs += 1
            assert raises_zero(frac_div, X, Y)
        c_eager += 4
        assert frac_eq(X, Y) == (FX == FY)
        c_eq += 1
        # lazy formulas (7a)-(7d): equal values, not equal representations
        for lazy, F in ((frac_add_lazy, FX + FY), (frac_sub_lazy, FX - FY),
                        (frac_mul_lazy, FX * FY)):
            L = lazy(X, Y)
            assert Fraction(*as_ratio(L)) == F
        if r != 0:
            L = frac_div_lazy(X, Y)
            assert Fraction(*as_ratio(L)) == FX / FY
        else:
            assert raises_zero(frac_div_lazy, X, Y)
        c_lazy += 4
        # independence of the representative: other encodings of the same values
        X2, Y2 = _reencode(X, rng), _reencode(Y, rng)
        assert frac_eq(X, X2) and frac_eq(Y2, Y)
        c_eq += 1
        assert frac_add(X2, Y2) == S and frac_sub(X2, Y2) == D and frac_mul(X2, Y2) == M
        if r != 0:
            assert frac_div(X2, Y2) == Q
        else:
            assert raises_zero(frac_div, X2, Y2)
        c_repr += 4
        # classical algorithm behind the register interface: normal-form inputs,
        # output identical (as register tuples) to Algorithm 3
        XN, YN = reduce_fraction(*X), reduce_fraction(*Y)
        assert frac_add_opt(XN, YN) == S and frac_sub_opt(XN, YN) == D
        assert frac_mul_opt(XN, YN) == M
        if r != 0:
            assert frac_div_opt(XN, YN) == Q
        else:
            assert raises_zero(frac_div_opt, XN, YN)
        c_opt += 4
    total = c_eager + c_eq + c_lazy + c_repr + c_opt
    _count(total)
    REPORT['e2'] = dict(trials=trials, eager=c_eager, eq=c_eq, lazy=c_lazy,
                        repr=c_repr, opt=c_opt, total=total, zero_divs=zero_divs)
    print(f"E2  randomized differential vs Fraction: {trials:,} operand sets "
          f"({'/'.join(map(str, bits))}-bit, zero numerators included)")
    print(f"      Algorithm 3, (7a)-(7d), normal-form output  {c_eager:>9,} checks")
    print(f"      equality predicate (6b)                     {c_eq:>9,} checks")
    print(f"      lazy (7a)-(7d), all four operations         {c_lazy:>9,} checks")
    print(f"      re-encoded operands give the same result    {c_repr:>9,} checks")
    print(f"      classical variants on normal-form inputs    {c_opt:>9,} checks")
    print(f"      E2 total                                    {total:>9,} checks "
          f"({zero_divs:,} of the sets have a zero divisor; all rejected)")

def e2_targeted(K=50):
    """Zero divisors, distinct representatives of zero, equal fractions
    with different encodings."""
    one, half = from_ratio(1, 1), from_ratio(1, 2)
    zero = from_ratio(0, 1)
    binary = (frac_add, frac_sub, frac_mul, frac_div,
              frac_add_lazy, frac_sub_lazy, frac_mul_lazy, frac_div_lazy,
              frac_add_opt, frac_sub_opt, frac_mul_opt, frac_div_opt)
    c_den = c_div = c_zero = c_enc = 0
    for k in list(range(-K, 0)) + list(range(1, K + 1)):
        Z = (3 * k, -2 * k)                    # a representative of zero, != (0, 0)
        assert Z != (0, 0) and value(Z) == 0
        # (a) zero-valued denominator: rejected by every entry point
        bad = ((-1, 1), Z)                     # "1 / 0"
        bad0 = ((0, 0), Z)                     # "0 / 0"
        assert raises_zero(reduce_fraction, (-1, 1), Z); c_den += 1
        assert raises_zero(to_ratio, bad); c_den += 1
        for B in (bad, bad0):
            assert raises_zero(frac_eq, B, one) and raises_zero(frac_eq, one, B); c_den += 2
            for op in binary:
                assert raises_zero(op, B, half) and raises_zero(op, half, B); c_den += 2
        # (b) division by a fraction whose numerator is a representative of zero
        ZF = (Z, (-1, 1))                      # 0/1 with a non-canonical zero
        for op in (frac_div, frac_div_lazy, frac_div_opt):
            assert raises_zero(op, half, ZF); c_div += 1
        # (c) distinct representatives of zero behave as zero
        for Qv in (1, -1, 7, -12):
            ZQ = (Z, canonicalize(Qv))
            assert reduce_fraction(*ZQ) == ((0, 0), (-1, 1)); c_zero += 1
            assert frac_eq(ZQ, zero) and not frac_eq(ZQ, one); c_zero += 2
            assert frac_add(half, ZQ) == half and frac_sub(half, ZQ) == half; c_zero += 2
            assert frac_mul(half, ZQ) == ((0, 0), (-1, 1)); c_zero += 1
    # (d) equal fractions with different encodings
    scales = (1, -1, 2, -3, 7)
    shifts = (-2, 0, 3)
    def encodings(p, q):
        for m in scales:
            for k1 in shifts:
                for k2 in shifts:
                    a, b = canonicalize(p * m); c, d = canonicalize(q * m)
                    yield ((a + 3 * k1, b - 2 * k1), (c + 3 * k2, d - 2 * k2))
    grid = [(p, q) for p in range(-6, 7) for q in range(1, 7)]
    for p, q in grid:
        F = Fraction(p, q)
        ref = reduce_fraction(canonicalize(p), canonicalize(q))
        assert as_ratio(ref) == F.as_integer_ratio()
        encs = list(encodings(p, q))
        for E in encs:
            assert reduce_fraction(*E) == ref; c_enc += 1
            assert frac_eq(E, ref) and frac_eq(encs[0], E); c_enc += 2
            assert frac_eq(E, from_ratio(p + 1, q)) is False; c_enc += 1
            assert frac_add(E, half) == frac_add(ref, half); c_enc += 1
            assert frac_mul(E, E) == frac_mul(ref, ref); c_enc += 1
    total = c_den + c_div + c_zero + c_enc
    _count(total)
    REPORT['e2t'] = dict(den=c_den, div=c_div, zero=c_zero, enc=c_enc, total=total)
    print(f"E2t targeted cases")
    print(f"      zero-valued denominators rejected           {c_den:>9,} checks")
    print(f"      zero divisors rejected                      {c_div:>9,} checks")
    print(f"      distinct representatives of zero            {c_zero:>9,} checks")
    print(f"      equal fractions, different encodings        {c_enc:>9,} checks")
    print(f"      E2t total                                   {total:>9,} checks")

# ---------------------------------------------------------------- E3
class Ctr:
    __slots__ = ('MUL', 'ADD', 'DIV', 'MOD')
    def __init__(self): self.MUL = self.ADD = self.DIV = self.MOD = 0

def counted_ops(C):
    def mul(x, y): C.MUL += 1; return x * y
    def add(x, y): C.ADD += 1; return x + y
    def sub(x, y): C.ADD += 1; return x - y
    def div(x, y): C.DIV += 1; return x // y
    def mod(x, y): C.MOD += 1; return x % y
    return mul, add, sub, div, mod

def e3_opcounts():
    C = Ctr(); mul, add, sub, div, mod = counted_ops(C)
    X, Y = (5, 1), (-3, 1)
    # reg_add / reg_sub: literal operator count of formula (4)
    C.__init__(); add(X[0], Y[0]); add(X[1], Y[1]); add_ADD = C.ADD
    # reg_mul: literal operator count of formula (5)
    #   a' = 2*(a1*a2)                         -> 2 MUL
    #   b' = 2*(a1*b2 + b1*a2) + 3*(b1*b2)     -> 5 MUL + 2 ADD
    C.__init__()
    ap = mul(2, mul(X[0], Y[0]))
    bp = add(mul(2, add(mul(X[0], Y[1]), mul(X[1], Y[0]))), mul(3, mul(X[1], Y[1])))
    assert (ap, bp) == reg_mul(X, Y)
    mul_MUL, mul_ADD = C.MUL, C.ADD          # 7 MUL (4 full-precision) + 2 ADD
    # value: 2a + 3b
    C.__init__(); v = add(mul(2, X[0]), mul(3, X[1])); assert v == value(X)
    val_MUL, val_ADD = C.MUL, C.ADD
    # canonicalize: N = 2q + r, then q - r
    C.__init__(); q = div(7, 2); r = mod(7, 2); k = (sub(q, r), r)
    assert k == canonicalize(7)
    can = (C.DIV, C.MOD, C.ADD)              # 1 DIV + 1 MOD + 1 SUB
    # normalize: (a, b) -> (a + 3*floor(b/2), b mod 2)
    C.__init__(); R = (-6, 7)
    k = div(R[1], 2); r = mod(R[1], 2); n = (add(R[0], mul(3, k)), r)
    assert n == normalize(R)
    nor = (C.DIV, C.MOD, C.MUL, C.ADD)       # 1 DIV + 1 MOD + 1 MUL + 1 ADD
    # gcd division count, random 64-bit pairs
    rng = random.Random(SEED + 3)
    tot = n = 0
    for _ in range(10000):
        x, y = rng.getrandbits(64), rng.getrandbits(64)
        k = 0
        while y: x, y = y, x % y; k += 1
        tot += k; n += 1
    gcd_avg = tot / n
    REPORT['e3'] = dict(add_ADD=add_ADD, mul_MUL=mul_MUL, mul_ADD=mul_ADD,
                        val=(val_MUL, val_ADD), can=can, nor=nor, gcd_avg=gcd_avg)
    print(f"E3  literal operator counts (unit-cost model): reg_add/sub = {add_ADD} ADD; "
          f"reg_mul = {mul_MUL} MUL + {mul_ADD} ADD (4 of the MUL are full-precision products); "
          f"value = {val_MUL} MUL + {val_ADD} ADD; "
          f"canonicalize = {can[0]} DIV + {can[1]} MOD + {can[2]} SUB; "
          f"normalize = {nor[0]} DIV + {nor[1]} MOD + {nor[2]} MUL + {nor[3]} ADD; "
          f"gcd on 64-bit operands: {gcd_avg:.1f} division steps on average")

# ---------------------------------------------------------------- E4
def e4_growth(ns=(16, 64, 256, 1024)):
    rows = []
    for n in ns:
        # lazy: H_n via formulas (7a) only
        L = from_ratio(1, 1)
        for k in range(2, n + 1):
            L = frac_add_lazy(L, from_ratio(1, k))
        lazy_bits = value(L[1]).bit_length()
        # Algorithm 3 (reduce after every op)
        A = from_ratio(1, 1)
        for k in range(2, n + 1):
            A = frac_add(A, from_ratio(1, k))
        red_bits = value(A[1]).bit_length()
        # cross-check against Fraction
        F = sum((Fraction(1, k) for k in range(2, n + 1)), Fraction(1, 1))
        assert as_ratio(A) == F.as_integer_ratio()
        assert Fraction(*as_ratio(L)) == F
        rows.append((n, lazy_bits, red_bits))
        print(f"E4  H_{n}: lazy denominator {lazy_bits:,} bits; "
              f"Algorithm 3 (reduced) {red_bits:,} bits; ratio {lazy_bits/red_bits:.1f}x")
    REPORT['e4'] = rows

# ---------------------------------------------------------------- E5
def bench(fn, args, min_time=0.1, repeat=3):
    """Seconds per call of fn(*a) over the argument list: the loop count is
    calibrated to last at least min_time, the loop is timed `repeat` times
    and the fastest timing is returned (as timeit recommends; slower
    repetitions measure interference from other processes)."""
    n = len(args)
    reps = 1
    while True:
        t0 = time.perf_counter()
        for _ in range(reps):
            for a in args:
                fn(*a)
        dt = time.perf_counter() - t0
        if dt >= min_time:
            break
        reps *= 2
    best = dt
    for _ in range(repeat - 1):
        t0 = time.perf_counter()
        for _ in range(reps):
            for a in args:
                fn(*a)
        best = min(best, time.perf_counter() - t0)
    return best / (reps * n)

def med_range(xs):
    return statistics.median(xs), min(xs), max(xs)

def cpu_model():
    try:
        with open('/proc/cpuinfo') as f:
            for line in f:
                if line.lower().startswith('model name'):
                    return line.split(':', 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or platform.machine() or 'unknown'

def _rand_width(rng, nb, parity=None):
    """Random integer of exactly nb bits with the requested parity."""
    x = rng.getrandbits(nb) | (1 << (nb - 1))
    if parity == 0:
        x &= ~1
    elif parity == 1:
        x |= 1
    return x

def e5a_products(widths=(64, 1024, 16384), npairs=200):
    """Table 10: isolated products, canonical inputs by parity class."""
    rng = random.Random(SEED + 5)
    rows = []
    print("E5a isolated products: median ratio to native multiplication "
          f"[min-max over {RUNS} runs]")
    for nb in widths:
        ints = {}
        for name, (px, py) in (('ee', (0, 0)), ('eo', (0, 1)), ('oo', (1, 1))):
            ints[name] = [(_rand_width(rng, nb, px), _rand_width(rng, nb, py))
                          for _ in range(npairs)]
        can = {k: [(canonicalize(x), canonicalize(y)) for x, y in v] for k, v in ints.items()}
        mixed_i = ints['ee'][::3] + ints['eo'][::3] + ints['oo'][::3]
        mixed = [(canonicalize(x), canonicalize(y)) for x, y in mixed_i]
        bal = [((rng.getrandbits(nb), rng.getrandbits(nb)),
                (rng.getrandbits(nb), rng.getrandbits(nb))) for _ in range(npairs)]
        nat_us, ratios = [], {k: [] for k in ('ee', 'eo', 'oo', 'norm', 'bal')}
        for _ in range(RUNS):
            t_nat = bench(lambda x, y: x * y, mixed_i)
            nat_us.append(t_nat * 1e6)
            for k in ('ee', 'eo', 'oo'):
                ratios[k].append(bench(reg_mul, can[k]) / t_nat)
            ratios['norm'].append(bench(reg_mul_norm, mixed) / t_nat)
            ratios['bal'].append(bench(reg_mul, bal) / t_nat)
        row = dict(bits=nb, native=med_range(nat_us),
                   **{k: med_range(v) for k, v in ratios.items()})
        rows.append(row)
        f = lambda t: f"{t[0]:5.2f} [{t[1]:.2f}-{t[2]:.2f}]"
        print(f"    {nb:>6}-bit: native {row['native'][0]:9.3f} us "
              f"[{row['native'][1]:.3f}-{row['native'][2]:.3f}]")
        print(f"        raw reg_mul, canonical inputs  b=(0,0) {f(row['ee'])}   "
              f"b=(0,1) {f(row['eo'])}   b=(1,1) {f(row['oo'])}")
        print(f"        reg_mul + normalize (mixed parity)     {f(row['norm'])}")
        print(f"        raw reg_mul, balanced non-canonical    {f(row['bal'])}")
    REPORT['e5a'] = rows

def e5b_chains():
    """Table 11: chains of products under three policies —
    native integers, retained raw outputs, normalisation after every product."""
    rng = random.Random(SEED + 6)
    print(f"E5b product chains: per-product time, median [min-max] over {RUNS} runs")
    rows = []

    # W1: running product of L fresh 64-bit factors
    L = 512
    fac_i = [_rand_width(rng, 64) for _ in range(L)]
    fac_r = [canonicalize(x) for x in fac_i]
    def w1_native():
        acc = 1
        for x in fac_i: acc = acc * x
        return acc
    def w1_raw():
        acc = (-1, 1)
        for X in fac_r: acc = reg_mul(acc, X)
        return acc
    def w1_norm():
        acc = (-1, 1)
        for X in fac_r: acc = reg_mul_norm(acc, X)
        return acc
    # W2: successive squarings of a 64-bit odd value
    S2 = 10
    x0 = _rand_width(rng, 64, 1)
    def w2_native():
        acc = x0
        for _ in range(S2): acc = acc * acc
        return acc
    def w2_raw():
        acc = canonicalize(x0)
        for _ in range(S2): acc = reg_mul(acc, acc)
        return acc
    def w2_norm():
        acc = canonicalize(x0)
        for _ in range(S2): acc = reg_mul_norm(acc, acc)
        return acc
    # W3: successive squarings of kappa(1) = (-1, 1); the value stays 1
    S3 = 16
    def w3_native():
        acc = 1
        for _ in range(S3): acc = acc * acc
        return acc
    def w3_raw():
        acc = (-1, 1)
        for _ in range(S3): acc = reg_mul(acc, acc)
        return acc
    def w3_norm():
        acc = (-1, 1)
        for _ in range(S3): acc = reg_mul_norm(acc, acc)
        return acc

    for name, nops, nat, raw, nrm in (
            (f"W1 running product of {L} 64-bit factors", L, w1_native, w1_raw, w1_norm),
            (f"W2 {S2} successive squarings of a 64-bit value", S2, w2_native, w2_raw, w2_norm),
            (f"W3 {S3} successive squarings of kappa(1)", S3, w3_native, w3_raw, w3_norm)):
        v, R, N = nat(), raw(), nrm()
        assert value(R) == v and value(N) == v and N == canonicalize(v)
        t = {k: [] for k in ('nat', 'raw', 'nrm')}
        for _ in range(RUNS):
            t['nat'].append(bench(nat, [()]) / nops * 1e6)
            t['raw'].append(bench(raw, [()]) / nops * 1e6)
            t['nrm'].append(bench(nrm, [()]) / nops * 1e6)
        row = dict(name=name, nops=nops, value_bits=v.bit_length(),
                   raw_bits=(abs(R[0]).bit_length(), abs(R[1]).bit_length()),
                   nrm_bits=(abs(N[0]).bit_length(), abs(N[1]).bit_length()),
                   nat=med_range(t['nat']), raw=med_range(t['raw']), nrm=med_range(t['nrm']))
        rows.append(row)
        g = lambda x: f"{x[0]:10.3f} us [{x[1]:.3f}-{x[2]:.3f}]"
        print(f"    {name}: value {row['value_bits']:,} bits")
        print(f"        native integers        {g(row['nat'])}")
        print(f"        retained raw outputs   {g(row['raw'])}  ({row['raw'][0]/row['nat'][0]:.2f}x)  "
              f"final registers {row['raw_bits'][0]:,} / {row['raw_bits'][1]:,} bits")
        print(f"        normalise each product {g(row['nrm'])}  ({row['nrm'][0]/row['nat'][0]:.2f}x)  "
              f"final registers {row['nrm_bits'][0]:,} / {row['nrm_bits'][1]:,} bits")
    # deterministic register growth of the squaring chain of kappa(1)
    acc, growth = (-1, 1), []
    for k in range(1, 21):
        acc = reg_mul(acc, acc)
        assert value(acc) == 1 and acc[0] == 2 ** (2 ** k - 1) and 3 * acc[1] == 1 - 2 * acc[0]
        growth.append((k, acc[0].bit_length(), abs(acc[1]).bit_length()))
    print("    kappa(1) squared k times with raw reg_mul (value stays 1), "
          "first-register bits: " + ", ".join(f"k={k}: {a:,}" for k, a, _ in growth if k in (1, 2, 3, 5, 10, 20)))
    print("    first three pairs:", reg_mul((-1, 1), (-1, 1)),
          reg_mul((2, -1), (2, -1)), reg_mul((8, -5), (8, -5)))
    REPORT['e5b'] = rows
    REPORT['e5b_growth'] = growth

def e5c_rational_chain(nops=2000):
    """Table 12: chain of mixed rational operations under matched invariants:
    every engine builds its operands in reduced form and keeps every
    intermediate result reduced with a positive denominator."""
    rng = random.Random(SEED + 7)
    pool = [(rng.getrandbits(64) + 1, rng.getrandbits(64) + 1) for _ in range(nops)]
    kinds = [i % 4 for i in range(len(pool))]
    def chain(engine_ops, mk):
        acc = mk(3, 7)
        for k, (p, q) in zip(kinds, pool):
            acc = engine_ops[k](acc, mk(p, q))
        return acc
    ra_verb = [frac_add, frac_mul, frac_sub, frac_div]
    ra_opt  = [frac_add_opt, frac_mul_opt, frac_sub_opt, frac_div_opt]
    fr_ops  = [Fraction.__add__, Fraction.__mul__, Fraction.__sub__, Fraction.__truediv__]
    r_verb = chain(ra_verb, from_ratio_reduced)
    r_opt  = chain(ra_opt,  from_ratio_reduced)
    r_fr   = chain(fr_ops, Fraction)
    # identical reduced representations, compared without any final reduction
    assert r_verb == r_opt and is_normal(r_verb)
    assert as_ratio(r_verb) == (r_fr.numerator, r_fr.denominator)
    den_bits = r_fr.denominator.bit_length()
    t = {k: [] for k in ('verb', 'opt', 'fr')}
    for _ in range(RUNS):
        t['fr'].append(bench(lambda: chain(fr_ops, Fraction), [()], min_time=0.3) / nops * 1e6)
        t['opt'].append(bench(lambda: chain(ra_opt, from_ratio_reduced), [()], min_time=0.3) / nops * 1e6)
        t['verb'].append(bench(lambda: chain(ra_verb, from_ratio_reduced), [()], repeat=1) / nops * 1e6)
    ratio_v = [a / b for a, b in zip(t['verb'], t['fr'])]
    ratio_o = [a / b for a, b in zip(t['opt'], t['fr'])]
    row = dict(nops=nops, den_bits=den_bits, verb=med_range(t['verb']), opt=med_range(t['opt']),
               fr=med_range(t['fr']), ratio_v=med_range(ratio_v), ratio_o=med_range(ratio_o))
    g = lambda x: f"{x[0]:9.2f} us [{x[1]:.2f}-{x[2]:.2f}]"
    h = lambda x: f"{x[0]:6.2f}x [{x[1]:.2f}-{x[2]:.2f}]"
    print(f"E5c rational chain ({nops} mixed ops, reduced operands, final denominator "
          f"{den_bits:,} bits), per operation, median [min-max] over {RUNS} runs:")
    print(f"      Algorithm 3, formulas (7a)-(7d) + reduce   {g(row['verb'])}   {h(row['ratio_v'])} Fraction")
    print(f"      classical algorithm, register interface    {g(row['opt'])}   {h(row['ratio_o'])} Fraction")
    print(f"      fractions.Fraction                         {g(row['fr'])}")
    REPORT['e5c'] = row

def e5_timing():
    print(f"E5  CPython {sys.version.split()[0]}, {cpu_model()}, {platform.system()} "
          f"{platform.machine()}, single core, {RUNS} independent runs per measurement")
    e5a_products()
    e5b_chains()
    e5c_rational_chain()

# ----------------------------------------------------------------
if __name__ == "__main__":
    t0 = time.perf_counter()
    phases = dict(e1=e1_exhaustive, e2=e2_differential, e2t=e2_targeted,
                  e3=e3_opcounts, e4=e4_growth, e5=e5_timing)
    todo = [a.lower() for a in sys.argv[1:]] or list(phases)
    unknown = [name for name in todo if name not in phases]
    if unknown:
        sys.exit("unknown experiment: %s\nvalid names are: %s"
                 % (", ".join(unknown), ", ".join(phases)))
    for name in todo:
        phases[name]()
    if TOTAL_CHECKS:
        print(f"TOTAL individual checks: {TOTAL_CHECKS:,}")
    print(f"ALL CHECKS PASSED in {time.perf_counter()-t0:.1f}s")
