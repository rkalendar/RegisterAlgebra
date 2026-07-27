"""reproduce_results.py — reproduces every number in Section 11 of the paper.

Usage:  python3 reproduce_results.py          (runs in ~1-2 minutes)

Experiments
  E1  exhaustive homomorphism check on small registers
  E2  randomized differential testing against fractions.Fraction
  E3  integer-operation counts (empirical confirmation of Table 5)
  E4  register growth: lazy formulas (7a)-(7d) vs Algorithm 3
  E5  microbenchmarks: register vs native operations

Supplementary Material S1. Pure Python, no dependencies.
"""
import random, sys, time
from fractions import Fraction

from rau import (value, canonicalize, reg_add, reg_sub, reg_mul,
                 normalize, gcd, reduce_fraction,
                 frac_add, frac_sub, frac_mul, frac_div, frac_eq,
                 frac_add_lazy, frac_mul_lazy, from_ratio, to_ratio)

random.seed(20260713)
REPORT = {}

# ---------------------------------------------------------------- E1
def e1_exhaustive(B=15, NRANGE=100000):
    checks = 0
    regs = [(a, b) for a in range(-B, B + 1) for b in range(-B, B + 1)]
    for X in regs:
        vx = value(X)
        for Y in regs:
            vy = value(Y)
            assert value(reg_add(X, Y)) == vx + vy
            assert value(reg_sub(X, Y)) == vx - vy
            assert value(reg_mul(X, Y)) == vx * vy
            checks += 3
    for N in range(-NRANGE, NRANGE + 1):
        K = canonicalize(N)
        assert K[1] in (0, 1) and value(K) == N and normalize(K) == K
        checks += 3
    REPORT['e1_pairs'] = len(regs) ** 2
    REPORT['e1_checks'] = checks
    print(f"E1  exhaustive |a|,|b|<={B}: {len(regs)**2:,} pairs, "
          f"{checks:,} identities hold; kappa canonical on [-10^5,10^5]")

# ---------------------------------------------------------------- E2
def e2_differential(trials=125000, bits=(32, 64, 128, 256)):
    def rnd(nb):
        x = random.getrandbits(nb) or 1
        return -x if random.getrandbits(1) else x
    ops = 0
    for i in range(trials):
        nb = bits[i % len(bits)]
        p, q, r, s = rnd(nb), rnd(nb), rnd(nb), rnd(nb)
        X, Y = from_ratio(p, q), from_ratio(r, s)
        FX, FY = Fraction(p, q), Fraction(r, s)
        assert to_ratio(frac_add(X, Y)) == (FX + FY).as_integer_ratio()
        assert to_ratio(frac_sub(X, Y)) == (FX - FY).as_integer_ratio()
        assert to_ratio(frac_mul(X, Y)) == (FX * FY).as_integer_ratio()
        assert to_ratio(frac_div(X, Y)) == (FX / FY).as_integer_ratio()
        assert frac_eq(X, Y) == (FX == FY)
        ops += 5
        # lazy formulas agree semantically as well
        L = frac_add_lazy(X, Y)
        assert Fraction(value(L[0]), value(L[1])) == FX + FY
        M = frac_mul_lazy(X, Y)
        assert Fraction(value(M[0]), value(M[1])) == FX * FY
        ops += 2
    REPORT['e2_trials'] = trials
    REPORT['e2_ops'] = ops
    print(f"E2  randomized differential vs Fraction: {trials:,} operand sets "
          f"({'/'.join(map(str,bits))}-bit), {ops:,} operation checks, 0 discrepancies")

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
    # canonicalize: one parity test, one exact division
    C.__init__(); mod(7, 2); div(7 - 3, 2); can_DIV, can_MOD = C.DIV, C.MOD
    # gcd division count, random 64-bit pairs
    tot = n = 0
    for _ in range(10000):
        x, y = random.getrandbits(64), random.getrandbits(64)
        k = 0
        while y: x, y = y, x % y; k += 1
        tot += k; n += 1
    gcd_avg = tot / n
    REPORT['e3'] = dict(add_ADD=add_ADD, mul_MUL=mul_MUL, mul_ADD=mul_ADD,
                        val=(val_MUL, val_ADD), can=(can_DIV, can_MOD), gcd_avg=gcd_avg)
    print(f"E3  literal operator counts: reg_add/sub = {add_ADD} ADD; "
          f"reg_mul = {mul_MUL} MUL + {mul_ADD} ADD (4 of the MUL are full-precision products); "
          f"value = {val_MUL} MUL + {val_ADD} ADD; canonicalize = {can_DIV} DIV + {can_MOD} MOD; "
          f"gcd on 64-bit operands: {gcd_avg:.1f} divisions on average")

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
        assert to_ratio(A) == F.as_integer_ratio()
        assert Fraction(value(L[0]), value(L[1])) == F
        rows.append((n, lazy_bits, red_bits))
        print(f"E4  H_{n}: lazy denominator {lazy_bits:,} bits; "
              f"Algorithm 3 (reduced) {red_bits:,} bits; ratio {lazy_bits/red_bits:.1f}x")
    REPORT['e4'] = rows

# ---------------------------------------------------------------- E5
def bench(fn, args_iter, min_time=0.4):
    args = list(args_iter)
    n = len(args)
    reps = 1
    while True:
        t0 = time.perf_counter()
        for _ in range(reps):
            for a in args:
                fn(*a)
        dt = time.perf_counter() - t0
        if dt >= min_time:
            return dt / (reps * n)
        reps *= 2

def e5_timing():
    from rau import frac_add_opt, frac_sub_opt, frac_mul_opt, frac_div_opt
    print(f"E5  CPython {sys.version.split()[0]}, single core")
    rows = []
    # (a) integer layer: reg_mul vs native multiplication
    for nb in (64, 1024, 16384):
        xs = [(random.getrandbits(nb) | 1, random.getrandbits(nb) | 1) for _ in range(200)]
        can = [(canonicalize(x), canonicalize(y)) for x, y in xs]          # b in {0,1}
        bal = [((random.getrandbits(nb), random.getrandbits(nb)),
                (random.getrandbits(nb), random.getrandbits(nb))) for _ in range(200)]
        t_nat = bench(lambda x, y: x * y, xs)
        t_can = bench(reg_mul, can)
        t_bal = bench(reg_mul, bal)
        rows.append((nb, t_nat*1e6, t_can*1e6, t_bal*1e6))
        print(f"    reg_mul {nb:>6}-bit: canonical {t_can*1e6:9.3f} us ({t_can/t_nat:4.2f}x), "
              f"balanced {t_bal*1e6:9.3f} us ({t_bal/t_nat:4.2f}x)  vs native {t_nat*1e6:9.3f} us")
    REPORT['e5a'] = rows
    # (b) rational chain: 2000 mixed ops, 64-bit operands, three engines
    pool = [(random.getrandbits(64) + 1, random.getrandbits(64) + 1) for _ in range(2000)]
    kinds = [i % 4 for i in range(len(pool))]
    def chain(engine_ops, mk, fin):
        acc = mk(3, 7)
        for k, (p, q) in zip(kinds, pool):
            acc = engine_ops[k](acc, mk(p, q))
        return fin(acc)
    ra_verb = [frac_add, frac_mul, frac_sub, frac_div]
    ra_opt  = [frac_add_opt, frac_mul_opt, frac_sub_opt, frac_div_opt]
    fr_ops  = [Fraction.__add__, Fraction.__mul__, Fraction.__sub__, Fraction.__truediv__]
    r_verb = chain(ra_verb, from_ratio, to_ratio)
    r_opt  = chain(ra_opt,  from_ratio, to_ratio)
    r_fr   = chain(fr_ops, Fraction, lambda f: f.as_integer_ratio())
    assert r_verb == r_opt == r_fr
    den_bits = r_fr[1].bit_length()
    t_v = bench(lambda: chain(ra_verb, from_ratio, to_ratio), [()], min_time=1.0)
    t_o = bench(lambda: chain(ra_opt,  from_ratio, to_ratio), [()], min_time=1.0)
    t_f = bench(lambda: chain(fr_ops, Fraction, lambda f: f.as_integer_ratio()), [()], min_time=1.0)
    n = len(pool)
    print(f"    rational chain ({n} ops, final denominator {den_bits:,} bits), per op:")
    print(f"      Algorithm 3 verbatim {t_v/n*1e6:9.2f} us   ({t_v/t_f:5.2f}x Fraction)")
    print(f"      optimised (Knuth)    {t_o/n*1e6:9.2f} us   ({t_o/t_f:5.2f}x Fraction)")
    print(f"      fractions.Fraction   {t_f/n*1e6:9.2f} us")
    REPORT['e5b'] = (n, den_bits, t_v/n*1e6, t_o/n*1e6, t_f/n*1e6)

# ----------------------------------------------------------------
if __name__ == "__main__":
    t0 = time.perf_counter()
    phases = dict(e1=e1_exhaustive, e2=e2_differential, e3=e3_opcounts,
                  e4=e4_growth, e5=e5_timing)
    todo = sys.argv[1:] or list(phases)
    for name in todo:
        phases[name]()
    print(f"ALL CHECKS PASSED in {time.perf_counter()-t0:.1f}s")
