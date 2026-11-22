"""bench_bits.py — what do shifts and masks gain over the reference operations?

Usage:  python3 bench_bits.py

For 64-, 1,024- and 16,384-bit operands it times the reference operations of
rau.py against the shift-and-mask variants of rau_bits.py and against native
integer multiplication.  Each measurement is the fastest of five timed loops,
repeated in RUNS independent runs; the median and the range are printed.
Results are wall-clock times and depend on the machine.
"""
import platform, random, statistics, sys, time

from rau import canonicalize, normalize, reg_mul, reg_mul_norm
from rau_bits import (canonicalize_bits, normalize_bits, reg_mul_bits,
                      reg_mul_canon, reg_mul_via_value)

SEED, RUNS = 20260713, 7

def bench(fn, args, min_time=0.1, repeat=5):
    n, reps = len(args), 1
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

def cpu_model():
    try:
        with open('/proc/cpuinfo') as f:
            for line in f:
                if line.lower().startswith('model name'):
                    return line.split(':', 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or platform.machine() or 'unknown'

def rand_width(rng, nb, parity):
    x = rng.getrandbits(nb) | (1 << (nb - 1))
    return (x & ~1) if parity == 0 else (x | 1)

def med(xs):
    return "%6.2f [%.2f-%.2f]" % (statistics.median(xs), min(xs), max(xs))

def main(widths=(64, 1024, 16384), npairs=198):
    rng = random.Random(SEED + 11)
    print("CPython %s, %s, %s %s, %d independent runs" % (
        sys.version.split()[0], cpu_model(), platform.system(), platform.machine(), RUNS))
    for nb in widths:
        ints = [(rand_width(rng, nb, px), rand_width(rng, nb, py))
                for px, py in ((0, 0), (0, 1), (1, 1)) for _ in range(npairs // 3)]
        can = [(canonicalize(x), canonicalize(y)) for x, y in ints]
        bal = [((rng.getrandbits(nb), rng.getrandbits(nb)),
                (rng.getrandbits(nb), rng.getrandbits(nb))) for _ in range(npairs)]
        raw = [(reg_mul(X, Y),) for X, Y in can]              # inputs for normalisation
        one = [(x,) for x, _ in ints]
        # all variants must agree before anything is timed
        for (X, Y), (x, y) in zip(can, ints):
            assert reg_mul_canon(X, Y) == reg_mul_norm(X, Y) == reg_mul_via_value(X, Y) == canonicalize(x * y)
            assert reg_mul_bits(X, Y) == reg_mul(X, Y)
        for U, V in bal:
            assert reg_mul_bits(U, V) == reg_mul(U, V)
        cases = [
            ("canonical in, raw out",       "reg_mul          (reference)",        reg_mul,           can),
            ("",                            "reg_mul_bits     (shifts)",           reg_mul_bits,      can),
            ("canonical in, canonical out", "reg_mul_norm     (reference)",        reg_mul_norm,      can),
            ("",                            "reg_mul_canon    (shifts and masks)", reg_mul_canon,     can),
            ("",                            "reg_mul_via_value (decode, multiply, encode)", reg_mul_via_value, can),
            ("balanced non-canonical, raw", "reg_mul          (reference)",        reg_mul,           bal),
            ("",                            "reg_mul_bits     (shifts)",           reg_mul_bits,      bal),
        ]
        ratios = [[] for _ in cases]
        nat_us, small = [], {k: [] for k in ("canonicalize", "canonicalize_bits", "normalize", "normalize_bits")}
        for _ in range(RUNS):
            t_nat = bench(lambda x, y: x * y, ints)
            nat_us.append(t_nat * 1e6)
            for i, c in enumerate(cases):
                ratios[i].append(bench(c[2], c[3]) / t_nat)
            small["canonicalize"].append(bench(canonicalize, one) * 1e6)
            small["canonicalize_bits"].append(bench(canonicalize_bits, one) * 1e6)
            small["normalize"].append(bench(normalize, raw) * 1e6)
            small["normalize_bits"].append(bench(normalize_bits, raw) * 1e6)
        print("\n%d-bit operands: native multiplication %.3f us [%.3f-%.3f]" % (
            nb, statistics.median(nat_us), min(nat_us), max(nat_us)))
        print("  register product, time relative to native multiplication: median [min-max]")
        for i, c in enumerate(cases):
            if c[0]:
                print("    " + c[0])
            print("      %-46s %s" % (c[1], med(ratios[i])))
        print("  linear-time steps, us: median [min-max]")
        for k, v in small.items():
            print("      %-46s %s" % (k, med(v)))

if __name__ == "__main__":
    main()
