/* rau_word.c — Register Algebra with registers held in one machine word.
 *
 * Question: how close to one native instruction can the register operations
 * get when every register is a 64-bit word and the constant scalings, the
 * halving and the products by b, d in {0, 1} are done with shifts and masks?
 *
 * Operand values are kept below 2^30 in magnitude, so that every product and
 * every intermediate register fits in a signed 64-bit word.  This is the
 * price of the word-size setting: the range is bounded, and a chain of raw
 * (unnormalised) products overflows after a few steps (Section 5.3 of the
 * paper), so only normalised chains are timed.
 *
 *   gcc -O2 -fno-tree-vectorize -o rau_word rau_word.c
 *   ./rau_word            benchmark
 *   ./rau_word --dump     print test vectors for cross-checking with rau.py
 *
 * Two figures are reported for each operation:
 *   throughput — independent operations on arrays of operands;
 *   latency    — a dependent chain, each result feeding the next operation.
 * Each figure is the fastest of RUNS runs (the estimate least affected by
 * other processes); the median is printed beside it.  The time-stamp counter
 * ticks at the nominal clock rate, which is not the actual core frequency, so
 * core cycles are estimated by calibration: a native 64-bit multiply has a
 * latency of 3 cycles on current x86-64 cores.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <x86intrin.h>

typedef struct { int64_t a, b; } reg;          /* value 2a + 3b */

#define NOINLINE __attribute__((noinline))

/* ---- reference: the formulas as written in the paper ------------------ */
static inline int64_t value(reg r)            { return 2 * r.a + 3 * r.b; }
static inline reg mul_raw(reg x, reg y)        /* formula (5), 4 products */
{ return (reg){ 2 * x.a * y.a, 2 * (x.a * y.b + x.b * y.a) + 3 * x.b * y.b }; }

/* ---- shifts and masks -------------------------------------------------- */
static inline reg canon(int64_t n)             /* (3c): q = n >> 1, r = n & 1 */
{ int64_t r = n & 1; return (reg){ (n >> 1) - r, r }; }
static inline reg normalize(reg r)             /* (3d) */
{ int64_t k = r.b >> 1; return (reg){ r.a + 3 * k, r.b & 1 }; }

static inline reg mul_norm_ref(reg x, reg y)   /* formula (5), then (3d) */
{ return normalize(mul_raw(x, y)); }

static inline reg mul_canon(reg x, reg y)      /* Table 3: canonical operands, 1 product */
{
    int64_t t = (x.a & -y.b) + (y.a & -x.b) + (x.b & y.b);   /* a*d + c*b + b*d */
    return (reg){ 2 * (x.a * y.a) + 3 * t, x.b & y.b };
}
static inline reg mul_via_value(reg x, reg y)  /* decode, one product, encode */
{ return canon(value(x) * value(y)); }

static inline reg add_norm(reg x, reg y)       /* (4), then (3d), canonical operands */
{ int64_t s = x.b + y.b; return (reg){ x.a + y.a + 3 * (s >> 1), s & 1 }; }

/* ---- data --------------------------------------------------------------- */
#define N 1024                                  /* operands stay in the L1 cache */
static int64_t XI[N], YI[N], UI[N], VI[N];
static reg     XR[N], YR[N], UR[N], VR[N];

static uint64_t rng_state = 20260713u;
static uint64_t rnd(void)                        /* xorshift64* */
{
    rng_state ^= rng_state >> 12; rng_state ^= rng_state << 25; rng_state ^= rng_state >> 27;
    return rng_state * 2685821657736338717ull;
}
static int64_t rnd_signed(int bits)
{
    int64_t v = (int64_t)(rnd() >> (64 - bits));
    return (rnd() & 1) ? -v : v;
}
#define BARRIER() __asm__ volatile("" ::: "memory")   /* arrays may have changed */

/* ---- throughput: independent operations --------------------------------- */
NOINLINE static uint64_t tp_native_mul(long reps)
{ uint64_t s = 0; for (long r = 0; r < reps; r++) { BARRIER();
    for (int i = 0; i < N; i++) s += (uint64_t)(XI[i] * YI[i]); } return s; }
NOINLINE static uint64_t tp_mul_raw(long reps)
{ uint64_t s = 0; for (long r = 0; r < reps; r++) { BARRIER();
    for (int i = 0; i < N; i++) { reg z = mul_raw(XR[i], YR[i]); s += (uint64_t)(z.a + z.b); } } return s; }
NOINLINE static uint64_t tp_mul_norm_ref(long reps)
{ uint64_t s = 0; for (long r = 0; r < reps; r++) { BARRIER();
    for (int i = 0; i < N; i++) { reg z = mul_norm_ref(XR[i], YR[i]); s += (uint64_t)(z.a + z.b); } } return s; }
NOINLINE static uint64_t tp_mul_canon(long reps)
{ uint64_t s = 0; for (long r = 0; r < reps; r++) { BARRIER();
    for (int i = 0; i < N; i++) { reg z = mul_canon(XR[i], YR[i]); s += (uint64_t)(z.a + z.b); } } return s; }
NOINLINE static uint64_t tp_mul_via_value(long reps)
{ uint64_t s = 0; for (long r = 0; r < reps; r++) { BARRIER();
    for (int i = 0; i < N; i++) { reg z = mul_via_value(XR[i], YR[i]); s += (uint64_t)(z.a + z.b); } } return s; }
NOINLINE static uint64_t tp_native_add(long reps)
{ uint64_t s = 0; for (long r = 0; r < reps; r++) { BARRIER();
    for (int i = 0; i < N; i++) s += (uint64_t)(XI[i] + YI[i]); } return s; }
NOINLINE static uint64_t tp_add_norm(long reps)
{ uint64_t s = 0; for (long r = 0; r < reps; r++) { BARRIER();
    for (int i = 0; i < N; i++) { reg z = add_norm(XR[i], YR[i]); s += (uint64_t)(z.a + z.b); } } return s; }

/* ---- latency: each result feeds the next operation ----------------------
 * The factors are +1 or -1, so the value stays +-1 and nothing overflows;
 * the multiply instruction takes the same time whatever its operands.     */
NOINLINE static uint64_t lat_native_mul(long reps)
{ int64_t acc = 1; for (long r = 0; r < reps; r++) { BARRIER();
    for (int i = 0; i < N; i++) acc = acc * UI[i]; } return (uint64_t)acc; }
NOINLINE static uint64_t lat_mul_norm_ref(long reps)
{ reg acc = canon(1); for (long r = 0; r < reps; r++) { BARRIER();
    for (int i = 0; i < N; i++) acc = mul_norm_ref(acc, UR[i]); } return (uint64_t)value(acc); }
NOINLINE static uint64_t lat_mul_canon(long reps)
{ reg acc = canon(1); for (long r = 0; r < reps; r++) { BARRIER();
    for (int i = 0; i < N; i++) acc = mul_canon(acc, UR[i]); } return (uint64_t)value(acc); }
NOINLINE static uint64_t lat_mul_via_value(long reps)
{ reg acc = canon(1); for (long r = 0; r < reps; r++) { BARRIER();
    for (int i = 0; i < N; i++) acc = mul_via_value(acc, UR[i]); } return (uint64_t)value(acc); }
NOINLINE static uint64_t lat_native_add(long reps)
{ int64_t acc = 0; for (long r = 0; r < reps; r++) { BARRIER();
    for (int i = 0; i < N; i++) acc = acc + VI[i]; } return (uint64_t)acc; }
NOINLINE static uint64_t lat_add_norm(long reps)
{ reg acc = canon(0); for (long r = 0; r < reps; r++) { BARRIER();
    for (int i = 0; i < N; i++) acc = add_norm(acc, VR[i]); } return (uint64_t)value(acc); }

/* ---- harness -------------------------------------------------------------- */
#define RUNS 15
static double now_ns(void)
{ struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec * 1e9 + ts.tv_nsec; }
static int cmp_d(const void *p, const void *q)
{ double a = *(const double *)p, b = *(const double *)q; return (a > b) - (a < b); }

typedef struct { double ns, tsc, ns_min, tsc_min; uint64_t sum; } result;
static result measure(uint64_t (*fn)(long), long reps)
{
    double ns[RUNS], tk[RUNS]; result r; r.sum = 0;
    fn(reps / 8 + 1);                                   /* warm-up */
    for (int k = 0; k < RUNS; k++) {
        double t0 = now_ns(); uint64_t c0 = __rdtsc();
        r.sum = fn(reps);
        uint64_t c1 = __rdtsc(); double t1 = now_ns();
        ns[k] = (t1 - t0) / ((double)reps * N); tk[k] = (double)(c1 - c0) / ((double)reps * N);
    }
    qsort(ns, RUNS, sizeof ns[0], cmp_d); qsort(tk, RUNS, sizeof tk[0], cmp_d);
    r.ns = ns[RUNS / 2]; r.tsc = tk[RUNS / 2]; r.ns_min = ns[0]; r.tsc_min = tk[0];
    return r;
}
static double tick_per_cycle = 1.0;              /* set by calibration in main */
static void line(const char *name, result r, result base)
{ printf("    %-42s %6.3f ns  ~%5.1f cycles  %5.2fx native   (median %6.3f ns)\n",
         name, r.ns_min, r.tsc_min / tick_per_cycle, r.tsc_min / base.tsc_min, r.ns); }

static void fail(const char *what) { fprintf(stderr, "SELF-TEST FAILED: %s\n", what); exit(1); }
static int same(reg x, reg y) { return x.a == y.a && x.b == y.b; }

static void self_test(void)
{
    long n = 0;
    for (int64_t m = -300; m <= 300; m++) {
        reg X = canon(m);
        if (!(X.b == 0 || X.b == 1) || value(X) != m) fail("canon");
        for (int64_t k = -300; k <= 300; k++) {
            reg Y = canon(k), want = canon(m * k), sum = canon(m + k);
            if (value(mul_raw(X, Y)) != m * k) fail("mul_raw value");
            if (!same(mul_norm_ref(X, Y), want)) fail("mul_norm_ref");
            if (!same(mul_canon(X, Y), want)) fail("mul_canon");
            if (!same(mul_via_value(X, Y), want)) fail("mul_via_value");
            if (!same(add_norm(X, Y), sum)) fail("add_norm");
            n += 5;
        }
    }
    for (int i = 0; i < 2000000; i++) {
        int64_t m = rnd_signed(30), k = rnd_signed(30);
        reg X = canon(m), Y = canon(k), want = canon(m * k);
        if (!same(mul_norm_ref(X, Y), want) || !same(mul_canon(X, Y), want) ||
            !same(mul_via_value(X, Y), want) || value(mul_raw(X, Y)) != m * k ||
            !same(add_norm(X, Y), canon(m + k))) fail("random 30-bit operands");
        n += 5;
        reg R = { rnd_signed(40), rnd_signed(40) };       /* arbitrary representative */
        if (!same(normalize(R), canon(value(R)))) fail("normalize");
        n += 1;
    }
    printf("self-test: %ld checks passed (all variants agree with the native result)\n", n);
}

int main(int argc, char **argv)
{
    if (argc > 1 && strcmp(argv[1], "--dump") == 0) {      /* vectors for rau.py */
        for (int64_t m = -25; m <= 25; m++) for (int64_t k = -25; k <= 25; k++) {
            reg X = canon(m), Y = canon(k), r = mul_raw(X, Y), c = mul_canon(X, Y),
                v = mul_via_value(X, Y), s = add_norm(X, Y);
            printf("%lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld\n",
                   (long long)m, (long long)k, (long long)X.a, (long long)X.b, (long long)Y.a, (long long)Y.b,
                   (long long)r.a, (long long)r.b, (long long)c.a, (long long)c.b,
                   (long long)v.a, (long long)v.b, (long long)s.a, (long long)s.b);
        }
        return 0;
    }
    self_test();
    for (int i = 0; i < N; i++) {
        XI[i] = rnd_signed(30); YI[i] = rnd_signed(30);
        XR[i] = canon(XI[i]);   YR[i] = canon(YI[i]);
        UI[i] = (rnd() & 1) ? 1 : -1;          UR[i] = canon(UI[i]);
        VI[i] = rnd_signed(20);                VR[i] = canon(VI[i]);
    }
    long reps = 200000;                          /* 2.0e8 operations per run */
    uint64_t sink = 0;

    result ln = measure(lat_native_mul, reps);        /* calibration: 3 cycles */
    tick_per_cycle = ln.tsc_min / 3.0;
    printf("calibration: native multiply latency %.3f ns = 3 cycles, so 1 cycle = %.3f ns "
           "(core at about %.2f GHz)\n", ln.ns_min, ln.ns_min / 3.0, 3.0 / ln.ns_min);

    printf("\nMULTIPLICATION, operand values below 2^30, canonical pairs\n");
    printf("  throughput (independent operations): fastest of %d runs\n", RUNS);
    result nm = measure(tp_native_mul, reps);          line("native  x * y", nm, nm);                       sink ^= nm.sum;
    result a  = measure(tp_mul_raw, reps);             line("raw product, formula (5)", a, nm);              sink ^= a.sum;
    result b  = measure(tp_mul_norm_ref, reps);        line("formula (5) + normalisation (3d)", b, nm);      sink ^= b.sum;
    result c  = measure(tp_mul_canon, reps);           line("shifts and masks, one product (Table 3)", c, nm); sink ^= c.sum;
    result d  = measure(tp_mul_via_value, reps);       line("decode, one product, encode", d, nm);           sink ^= d.sum;
    if (b.sum != c.sum || c.sum != d.sum) fail("normalised variants disagree on the benchmark data");
    printf("  latency (each result feeds the next product)\n");
    line("native  acc * x", ln, ln);                                                                         sink ^= ln.sum;
    result lb = measure(lat_mul_norm_ref, reps);       line("formula (5) + normalisation (3d)", lb, ln);     sink ^= lb.sum;
    result lc = measure(lat_mul_canon, reps);          line("shifts and masks, one product (Table 3)", lc, ln); sink ^= lc.sum;
    result ld = measure(lat_mul_via_value, reps);      line("decode, one product, encode", ld, ln);          sink ^= ld.sum;
    if (ln.sum != lb.sum || lb.sum != lc.sum || lc.sum != ld.sum) fail("chains disagree");

    printf("\nADDITION, canonical pairs, normalised result\n");
    printf("  throughput\n");
    result na = measure(tp_native_add, reps);          line("native  x + y", na, na);                        sink ^= na.sum;
    result ra = measure(tp_add_norm, reps);            line("register sum (4) + normalisation (3d)", ra, na); sink ^= ra.sum;
    printf("  latency\n");
    result la = measure(lat_native_add, reps);         line("native  acc + x", la, la);                      sink ^= la.sum;
    result lr = measure(lat_add_norm, reps);           line("register sum (4) + normalisation (3d)", lr, la); sink ^= lr.sum;
    if (la.sum != lr.sum) fail("addition chains disagree");
    printf("\nchecksum %016llx\n", (unsigned long long)sink);
    return 0;
}
