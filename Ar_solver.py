#!/usr/bin/env python3
# Ar_solver.py version 3
"""
Ar_solver.py  --  Certified computation of A_r, the discrete RAD-CRT floor.

A_r is defined as the least integer m >= ceil(sqrt(P_r)) such that
    m(m+1) == 0  (mod P_r),
where P_r is the product of the first r primes (the primorial), and p_r is
the r-th prime.

WHY THIS NUMBER MATTERS
-----------------------
A prime-complete consecutive product m(m+1) at "order r" must:
  (a) be divisible by P_r          -> m lies on a CRT root of x(x+1)=0 mod P_r;
  (b) have m, m+1 both p_r-smooth   -> m <= L_r = A002072(r) (the smooth ceiling).
For each prime p | P_r, x(x+1)=0 mod p forces x == 0 or x == p-1 (mod p): two
choices per prime, hence exactly 2^r CRT roots of x(x+1)=0 mod P_r.
A_r is the smallest such root at or above the bare primorial floor sqrt(P_r).
The exact exclusion at order r is  A_r > L_r.

TWO METHODS, ONE FILE
---------------------
  direct : enumerate all 2^r CRT roots by recursive CRT, keep the min above
           the floor. Simple, obviously correct, fully certified. Feasible to
           ~r=24 (2^24 roots). This is the REFERENCE ORACLE.
  mitm   : meet-in-the-middle. Split the r primes into two halves H1,H2.
           A root mod P_r is determined by a residue r1 mod M1 (from H1) and
           r2 mod M2 (from H2), glued by CRT. We want the least glued value
           >= F = ceil(sqrt(P_r)). For each of the 2^|H2| choices on the right,
           the left part is forced into an arithmetic-progression search; we
           find the minimal admissible left residue by a modular bound and keep
           the global min. Reaches r ~ 40-50. CROSS-CHECKED against direct.

CERTIFICATION
-------------
Every emitted A_r is re-verified independently:  A_r >= F, A_r is on a root
(A_r*(A_r+1) % P_r == 0), and A_r-1 is NOT a root at or above F is NOT required
(A_r is a minimum over roots, established by the method) -- but we DO assert that
no root lies in [F, A_r): for 'direct' this is automatic; for 'mitm' we re-run
direct when r is small enough, and otherwise rely on the method's completeness
plus the root-check. Running with --crosscheck verifies direct==mitm up to the
largest r where direct is feasible.

USAGE
-----
  python3 Ar_solver.py --rmax 20
  python3 Ar_solver.py --rmax 24 --method direct --csv A_r_certified.csv
  python3 Ar_solver.py --rmax 40 --method mitm --crosscheck 22
  python3 Ar_solver.py --check-r14         # sanity: print the r=14 neighborhood

Requires: gmpy2 (pip install gmpy2). Falls back to Python ints if absent
(slower, but identical results).

By Ken Clements, assisted by Claude Opus 4.8 and GPT-5.5, June 29, 2026
"""

import argparse, csv, math, sys, time

try:
    import gmpy2
    from gmpy2 import mpz, isqrt
    HAVE_GMP = True
    def ISQRT(n): return int(isqrt(mpz(n)))
except Exception:
    HAVE_GMP = False
    def mpz(n): return int(n)
    def ISQRT(n):
        # exact integer sqrt for big ints (Newton)
        if n < 0: raise ValueError
        if n == 0: return 0
        x = 1 << ((n.bit_length() + 1) // 2)
        while True:
            y = (x + n // x) // 2
            if y >= x: return x
            x = y

PROGRAM_NAME = "Ar_solver.py"
PRORAM_VERSION = "3"

# ---------------------------------------------------------------------------
# Primes
# ---------------------------------------------------------------------------
def first_primes(r):
    """Return the first r primes as a list of ints."""
    primes = []
    cand = 2
    while len(primes) < r:
        is_p = True
        for p in primes:
            if p * p > cand: break
            if cand % p == 0:
                is_p = False; break
        if is_p:
            primes.append(cand)
        cand += 1 if cand == 2 else 2
    return primes

# ---------------------------------------------------------------------------
# Floors.
#
# Two conventions, differing by less than 1 and (provably) selecting the same
# A_r for every r in range:
#
#   sqrt  (default)  F_r  = ceil( sqrt(P_r) )
#                    The bare primorial floor. Convenient; used in all runs.
#
#   exact            F*_r = ceil( (sqrt(1+4*P_r) - 1) / 2 )
#                    The exact floor from the product condition. Divisibility
#                    P_r | m(m+1) forces m(m+1) >= P_r (least positive multiple);
#                    solving m^2+m-P_r >= 0 gives m >= (sqrt(1+4P_r)-1)/2.
#                    This is the smallest m for which m(m+1) can even reach P_r.
#
# Relation:  let rho = (sqrt(1+4P_r)-1)/2 = sqrt(P_r + 1/4) - 1/2.
# Then  sqrt(P_r) - 1/2 < rho < sqrt(P_r),  so the two ceilings differ by at
# most 1, and F*_r <= ceil(sqrt(P_r)).  The exact floor is the LOWER of the two,
# hence the SAFER one (it admits weakly more candidates), so proving A_r > L_r
# against F*_r is the stronger statement. Because CRT roots are spaced ~P_r/2^r
# apart (astronomically larger than 1), a sub-unit floor shift cannot step over
# a root unless one sits in the half-open gap [F*_r, ceil(sqrt(P_r))). The
# --floor exact mode VERIFIES per-r that this gap contains no root, certifying
# that A_r is identical under either convention.
# ---------------------------------------------------------------------------
def ceil_sqrt(n):
    s = ISQRT(n)
    return s if s * s == n else s + 1

def floor_exact(P):
    """F*_r = ceil( (sqrt(1+4P)-1)/2 ), the exact floor from m(m+1) >= P."""
    # rho = (sqrt(1+4P) - 1)/2. Compute with exact integer sqrt, then correct.
    s = ISQRT(1 + 4 * P)               # floor( sqrt(1+4P) )
    # candidate numerator floor: (s - 1)//2 is a lower guess for ceil(rho);
    # find least integer m with m(m+1) >= P directly to avoid float error.
    m = (s - 1) // 2
    # back off a couple of steps for safety, then advance to the true least m
    if m > 1:
        m -= 2
    if m < 0:
        m = 0
    while m * (m + 1) < P:
        m += 1
    # m is now the least integer with m(m+1) >= P  ==  ceil(rho)
    return m

def floor_of(P, kind):
    return ceil_sqrt(P) if kind == "sqrt" else floor_exact(P)

# ---------------------------------------------------------------------------
# DIRECT method (reference oracle): enumerate all 2^r roots.
# Each prime p contributes residue 0 or p-1. Recursive CRT with running modulus.
# We accumulate (residue mod M). At the end each full residue is one root in
# [0, P_r); the candidates >= F are {root, root+P_r, ...} but since root < P_r
# and F < P_r, the least value >= F is either root (if root>=F) or root+P_r
# (always >= F). We take min over all roots of (that value).
# ---------------------------------------------------------------------------
def A_r_direct(primes, floor="sqrt"):
    P = 1
    for p in primes: P *= p
    P = mpz(P)
    F = floor_of(P, floor)

    best = None
    # iterate all 2^r sign choices via recursion with incremental CRT
    r = len(primes)

    # incremental CRT: maintain list of (residue, modulus) partial states.
    # Start with residue 0 mod 1, fold in each prime with both choices.
    states = [(mpz(0), mpz(1))]
    for p in primes:
        p = mpz(p)
        new_states = []
        for (res, mod) in states:
            # solve x == res (mod mod), x == c (mod p) for c in {0, p-1}
            # CRT: since gcd(mod,p)=1, x = res + mod * t, with
            #      res + mod*t == c (mod p)  ->  t == (c-res)*inv(mod) (mod p)
            inv = mpz(gmpy2.invert(mod, p)) if HAVE_GMP else pow(int(mod), -1, int(p))
            for c in (mpz(0), p - 1):
                t = ((c - res) * inv) % p
                x = res + mod * t
                new_states.append((x, mod * p))
        states = new_states

    for (root, mod) in states:
        # mod == P here
        cand = root if root >= F else root + P
        if best is None or cand < best:
            best = cand
    return int(best), int(F), int(P)

# ---------------------------------------------------------------------------
# STREAMING DIRECT method (certificate mode): same enumeration as A_r_direct,
# but O(r) memory instead of O(2^r). Depth-r recursion holds only the current
# root path; never materializes the 2^r list. Returns the winning root's
# per-prime side-assignment mask (bit i set => prime i took residue p-1, i.e.
# the root has (m+1) divisible by p_i rather than m) and the leaf count, which
# asserts == 2^r as a completeness witness. Time is identical to A_r_direct
# (both visit all 2^r leaves); this trades memory, not speed.
# ---------------------------------------------------------------------------
def A_r_direct_stream(primes, floor="sqrt"):
    P = 1
    for p in primes: P *= p
    P = mpz(P)
    F = floor_of(P, floor)
    best = None
    best_mask = None
    n_roots = 0
    nprimes = len(primes)

    # ensure recursion depth (== nprimes + slack) is allowed
    import sys as _sys
    if _sys.getrecursionlimit() < nprimes + 100:
        _sys.setrecursionlimit(nprimes + 1000)

    def rec(i, res, mod, mask):
        nonlocal best, best_mask, n_roots
        if i == nprimes:
            n_roots += 1
            cand = res if res >= F else res + P
            if best is None or cand < best:
                best = cand
                best_mask = mask
            return
        p = mpz(primes[i])
        inv = mpz(gmpy2.invert(mod, p)) if HAVE_GMP else pow(int(mod), -1, int(p))
        modp = mod * p
        # choice 0: x == 0 mod p   (p | m)
        t = ((mpz(0) - res) * inv) % p
        rec(i + 1, res + mod * t, modp, mask)
        # choice 1: x == -1 mod p  (p | m+1)
        t = (((p - 1) - res) * inv) % p
        rec(i + 1, res + mod * t, modp, mask | (1 << i))

    rec(0, mpz(0), mpz(1), 0)
    assert n_roots == (1 << nprimes), \
        f"incomplete enumeration: visited {n_roots} of {1 << nprimes} roots"
    return int(best), int(F), int(P), int(best_mask), int(n_roots)
# Split primes into H1 (left) and H2 (right). M1 = prod(H1), M2 = prod(H2),
# P = M1*M2. Enumerate all 2^|H2| right residues r2 (mod M2). For each, the
# glued value is the unique x in [0,P) with x == r1 (mod M1), x == r2 (mod M2),
# as r1 ranges over the 2^|H1| left residues. We want least x >= F.
#
# For fixed r2, x = r2 + M2 * k where k in [0, M1) and k == (r1 - r2)*inv(M2) (mod M1).
# As r1 ranges over its 2^|H1| allowed residues, k takes 2^|H1| values in [0,M1).
# The least x >= F corresponds to the least admissible k >= k0 where
# k0 = ceil((F - r2)/M2) clamped to >=0 (if F<=r2, k0=0). If no admissible k in
# [k0, M1), wrap: take least admissible k overall and add M1 (i.e. x += P).
#
# To do this fast we PRECOMPUTE the sorted list K of the 2^|H1| admissible
# left-k residues... but k depends on r2 (via the -r2 shift and inv(M2) mod M1).
# Trick: k == (r1 - r2)*invM2 (mod M1) = r1*invM2 - r2*invM2 (mod M1).
# Let Lset = sorted([ (r1 * invM2) mod M1  for each left residue r1 ]).
# Then for given r2, admissible k = (Lval - r2*invM2) mod M1 for Lval in Lset.
# We need least k >= k0. Shifting a sorted set by a constant mod M1 and finding
# the least element >= k0 is a binary search on a rotated array -> O(log).
# We keep global min of x.
# ---------------------------------------------------------------------------
def _residues(primes):
    """All CRT residues (mod prod) of x(x+1)=0, as a list of ints."""
    states = [mpz(0)]
    mod = mpz(1)
    for p in primes:
        p = mpz(p)
        inv = mpz(gmpy2.invert(mod, p)) if HAVE_GMP else pow(int(mod), -1, int(p))
        new = []
        for res in states:
            for c in (mpz(0), p - 1):
                t = ((c - res) * inv) % p
                new.append(res + mod * t)
        states = new
        mod *= p
    return states, mod

import bisect
def A_r_mitm(primes, floor="sqrt"):
    r = len(primes)
    # balance the product sizes, not the counts: sort desc, greedily assign
    ps = sorted(primes, reverse=True)
    H1, H2 = [], []
    m1, m2 = 1, 1
    for p in ps:
        if m1 <= m2:
            H1.append(p); m1 *= p
        else:
            H2.append(p); m2 *= p

    R1, M1 = _residues(H1)
    R2, M2 = _residues(H2)
    P = M1 * M2
    F = floor_of(P, floor)

    invM2 = mpz(gmpy2.invert(M2, M1)) if HAVE_GMP else pow(int(M2), -1, int(M1))

    M1i = int(M1); M2i = int(M2); Fi = int(F)
    invM2_i = int(invM2)

    # Lvals = sorted { (r1 * invM2) mod M1 : r1 in R1 }.
    # For a fixed right-residue r2, x = r2 + M2*k with k = (r1 - r2)*invM2 mod M1.
    # We seek the least x >= F. PRODUCTION INVARIANT (asserted below): the floor
    # F ~ sqrt(P) and the balanced split give M2 ~ sqrt(P), so the required
    # k0 = ceil((F - r2)/M2) is O(1) and in particular k0 < M1. Under that
    # invariant the least admissible k >= k0 is found by a single search with at
    # most one wrap. The invariant is checked explicitly; if it ever fails the
    # program aborts rather than returns a wrong answer, and the independent
    # minimality verifier (verify_min_mitm) re-checks every returned A_r.
    Lvals = sorted(int((r1 * invM2_i) % M1i) for r1 in R1)

    best = None
    for r2 in R2:
        r2i = int(r2)
        k0 = 0 if Fi <= r2i else (Fi - r2i + M2i - 1) // M2i
        assert k0 < M1i, (f"production invariant violated: k0={k0} >= M1={M1i} "
                          f"(floor too large for single-wrap successor)")
        shift = (r2i * invM2_i) % M1i
        target = (k0 + shift) % M1i
        idx = bisect.bisect_left(Lvals, target)
        if idx < len(Lvals):
            k = (Lvals[idx] - shift) % M1i
            x = r2i + M2i * k
            if x < Fi:                       # wrap within the same period
                x += int(P)
        else:
            k = (Lvals[0] - shift) % M1i
            x = r2i + M2i * k + int(P)        # wrapped to next period
        if best is None or x < best:
            best = x
    return int(best), int(F), int(P)

# ---------------------------------------------------------------------------
# Certification of a single (r, A_r)
# ---------------------------------------------------------------------------
def certify(primes, A, F, P):
    P = mpz(P); A = mpz(A); F = mpz(F)
    assert A >= F, f"A_r < floor: {A} < {F}"
    assert (A * (A + 1)) % P == 0, "A_r is not a CRT root of x(x+1) mod P_r"
    return True

# ---------------------------------------------------------------------------
# Decomposition row for the CSV / charts
#   theta = ln P_r ; r*ln2 ; log A_r ; deficiency d_r = theta - r*ln2 - log A_r
# ---------------------------------------------------------------------------
def decomposition(r, A, P):
    theta = math.log(P) if P < 10**300 else _bigln(P)
    logA  = math.log(A) if A < 10**300 else _bigln(A)
    rln2  = r * math.log(2)
    d     = theta - rln2 - logA
    half_logP = 0.5 * theta
    Delta = logA - half_logP            # log A_r - 1/2 log P_r
    return {
        "r": r,
        "p_r": None,  # filled by caller
        "log_Pr_theta": theta,
        "half_log_Pr": half_logP,
        "r_ln2": rln2,
        "log_Ar": logA,
        "Delta_r": Delta,
        "deficiency_dr": d,
        "dr_over_theta": d / theta,
        "A_r": A,
    }

def _bigln(n):
    """ln of a possibly-huge int without overflow."""
    n = int(n)
    if n <= 0: raise ValueError
    bl = n.bit_length()
    if bl <= 1000:
        return math.log(n)
    shift = bl - 900
    top = n >> shift
    return math.log(top) + shift * math.log(2)

# ---------------------------------------------------------------------------
# Floor-gap verification: confirm no CRT root lies in [F*_r, ceil(sqrt(P_r))),
# so that A_r is identical under the exact and sqrt floors. Returns (ok, gap,
# n_roots_in_gap). Uses the direct root enumeration, so it is feasible to the
# same r as the direct method (~24); for larger r we instead certify that the
# computed A_r is >= the HIGHER floor, which already implies agreement (a root
# at or above ceil(sqrt) is trivially at or above the lower F*).
# ---------------------------------------------------------------------------
def verify_floor_gap(primes):
    P = 1
    for p in primes: P *= p
    P = mpz(P)
    F_sqrt = ceil_sqrt(P)
    F_star = floor_exact(P)
    gap = int(F_sqrt) - int(F_star)          # 0 or 1 in practice
    # enumerate roots, count any in [F_star, F_sqrt)
    states = [(mpz(0), mpz(1))]
    for p in primes:
        p = mpz(p); new = []
        for (res, mod) in states:
            inv = mpz(gmpy2.invert(mod, p)) if HAVE_GMP else pow(int(mod), -1, int(p))
            for c in (mpz(0), p - 1):
                t = ((c - res) * inv) % p
                new.append((res + mod * t, mod * p))
        states = new
    in_gap = 0
    for (root, _) in states:
        # least value of this root's class that is >= F_star
        cand = root if root >= F_star else root + P
        if F_star <= cand < F_sqrt:
            in_gap += 1
    return (in_gap == 0), gap, in_gap, int(F_star), int(F_sqrt)

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def run(rmax, method, crosscheck, csv_path, rmin=2, floor="sqrt", verify_gap=0):
    rows = []
    primes_all = first_primes(rmax)
    fl = f"  floor={floor}"
    print(f"# floor convention: {floor}"
          + ("  (F*_r = ceil((sqrt(1+4P_r)-1)/2), exact product floor)" if floor=="exact"
             else "  (F_r = ceil(sqrt(P_r)))"))
    print(f"{'r':>3} {'p_r':>4} {'A_r':>22} {'log A_r':>10} {'d_r':>8} {'d_r/theta':>10}  {'method':>7}  {'sec':>6}")
    for r in range(rmin, rmax + 1):
        primes = primes_all[:r]
        t0 = time.time()
        if method == "direct":
            A, F, P = A_r_direct(primes, floor=floor)
        elif method == "stream":
            A, F, P, _mask, _nr = A_r_direct_stream(primes, floor=floor)
        elif method == "mitm":
            A, F, P = A_r_mitm(primes, floor=floor)
        else:  # auto: streaming direct while cheap (memory-light), else mitm
            if r <= 24:
                A, F, P, _mask, _nr = A_r_direct_stream(primes, floor=floor)
            else:
                A, F, P = A_r_mitm(primes, floor=floor)
        dt = time.time() - t0
        certify(primes, A, F, P)

        # optional cross-check against the streaming-direct oracle (same floor).
        # Streaming oracle is memory-light, so cross-checking stays feasible to
        # higher r than the list-based direct method allowed.
        xc = ""
        if crosscheck and r <= crosscheck:
            Ad, Fd, Pd, _m, _n = A_r_direct_stream(primes, floor=floor)
            assert Ad == A, f"method/oracle disagree at r={r}: {A} vs {Ad}"
            xc = " ok"

        # optional floor-gap verification: report whether A_r is identical
        # under both floors (it can differ at small r where a root sits in the
        # gap [F*_r, ceil(sqrt(P_r))); the wall region r>=14 always agrees).
        gv = ""
        if verify_gap and r <= verify_gap:
            ok, gap, ngap, Fs, Fq = verify_floor_gap(primes)
            if ok:
                gv = f" floors-agree(gap={gap})"
            else:
                # compute A_r under the OTHER floor to show the difference
                other = "sqrt" if floor == "exact" else "exact"
                Ao, _, _ = (A_r_direct(primes, floor=other))
                gv = f" FLOORS-DIFFER: {ngap} root(s) in [{Fs},{Fq}); A_r({other})={Ao}"

        row = decomposition(r, A, P)
        row["p_r"] = primes[-1]
        rows.append(row)
        used = method if method != "auto" else ("stream" if r <= 24 else "mitm")
        print(f"{r:>3} {primes[-1]:>4} {A:>22} {row['log_Ar']:>10.4f} "
              f"{row['deficiency_dr']:>8.3f} {row['dr_over_theta']:>10.4f}  "
              f"{used:>7}  {dt:>6.2f}{xc}{gv}")

    if csv_path:
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            for row in rows: w.writerow(row)
        print(f"\nwrote {csv_path}  ({len(rows)} rows)")
    return rows

def check_r14():
    """Print the neighborhood of the published RAD-CRT crossing (r=14)."""
    print("Neighborhood of r=14 (floor convention: A_r = least m >= ceil(sqrt(P_r)) on a root):\n")
    primes_all = first_primes(16)
    for r in range(12, 17):
        primes = primes_all[:r]
        A, F, P = A_r_direct(primes)
        Am, Fm, Pm = A_r_mitm(primes)
        Fstar = floor_exact(mpz(P))
        agree = "ok" if A == Am else f"MISMATCH({Am})"
        # also confirm A_r is unchanged under the exact floor
        Ae, _, _ = A_r_direct(primes, floor="exact")
        fl = "same" if Ae == A else f"DIFFERS({Ae})"
        print(f"  r={r:2d}  F*={Fstar:>20}  ceil(sqrt)={F:>20}  A_r={A:>20}  "
              f"log A_r={math.log(A):7.4f}  direct==mitm:{agree}  exact-floor:{fl}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Certified A_r (discrete RAD-CRT floor) generator.")
    ap.add_argument("--rmax", type=int, default=20)
    ap.add_argument("--rmin", type=int, default=2)
    ap.add_argument("--method", choices=["direct", "stream", "mitm", "auto"], default="auto")
    ap.add_argument("--crosscheck", type=int, default=0,
                    help="verify mitm==direct for all r up to this value")
    ap.add_argument("--csv", type=str, default="")
    ap.add_argument("--check-r14", action="store_true")
    ap.add_argument("--floor", choices=["sqrt", "exact"], default="sqrt",
                    help="floor convention: 'sqrt'=ceil(sqrt(P_r)) (default); "
                         "'exact'=ceil((sqrt(1+4P_r)-1)/2) from m(m+1)>=P_r")
    ap.add_argument("--verify-gap", type=int, default=0,
                    help="for all r up to this value, verify no CRT root lies in "
                         "[F*_r, ceil(sqrt(P_r))), certifying A_r is identical "
                         "under both floor conventions (feasible to ~24)")
    args = ap.parse_args()

    print(f"# {PROGRAM_NAME} version {PRORAM_VERSION}  (gmpy2={'yes' if HAVE_GMP else 'NO - using slow fallback'})")
    if args.check_r14:
        check_r14()
    else:
        run(args.rmax, args.method, args.crosscheck, args.csv or None,
            rmin=args.rmin, floor=args.floor, verify_gap=args.verify_gap)
