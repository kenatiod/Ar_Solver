# Ar_Solver

Certified computation of the **discrete RAD-CRT floor** `A_r` (OEIS A118478(r)) for the
prime-complete consecutive-product problem.

`A_r` is the smallest integer `m` at or above the primorial floor such that
`m(m+1)` is divisible by the product of the first `r` primes. It is the central
quantity in the **RAD-CRT wall**, an exact exclusion mechanism for prime-complete
products of two consecutive integers. This program computes `A_r` two ways —
exhaustive enumeration and meet-in-the-middle — and cross-checks them.

## Background

A product `m(m+1)` is called **prime-complete** when its radical (its set of
distinct prime factors) is a primorial: exactly the first `r` primes for some
`r`. The associated conjecture is that

> `633555 × 633556` is the **last** prime-complete product of two consecutive
> integers.

A prime-complete candidate at order `r` is squeezed between two bounds:

- a **floor**: `m(m+1)` must be divisible by `P_r` (the product of the first `r`
  primes), which forces `m(m+1) ≥ P_r`;
- a **ceiling**: `m` and `m+1` must both be `p_r`-smooth, so `m ≤ L_r`, where
  `L_r = A002072(r)` is the largest such smooth number.

For each prime `p | P_r`, divisibility of `m(m+1)` by `p` forces `m ≡ 0` or
`m ≡ p-1 (mod p)` — two choices per prime. By the Chinese Remainder Theorem
there are exactly `2^r` admissible residues ("CRT roots") of `x(x+1) ≡ 0
(mod P_r)`. A prime-complete `m` cannot sit just anywhere above the floor; it
must land on one of these discrete roots.

**`A_r` is the smallest CRT root at or above the floor.** The exact exclusion
condition at order `r` is therefore

> `A_r > L_r` ⟹ no prime-complete product exists at order `r`.

Computed exactly, this condition turns on at `r = 14` and the gap
`log A_r − log L_r` grows monotonically thereafter — the **RAD-CRT wall**,
many levels earlier than the older `L_r` vs `√P_r` crossing at `r = 30`.

## What this program computes

For each `r`, the program reports `A_r`, `log A_r`, and the deficiency
decomposition

```
log A_r = θ(p_r) − r·ln2 − d_r
```

where `θ(p_r) = log P_r`. The term `θ − r·ln2` is the parameter-free spacing
prediction (the typical gap `P_r / 2^r` between neighboring CRT roots, on the
log scale), and `d_r` is the measured deficiency. The ratio `d_r / θ` is the
diagnostic for the open floor lemma `d_r = o(θ)`.

Output columns: `r`, `p_r`, `log_Pr_theta`, `half_log_Pr`, `r_ln2`, `log_Ar`,
`Delta_r`, `deficiency_dr`, `dr_over_theta`, `A_r`.

Each `A_r` is the **smallest integer that could be a prime-complete consecutive 
product at order `r`** — the lowest admissible CRT root above the primorial floor. 
The proof works by showing this lower bound has risen _above_ the smooth ceiling 
`L_r` (the largest such product that smoothness permits). Once `A_r > L_r`, the 
window has closed: the cheapest possible candidate already exceeds the most 
expensive allowed one, so no prime-complete product can exist at that order. 
The values in these tables are the certified left-hand side of that inequality. 
Their growth — far outpacing `L_r` — is the wall that seals off every order 
beyond the last known solution.

## Methods

The program offers three computation methods plus an automatic dispatcher.

| Method | Memory | Reach | Certification status |
|---|---|---|---|
| `stream` | `O(r)` | bounded by time (`2^r` leaves) | **exhaustive** — every root visited; minimality automatic |
| `direct` | `O(2^r)` | ~`r=22` before memory pressure | exhaustive, but materializes the full root list |
| `mitm` | `O(2^{r/2})` | far higher (seconds per step) | algorithmic; values cross-checked against the exhaustive oracle |
| `auto` | — | — | `stream` through `r=24`, then `mitm` |

**`stream`** is the recommended exhaustive method. It enumerates all `2^r` CRT
roots by depth-`r` recursion, holding only the current root path in memory
(never the full list). It returns a completeness witness (`n_roots == 2^r`) and
the winning root's per-prime side-assignment mask. Despite being exhaustive it
is both lower-memory and faster than the list-based `direct` method.

**`mitm`** (meet-in-the-middle) splits the primes into two halves and matches
partial residues against the floor. It reaches much larger `r` cheaply. It is
deterministic and exact; its outputs agree with the exhaustive oracle wherever
both are run (verified through `r = 24`). What it does not do by itself is prove
*minimality* by exhaustion — that guarantee comes from the algorithm plus the
oracle agreement, a slightly weaker tier than `stream`.

### Two certification tiers

- **Exhaustively certified** (`stream`/`direct` range): `A_r` is the least root
  above the floor because every one of the `2^r` roots was examined.
- **Algorithmically certified, oracle-validated** (`mitm` tail): `A_r` is exact
  and deterministic, validated against the exhaustive oracle on the overlap
  range. Suitable as supporting evidence for the wall's behavior at large `r`.

Both tiers are exact integer computations — neither is heuristic.

## Floor conventions

Two floor definitions are supported, differing by less than 1:

- `--floor sqrt` (default): `F_r = ⌈√P_r⌉`, the bare primorial floor.
- `--floor exact`: `F*_r = ⌈(√(1+4·P_r) − 1)/2⌉`, the exact floor from the
  product condition `m(m+1) ≥ P_r`.

The exact floor is the correct one for proof purposes (it is what divisibility
actually requires) and is the *lower* of the two, hence the safer one — it admits
weakly more candidates. The two conventions yield the **same `A_r` for all
`r ≥ 8`** in the verified range; they differ only at `r ∈ {2, 3, 4, 7}`, where
the exact floor correctly admits boundary prime-complete products (for example
`714 × 715 = P_7` at `r = 7`) and where the finite base is established by direct
enumeration in any case. The wall region (`r ≥ 14`) is identical under either
convention.

The `--verify-gap N` option checks, for all `r ≤ N`, whether any CRT root lies
in the gap `[F*_r, ⌈√P_r⌉)`, reporting exactly where the two floors differ.

## Requirements

- Python 3.8+
- [`gmpy2`](https://pypi.org/project/gmpy2/) (recommended): `pip install gmpy2`
  or `apt install python3-gmpy2`. The program runs without it via a pure-Python
  fallback, but big-integer arithmetic is several times slower.

## Usage

```bash
# Exhaustive certified run with the exact (proof-clean) floor, write CSV
python3 Ar_solver.py --rmin 2 --rmax 24 --method stream --floor exact --csv Ar_certified.csv

# Fast extension into the large-r tail (meet-in-the-middle)
python3 Ar_solver.py --rmin 25 --rmax 54 --method mitm --floor exact

# Cross-check mitm against the exhaustive oracle on the overlap range
python3 Ar_solver.py --rmin 2 --rmax 24 --method mitm --crosscheck 24 --floor exact

# Verify the two floor conventions agree (and see exactly where they differ)
python3 Ar_solver.py --rmin 2 --rmax 16 --method stream --floor exact --verify-gap 16

# Sanity print of the wall-crossing neighborhood (r=12..16)
python3 Ar_solver.py --check-r14
```

### Command-line options

| Option | Description |
|---|---|
| `--rmin`, `--rmax` | order range (inclusive) |
| `--method` | `stream`, `direct`, `mitm`, or `auto` (default) |
| `--floor` | `sqrt` (default) or `exact` |
| `--crosscheck N` | for all `r ≤ N`, assert the chosen method agrees with the exhaustive oracle |
| `--verify-gap N` | for all `r ≤ N`, check no root lies in `[F*_r, ⌈√P_r⌉)` |
| `--csv PATH` | write the full decomposition table to CSV |
| `--check-r14` | print the wall-crossing neighborhood and confirm method agreement |

## Performance notes

The exhaustive methods scale as `2^r` in time (each added prime doubles the leaf
count); expect roughly a doubling of runtime per step. The `stream` method's
memory stays flat regardless of `r` — a run at `r = 24` (enumerating
`2^24 ≈ 1.68 × 10^7` roots) completes in well under 100 MB of resident memory,
including the interpreter and `gmpy2`. The `mitm` method scales as `2^{r/2}` and
reaches far larger `r` in seconds per step, bounded in practice by the memory of
its sorted half-residue table.

## Files

- `Ar_solver.py` — the solver.
- `*.csv` — certified output tables (decomposition columns plus exact `A_r`).

## Related work

- OEIS [A002072](https://oeis.org/A002072) — the smooth ceiling `L_r`.
- Companion repositories: `LC_Solver`, `Delta_min`, `CRT_Pruning_Survey`,
  `A002072_Solver`.

## License

MIT License — see [LICENSE](LICENSE).
