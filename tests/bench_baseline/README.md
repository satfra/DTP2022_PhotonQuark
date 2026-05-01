# Frozen Bench Baseline

Per-tier numerical-regression baseline for the performance work.
Captured once on `main`-equivalent code at the bench grid sizes
(see `Quark_Photon_Vertex/include/parameters_bench.hh`).

## Capture

```bash
bash bench/run_bench.sh tests/bench_baseline
```

The four sub-directories `nodse_nopv/`, `nodse_pv/`, `dse_nopv/`, `dse_pv/`
correspond to the executable's flag combinations `""`, `p`, `d`, `dp`.

## Diff

```bash
bash bench/run_bench.sh                     # writes to bench/results/current/
bash bench/diff_baseline.sh                 # tolerance-aware diff vs this dir
```

Defaults: `RTOL=1e-6 ATOL=1e-10`, well above the iteration target accuracy
of `1e-5` defined in `parameters.hh`.

## Why this is git-ignored

Outputs are platform-dependent under `-ffast-math` (CPU rounding) and large
(~56 MB total). Each developer / CI worker captures their own baseline
before starting performance work, then keeps it frozen across tiers.
