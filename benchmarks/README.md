# Benchmarks Layout

This directory is split by purpose:

- `perf/`: performance benchmarks and timing harnesses.
- `validation/`: correctness/parity probes and option sweeps.
- `common/`: shared helpers used by benchmark scripts.

## Quick Start

Single one-off run:

```bash
FI_TCP_IFACE=en0 Rscript /Users/jracine/Development/np-npRmpi/benchmarks/perf/oneoff/bench_oneoff_param_nprmpi.R \
  --fun=npunitest --n=500 --nslaves=1 --times=5 --seed_policy=fixed --base_seed=42 \
  --out_raw=/tmp/nprmpi_oneoff_raw.csv --out_summary=/tmp/nprmpi_oneoff_summary.csv
```

Method-combo run (`npreg`):

```bash
FI_TCP_IFACE=en0 Rscript /Users/jracine/Development/np-npRmpi/benchmarks/perf/methods/npreg/run_npreg_combos.R
```

## MPI Environment Notes

- Preferred interface hint: `FI_TCP_IFACE=en0`.
- If unavailable, retry with `FI_TCP_IFACE=lo0`.
- If MPI provider init fails in a restricted environment, rerun outside sandbox.

## Reproducibility Policy

- Use fixed-seed and varying-seed runs for performance comparisons.
- Pre/post comparisons should use matched seeds and identical arguments.
- Persist artifacts in `/tmp` and keep input/raw/summary outputs together.

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
