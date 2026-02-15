# npreg Benchmark Harness (np)

This folder provides a parameterized benchmark harness for `npregbw()` + `npreg()` in serial `np`.

## Files

- `bench_npreg_param.R`: run one benchmark configuration (possibly multiple seeds).
- `run_npreg_combos.R`: run all 32 option combinations (`2x2x2x2x2`) with a shared setup.

## One-Configuration Run

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/npreg/bench_npreg_param.R \
  --n=100 --times=1 --base_seed=42 \
  --regtype=lc --bwmethod=cv.ls --nmulti=1 \
  --ckertype=gaussian --np_tree=FALSE --seed_policy=fixed \
  --out_raw=/tmp/npreg_one_raw.csv --out_summary=/tmp/npreg_one_summary.csv
```

## Full Combination Run (32 combos)

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/npreg/run_npreg_combos.R \
  --n=100 --times=1 --base_seed=42 --nmulti=1 --tag=myrun
```

Outputs are written to `/tmp` with run IDs in filenames:

- `npreg_combo_raw_...csv`
- `npreg_combo_summary_...csv`
- `npreg_combo_manifest_...csv`

## Important Defaults

- `n=100`
- `times=1`
- `base_seed=42`
- `nmulti=1`
- `regtype=ll`
- `bwmethod=cv.ls`
- `ckertype=gaussian`
- `np_tree=FALSE`
- `seed_policy=varying`

Change sample size and repetitions via `--n=` and `--times=`.

## Compatibility Note

`num_fval` is extracted from `bw$num.fval` when available. If missing (older versions), it is reported as `NA`.
