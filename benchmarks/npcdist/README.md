# npcdist Benchmark Harness (np)

This folder provides a parameterized benchmark harness for `npcdistbw()` + `npcdist()` in serial `np` using `microbenchmark`.

## Files

- `bench_npcdist_param.R`: run one benchmark configuration with `microbenchmark` repetitions.
- `run_npcdist_combos.R`: run all 48 option combinations (`2x3x2x2x2`) with a shared setup.
- `compare_npcdist_versions.R`: compare two raw outputs and emit overall and combo-specific comparison tables.
- `make_npcdist_report.R`: generate a markdown report from comparison CSV outputs.
- `make_npcdist_combined_report.R`: generate a single markdown report combining `np` and `npRmpi` comparison CSV outputs.
- `REPRODUCE.md`: end-to-end current-vs-CRAN replication steps.

## One-Configuration Run

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/npcdist/bench_npcdist_param.R \
  --n=100 --times=50 --base_seed=42 \
  --bwmethod=cv.ls --bwtype=fixed --nmulti=1 \
  --cxkertype=gaussian --cykertype=gaussian \
  --np_tree=FALSE --seed_policy=fixed \
  --out_raw=/tmp/npcdist_one_raw.csv --out_summary=/tmp/npcdist_one_summary.csv
```

## Full Combination Run (48 combos)

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/npcdist/run_npcdist_combos.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --tag=myrun
```

## Important Defaults

- `n=100`
- `times=50`
- `base_seed=42`
- `nmulti=1`
- `bwmethod=cv.ls`
- `cxkertype=gaussian`
- `cykertype=gaussian`
- `np_tree=FALSE`
- `seed_policy=varying`

Notes:

- `bwmethod` choices used in combos are `cv.ls` and `normal-reference`.
- paired-kernel rule is enforced: `cxkertype == cykertype`, with pair choices `gaussian` or `epanechnikov`.

## DGP

The benchmark uses:

```r
set.seed(seed)
x <- runif(n)
z <- factor(rbinom(n, 3, .5))
y <- x + as.numeric(as.character(z)) + rnorm(n, sd = .25)
```

and fits `npcdistbw(y ~ x + z, ...)` then `npcdist(bws=bw, data=...)`.
