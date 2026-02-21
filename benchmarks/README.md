# Benchmarks Layout

This directory is split by purpose:

- `perf/`: performance benchmarks and timing harnesses.
- `validation/`: correctness/parity probes and option sweeps.
- `common/`: shared helpers used by benchmark scripts.

## Quick Start

Single one-off run:

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/oneoff/bench_oneoff_param.R \
  --fun=npunitest --n=500 --times=5 --seed_policy=fixed --base_seed=42 \
  --out_raw=/tmp/np_oneoff_raw.csv --out_summary=/tmp/np_oneoff_summary.csv
```

Method-combo run (`npreg`):

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npreg/run_npreg_combos.R
```

## Reproducibility Policy

- Use fixed-seed and varying-seed runs for performance comparisons.
- Pre/post comparisons should use matched seeds and identical arguments.
- Persist artifacts in `/tmp` and keep input/raw/summary outputs together.
