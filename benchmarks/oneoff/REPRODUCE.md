# Reproduce: one-off benchmarks (np)

Example synthetic run (`npunitest`):

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/oneoff/bench_oneoff_param.R \
  --fun=npunitest --n=1000 --times=5 --seed_policy=varying --base_seed=42 \
  --out_raw=/tmp/np_oneoff_npunitest_raw.csv \
  --out_summary=/tmp/np_oneoff_npunitest_summary.csv
```

Example real-data run (`npqreg`, ignores `--n`):

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/oneoff/bench_oneoff_param.R \
  --fun=npqreg --n=1000 --times=5 --seed_policy=fixed --base_seed=42 \
  --out_raw=/tmp/np_oneoff_npqreg_raw.csv \
  --out_summary=/tmp/np_oneoff_npqreg_summary.csv
```
