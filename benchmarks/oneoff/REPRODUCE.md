# Reproduce: one-off benchmarks (npRmpi)

Example synthetic run (`npunitest`):

```bash
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/oneoff/bench_oneoff_param_nprmpi.R \
  --fun=npunitest --n=1000 --nslaves=1 --times=5 --seed_policy=varying --base_seed=42 \
  --out_raw=/tmp/nprmpi_oneoff_npunitest_raw.csv \
  --out_summary=/tmp/nprmpi_oneoff_npunitest_summary.csv
```

Example real-data run (`npqreg`, ignores `--n`):

```bash
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/oneoff/bench_oneoff_param_nprmpi.R \
  --fun=npqreg --n=1000 --nslaves=1 --times=5 --seed_policy=fixed --base_seed=42 \
  --out_raw=/tmp/nprmpi_oneoff_npqreg_raw.csv \
  --out_summary=/tmp/nprmpi_oneoff_npqreg_summary.csv
```
