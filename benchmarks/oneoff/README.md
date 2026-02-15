# One-Off Benchmarks (np)

This directory provides a single benchmark entry-point for one-off `np` methods where we do not vary kernel/method options.

## Covered Functions

- `npcmstest`
- `npconmode`
- `npcopula`
- `npdeneqtest`
- `npdeptest`
- `npqreg`
- `npregiv`
- `npsdeptest`
- `npsigtest`
- `npsymtest`
- `npunitest`

## Inputs

- `--fun=` required function name from the list above.
- `--n=` sample size for synthetic DGP functions.
  - Ignored for real-data functions: `npcmstest`, `npconmode`, `npqreg`.
- `--times=` microbenchmark repetitions (default `50`).
- `--seed_policy=fixed|varying` (default `varying`).
- `--base_seed=` base seed (default `42`).
- `--seeds=` optional comma-separated seed list.

## Output

Two CSVs are written:

- raw run-level output (`out_raw`)
- summarized output (`out_summary`)

Default output paths are under `/tmp`.

## Example

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/oneoff/bench_oneoff_param.R \
  --fun=npunitest --n=1000 --times=5 --seed_policy=varying \
  --out_raw=/tmp/np_oneoff_npunitest_raw.csv \
  --out_summary=/tmp/np_oneoff_npunitest_summary.csv
```

## Suite Runner

To run all one-off functions in a single pass:

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/oneoff/run_oneoff_suite.R \
  --n_values=100,250,500 --times=50 --seed_policy=varying --base_seed=42 \
  --out_manifest=/tmp/np_oneoff_suite_manifest.csv
```
