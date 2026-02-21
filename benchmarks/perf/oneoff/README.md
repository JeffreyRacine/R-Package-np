# One-Off Benchmarks (npRmpi)

MPI work-alike benchmark entry-point for one-off `npRmpi` methods where we do not vary kernel/method options.

This `benchmarks/perf/oneoff` location is canonical. The top-level `benchmarks/oneoff` directory is legacy and retained for compatibility.

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
- `--nslaves=` (or `--rslaves=`) number of slaves (default `1`).
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
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/perf/oneoff/bench_oneoff_param_nprmpi.R \
  --fun=npunitest --n=1000 --nslaves=1 --times=5 --seed_policy=varying \
  --out_raw=/tmp/nprmpi_oneoff_npunitest_raw.csv \
  --out_summary=/tmp/nprmpi_oneoff_npunitest_summary.csv
```

## Suite Runner

To run all one-off functions in a single pass:

```bash
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/perf/oneoff/run_oneoff_suite.R \
  --n_values=100,250,500 --nslaves=1 --iface=en0 --times=50 --seed_policy=varying --base_seed=42 \
  --out_manifest=/tmp/nprmpi_oneoff_suite_manifest.csv
```
