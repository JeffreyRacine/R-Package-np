# `nmulti` Core Surface Study

This is the focused study harness for the default-`nmulti` question.

It is intentionally narrower than the broader `nmulti_default_study.R`:

- only core serial estimators are used: `npreg` and `npcdens`,
- only the default-sensitive fixed-bandwidth selectors are studied,
- the design uses one continuous and one unordered categorical regressor,
- the method grid spans `lc`, `ll`, and `lp` with degrees `1` and `2`,
- the outputs emphasize fitted-surface stability, not just raw timings.

## Why This Narrower Harness Is Better

For a default-choice question, breadth can become noise. This harness avoids:

- over-spending CPU on estimators not needed for the decision,
- mixing nearest-neighbor and fixed-bandwidth conclusions too early,
- drawing conclusions from bandwidth-vector movement alone.

Instead it asks:

1. How much time does smaller `nmulti` save?
2. How much do the fitted surfaces move at representative evaluation points?
3. Do those moves matter against a known truth?
4. Does the cross-validation objective move meaningfully?

## Default Run

```sh
Rscript /Users/jracine/Development/np-master/benchmarks/validation/nmulti_core_surface_study.R \
  --n_grid=125,250 \
  --seeds=101,102,103,104,105,106 \
  --nmulti_grid=1,2,5 \
  --out_dir=/tmp/nmulti_core_surface
```

## Outputs

- `nmulti_core_surface_runs.csv`
- `nmulti_core_surface_eval.csv`
- `nmulti_core_surface_summary.csv`
- `nmulti_core_surface_drift.csv`
- `nmulti_core_surface_report.md`
