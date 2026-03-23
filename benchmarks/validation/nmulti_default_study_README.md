# `nmulti` Default Study

This validation harness is aimed at one question: how much do extra
multistarts buy us, and when?

The study is deliberately paired by dataset seed. For each scenario and
sample size, it fits the same synthetic dataset with each requested
`nmulti` value and records:

- bandwidth-selection runtime,
- final bandwidth vector,
- objective diagnostics (`fval`, `ifval`, `num_fval` when available),
- truth-based prediction error on a fixed evaluation design,
- prediction-surface drift relative to the largest tested `nmulti`.

The goal is not to prove a true global optimum. The goal is to see whether
smaller `nmulti` values are practically indistinguishable from a larger
reference value on the same data while saving meaningful time.

## Why This Is Scientifically Better Than Timing Alone

Timing alone can make `nmulti=1` look attractive even if it occasionally
lands on materially different bandwidths. This harness separates three
questions:

1. Does lower `nmulti` save time?
2. Does it materially change the selected bandwidths or fitted surface?
3. Does any change show up in truth-based prediction error?

Only if the answer is "yes, no, no" with enough repetition should a lower
default be taken seriously.

## Presets

- `smoke`
  - scenarios: one representative case each for `npreg`, `npcdens`, and
    `npcdist`
  - `n`: `100,250`
  - seeds: `101:103`
- `full`
  - broader mix of fixed/generalized-NN and LL/LP or LC/LP scenarios
  - `n`: `100,250,500`
  - seeds: `101:120`

## Example Runs

```sh
Rscript /Users/jracine/Development/np-master/benchmarks/validation/nmulti_default_study.R \
  --preset=smoke \
  --nmulti_grid=1,2,5 \
  --out_dir=/tmp/nmulti_default_smoke
```

```sh
Rscript /Users/jracine/Development/np-master/benchmarks/validation/nmulti_default_study.R \
  --preset=full \
  --nmulti_grid=1,2,5 \
  --seeds=101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120 \
  --out_dir=/tmp/nmulti_default_full
```

## Outputs

- `nmulti_default_study_runs.csv`
  - one row per fit
- `nmulti_default_study_eval.csv`
  - one row per evaluation point
- `nmulti_default_study_summary.csv`
  - scenario-level summaries by `n` and `nmulti`
- `nmulti_default_study_drift.csv`
  - seed-matched prediction drift relative to the largest tested `nmulti`
- `nmulti_default_study_report.md`
  - compact markdown report

## Decision Frame

For each estimator family, a lower `nmulti` is a serious default candidate
only if:

- the median bandwidth-selection speedup is meaningful,
- prediction-surface drift versus the reference `nmulti` is negligible,
- truth-based primary error is within 1% of the best tested value for the
  vast majority of seeds,
- no hard scenario shows a long-tail failure pattern hidden by the mean.
