# Delegated numerical contracts

These nine external test blocks replace an undocumented promise of separate
MPI coverage. They deliberately stay outside the default CRAN/testthat lane.
The package implementation is unchanged. The runner loads the installed npRmpi;
set R_LIBS explicitly to the library being qualified. Use a fresh process for
each block, both one and three workers, a bounded external timeout, and an
ownership-aware observer to verify raw exit status and process teardown.

Example, from a source checkout (create the output directory first):

```sh
R_LIBS=/absolute/private/library TMPDIR=/absolute/Development/tmp/scratch \
  /opt/local/bin/mpiexec -n 1 Rscript --vanilla \
  tests/validation/run_delegated_numerical_contract.R 1 3 \
  /absolute/Development/tmp/results/gap01
```

The launcher is host-specific; use the locally qualified MPI executable.
Do not simply remove skip() calls from the default tests: some original blocks
were placeholders without the delegated assertions. This runner has actual
API/numerical assertions, rejects skips, saves per-block results and exits nonzero
on failed assertions. Test inputs/settings are retained from the serial contract.
The two core tests retain their original numerical tolerances.

| Block | Contract |
|---|---|
| 01 | ANN conditional density CVML/CVLS and conditional CDF CVLS, orders4/6/8, degrees0/2, scalar versus accelerated |
| 02 | Unconditional npudistbw Gaussian CVCDF scalar versus accelerated |
| 03 | npregiv starting values and first-state replay |
| 04 | npregiv coordinate-specific multivariate derivatives and derivative weight identity |
| 05 | npregiv categorical kernel/bandwidth scaling controls |
| 06 | npregivderiv output-only evaluation, state/RNG and S3 alignment |
| 07 | npregiv formula/degree/accessor/NA contracts |
| 08 | npregivderiv formula/default smoothing/accessor contracts |
| 09 | Public IV summary formatting on real fitted objects |

Conditional objective evaluations explicitly use collective calls with
force.local=FALSE; workers receive the acceleration option. Numerical parity
alone is not a universal formula/coverage proof. IV state identities and summary
assertions make their narrower API claims. No new statistical simulation is run.

P1 is separate: adaptive_regression_fullrow_blas_installed.R uses a maintainer-
approved absolute 1e-7 budget only for its fixed signed Gaussian6 fixture.
Positive-kernel and exact below-threshold controls are unchanged. See the dated
P1/delegated qualification report in Development for primitive/high-precision
adjudication and retained original failures. This does not loosen package-wide
arithmetic tolerances.
