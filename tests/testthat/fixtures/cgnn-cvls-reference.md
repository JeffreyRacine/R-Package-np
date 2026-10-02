# Independent conditional GNN CVLS references (2026-09-30)

These two cached values replace pre-CF167 training-radius convolution values.
They do not change the test tolerances or fixed/ANN/CVML/CDF obligations.
Small dynamic criterion oracles remain in test-conditional-gnn-general.R.

For each held-out identity i, recompute predictor radii on the n-1 donors.
Form the literal Gaussian-order-8 weighted raw polynomial basis
(1,x1,x2,x1^2,x1*x2,x2^2). Solve the deleted 6x6 Gram system at 80 digits,
then obtain signed donor weights. Bernstein and raw representations span the
same GLP space; the oracle does not use the package's basis/solve routines.

For every real response query q, recompute its ky-th distance on those same
n-1 response identities. Evaluate the Gaussian-order-2 response density with
that query radius. Integrate the mean squared density over the whole real line,
partitioned at all response values and pairwise midpoints. The score is
2*mean(f_minus_i(y_i|x_i)) - integral(mean(f_minus_i(q|x_i)^2)) dq.
No package radius, smoother, objective or integration function is used.

| Fixture in test | Independent score | Quadrature estimated absolute error |
|---|---:|---:|
| n54, seed2026080101, ky21/kx24,25 | 0.94383362371376034 | 3.20e-13 |
| n96, seed2026072701, ky31/kx34,29 | -30.148105184524063 | 6.79e-13 |

Reference environment: NumPy2.4.6, SciPy1.17.1 (QUADPACK), mpmath1.3.0.
The Gram solves use 80-digit arithmetic on double-precision input/kernel
weights; this is not an 80-digit evaluation of the complete criterion.
Quadrature error estimates exclude floating-point kernel/input error.
Largest Gram condition estimates: 3.76e5 and 3.61e7 respectively.
Row-sum errors: 5.56e-16 and 7.11e-15. Native absolute discrepancies after
independent computation: 4.95e-14 and 1.90e-10, within unchanged 2e-8 tolerance.

Reproduction/evidence archive:
Development/tmp/r24_actionable_repairs_20260930/harness/t1_oracle.py,
harness/t1_fixtures.R, t1-inputs, attempts/t1-oracle54-r1,
attempts/t1-oracle96-r1 and attempts/t1-compare-r1.
The 17-digit CSV input SHA256 values are:

- n54: 8c1fb0caa4c4cbafd82147cce5ebd9c253350f42995fb850386193d284f85ac6
- n96: 813c400360a21ca0688093e13f3ccb2c08ded4de5a1b25b6639e662a65442903

Native agreement is a subsequent comparison, not how the references were set.

## Bounded explanatory beta reference (R25 repair, 2026-10-01)

The beta-X order-4, LP(2,2), seed-2026080112, n=54, k=(21,24,25)
fixture in `test-conditional-density-cvls-delete-one-contract.R` has score
`0.54294422783929308`. Independent beta density weights, signed deleted WLS,
and response query-radius squared-density integration over every midpoint
interval plus both infinite tails give this value at 32 Gauss nodes. The
20-node result agrees to the existing test tolerance. The prior training-radius
convolution is not an oracle for this GNN case.

Retained derivation, reviewed oracle and raw results:
`tmp/r25_release_repairs_20261001_r1/harness/bounded_x_oracle.R` and
`attempts/bounded-x-source-2-{20,32}/raw.log` under the Development workspace.
