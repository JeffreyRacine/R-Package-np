# Issue #44: `Bvcov` in `npindex.sibandwidth`

**Status:** Resolved for 0.70-5

## Defect

For an Ichimura model with more than one free index coefficient, the previous
expressions multiplied an observation-length vector directly by an
`r`-by-`n` score matrix. R's column-major recycling therefore applied
link-gradient and residual weights to the wrong score elements when `r > 1`.
The special case `r = 1` was already correct.

## Resolution

For observation `i`, let `z_i` be the conditionally centered free-predictor
vector, `g_i` the link gradient, and `u_i` the residual. The implementation now
constructs score column `s_i = g_i z_i` and uses

- `A = sum_i s_i s_i'`;
- `B = sum_i u_i^2 s_i s_i'`;
- `A^-1 B A^-1` for the free-coefficient covariance.

The gradient and residual vectors are applied observation by observation with
column-wise `sweep()`. The existing matrix multiplication and direct Cholesky
contract are retained. No ridge, generalized inverse, silent fallback,
training-data retention, optimizer change, or public default was introduced.

## Validation

- the two-predictor case is bit-for-bit unchanged;
- installed three-, four-, and fifteen-predictor fits match an independent
  literal observation-loop oracle;
- beta, bandwidth, objective, function-evaluation count, fitted values,
  residuals, and gradients are exactly unchanged;
- the fifteen-predictor information matrix remains full rank in the retained
  provenance example;
- deliberately rank-deficient information matrices continue to fail Cholesky
  rather than being masked;
- the corresponding `npRmpi` implementation is exactly equivalent with one
  and three slaves and changes no MPI payload.

The GitHub issue closeout records the exact resolving commit hashes.
