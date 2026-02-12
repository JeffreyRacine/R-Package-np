# Optimizations Log

## 2026-02-12

### Summary
Implemented serial performance optimizations in `src/jksum.c` targeting high-frequency kernel evaluation paths.

### Code Changes
- Unordered kernel fast path:
  - Cache `same/different` kernel values once and reuse in `np_ukernelv` and `np_p_ukernelv`.
- Ordered kernel fast path:
  - Cache `lambda^d` values and use direct ordered-kernel evaluation for kernels `0/1/2` in `np_okernelv` and `np_p_okernelv`.
- Tree-path early exit:
  - If tree search finds no continuous-support interactions for a query point, skip downstream kernel work immediately.
- Persistent discrete profile cache:
  - Build discrete profile ids once per `kernel_weighted_sum_np` call and reuse across all query points.
  - Apply cached profile products in both tree and non-tree paths under conservative activation conditions.

### Runtime Summary (pre vs post)
- Real-world `wage1` example (`npreg`, `ckertype="epanechnikov"`, `times=8`):
  - `np.tree=FALSE`: `10927.178 ms` -> `8466.990 ms` (mean, `+22.5%` speedup)
  - `np.tree=TRUE`: `6033.697 ms` -> `5749.531 ms` (mean, `+4.7%` speedup)
- Non-tree categorical-only benchmark (`npreg`, `y~z1+z2+z3`, `n=4000`, `times=6`):
  - `31371.77 ms` -> `20732.23 ms` (mean, `+33.9%` speedup)
- Non-tree mixed benchmark (`npreg`, `y~x1+x2+z1+z2`, `n=2000`, `times=6`):
  - `4677.014 ms` -> `3921.773 ms` (mean, `+16.1%` speedup)
- Tree mixed benchmark (`npreg`, `n=2000`, `times=6`):
  - `1189.502 ms` -> `1099.57 ms` (mean, `+7.6%` speedup)

### Numerical Stability
- Real-world `wage1`, pre vs post:
  - non-tree: `max |Δ bw| = 2.759011e-10`, `max |Δ fitted| = 4.730198e-10`
  - tree: `max |Δ bw| = 1.603712e-04`, `max |Δ fitted| = 1.500707e-05`
- Additional synthetic checks show only very small floating-point-level differences from operation reordering.
