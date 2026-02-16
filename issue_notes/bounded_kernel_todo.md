# Bounded Continuous Kernel Rollout TODO (`np-master`)

## Done (R-side plumbing/checks)

- Added continuous bound args to constructors and persisted into bw objects:
  - `ckerbound/ckerlb/ckerub`
  - `cxkerbound/cxkerlb/cxkerub`
  - `cykerbound/cykerlb/cykerub`
- Added bound resolution helper (`none`/`range`/`fixed`) with strict validation:
  - `/Users/jracine/Development/np-master/R/util.R` (`npKernelBoundsResolve`)
- Added eval-time hard-stop helper for out-of-bounds prediction/evaluation:
  - `/Users/jracine/Development/np-master/R/util.R` (`npKernelBoundsCheckEval`)
- Enforced compatibility rule:
  - finite bounds require `bwtype = "fixed"` in constructors.
- Applied eval hard-stop checks in:
  - `npudens`, `npudist`, `npreg`, `npcdens`, `npcdist`, `npksum`.
- Completed semiparametric propagation in:
  - `np.plregression.bw.R`, `np.singleindex.bw.R`.

## Remaining (core implementation)

1. C kernel normalization path
- Implement bounded normalization for continuous kernels in `src/jksum.c`:
  - `K((x-x_i)/h) / (G((b-x)/h)-G((a-x)/h))`
- Ensure the path is active in all relevant operator modes where intended.
- Keep defaults exact when bounds are effectively unbounded.

2. C API plumbing
- Pass resolved bounds from R into C entry points that evaluate kernels directly:
  - `np_kernelsum`, `np_regression`, `np_density`, `np_density_conditional`
  - BW-selection entry points as needed for CV consistency.
- Update `src/np_init.c` registrations/signatures and all `.C(...)` call sites.

3. Convolution/CV behavior
- Confirm convolution-dependent objectives (`cv.ls` paths) under finite bounds.
- If bounded convolution is unavailable for a path, define/implement explicit fallback behavior and document it.

4. Denominator micro-optimization (`NZD_pos` in C)
- Add/use a positive-only fast denominator guard for known-positive quantities in `jksum.c` (C analogue of R-side `NZD_pos` usage in `util.R`).
- Apply this in hot paths where denominator sign is known positive (including bounded-kernel normalization), instead of generic signed guards.
- Benchmark before/after because this path is called very frequently.

5. Predict/eval policy
- Verify hard-stop behavior is consistent across all user-facing predict/eval paths.
- Add clear error messages for offending variable(s).

6. Testing
- Add/extend tests for:
  - `none`/`range`/`fixed`
  - scalar vs vector bounds
  - invalid bounds (`a >= b`, training support violation)
  - eval support violation
  - finite bounds + non-fixed `bwtype` hard-stop
  - parity check: `(-Inf, Inf)` equals baseline behavior.

7. Documentation
- Update relevant `.Rd` files for new args and semantics:
  - bw constructors and front-end bw wrappers
  - estimator docs where eval hard-stop applies.
- Add note in package-level docs/changelog about bounded-kernel feature and caveats.
- Confirm no `options()`-based interface is required for this feature.

8. Porting to `np-npRmpi`
- After `np-master` stabilizes and tests pass:
  - port shared core changes
  - merge repo-specific R-layer docs/options carefully
  - run smoke benchmark/parity check in MPI path.
