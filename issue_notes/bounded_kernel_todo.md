# Bounded Continuous Kernel Rollout TODO (`np-master`)

## Status (2026-02-24)

- Deferred to future modernization tranche by decision.
- This file remains the canonical backlog for bounded convolution/native completion.
- It is not a blocker for concluding the current modernization push.

## CRAN Example Policy

- `np` is submitted to CRAN with heavy examples in `\dontrun{}`.
- Do not treat full `\dontrun{}` example execution as a release gate.
- Run those examples only when explicitly validating behavior after code changes.

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
- Added positive-only denominator guards (`NZD_POS`) in bounded/native hot paths where denominator sign is known nonnegative:
  - kernel marginal denominator guards in convolution CV loops,
  - ridge-adjusted `KWM[0][0]` correction terms,
  - all-large-h CV leave-one-out denominator guard (`1 - hii`).

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
- Update `src/np_init.c` registrations/signatures and all `.Call(...)` call sites.

3. Convolution/CV behavior
- Confirm convolution-dependent objectives (`cv.ls` paths) under finite bounds.
- If bounded convolution is unavailable for a path, define/implement explicit fallback behavior and document it.

4. Denominator micro-optimization (`NZD_pos` in C)
- Expand/verify positive-only guard coverage for any remaining known-nonnegative denominator sites in `jksum.c` beyond the current rollout.
- Add focused pre/post perf artifacts for this micro-optimization slice (fixed and varying seed policy where applicable).

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

9. Conditional density/distribution `regtype` extension (`lc` -> `ll`)
- Add `regtype = c("lc","ll")` to:
  - `conbandwidth`/`npcdens*`
  - `condbandwidth`/`npcdist*`
- Keep default `regtype = "lc"` (current behavior unchanged).
- Thread `regtype` through R option lists and C indices (`headers.h`/`np.c`).
- Implement `ll` branch in conditional estimator core (`src/jksum.c`, `np_kernel_estimate_con_dens_dist_categorical`) rather than only wiring in R.
- Initial guardrails (phase 1):
  - hard-stop `gradients=TRUE` with `regtype="ll"` until ll-gradient formulas are implemented,
  - hard-stop ll bandwidth-selection paths (`cv.ml`/`cv.ls`) until ll CV objectives are validated/implemented.
- Evaluate shape constraints for conditional CDF under ll (`npcdist`):
  - potential finite-sample non-monotonicity and out-of-[0,1] behavior; decide policy (raw vs optional projection).
