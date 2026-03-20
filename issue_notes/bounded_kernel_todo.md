# Bounded Continuous Kernel Rollout TODO (`np-npRmpi`)

## Status (2026-03-20)

- Fixed-bandwidth bounded continuous-kernel support is the active supported
  public slice in serial `np-master` and is now mirrored here.
- The current `npRmpi` mirror builds cleanly after the bounded `npcdens cv.ls`
  repair tranche; public bounded support remains `bwtype = "fixed"` only.
- Public bounded continuous-kernel use with
  `bwtype %in% c("generalized_nn", "adaptive_nn")` remains intentionally
  unsupported and still fails fast in the R constructors.
- This file remains the canonical backlog for widening bounded support beyond the
  current fixed-bandwidth slice in MPI-aware routes.

## CRAN Example Policy

- `npRmpi` is submitted to CRAN with heavy examples in `\dontrun{}`.
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

1. Public nonfixed bounded activation
- Evaluate whether bounded continuous kernels can be activated publicly for:
  - `bwtype = "generalized_nn"`
  - `bwtype = "adaptive_nn"`
- Current public policy remains:
  - finite continuous bounds require `bwtype = "fixed"`
- Activation requirements before removing that hard-stop:
  - serial proof already complete,
  - mirrored `npRmpi` session/attach/profile validation,
  - fixed-point objective certification on bounded anchors,
  - end-to-end fit sanity for `cv.ml` and `cv.ls`,
  - explicit MPI route performance evidence.

2. C kernel normalization path
- Verify the bounded normalization path remains correct/active for all supported
  continuous kernels/orders in `src/jksum.c`:
  - Gaussian `2/4/6/8`
  - Epanechnikov `2/4/6/8`
  - rectangular
  - truncated Gaussian `2`
- Keep defaults exact when bounds are effectively unbounded.

3. C API plumbing
- Pass resolved bounds from R into C entry points that evaluate kernels directly:
  - `np_kernelsum`, `np_regression`, `np_density`, `np_density_conditional`
  - BW-selection entry points as needed for CV consistency.
- Update `src/np_init.c` registrations/signatures and all `.Call(...)` call sites.

4. Convolution/CV behavior
- Fixed-bandwidth bounded convolution is now the target mirrored behavior for
  public `npcdens cv.ls`.
- Remaining widening work:
  - confirm bounded convolution semantics for public `generalized_nn`,
  - determine whether bounded `adaptive_nn` convolution is correct and supportable,
  - decide whether additional analytic bounded-convolution kernels beyond
    Gaussian-2 are worth implementing for speed.

5. Denominator micro-optimization (`NZD_pos` in C)
- Expand/verify positive-only guard coverage for any remaining known-nonnegative denominator sites in `jksum.c` beyond the current rollout.
- Add focused pre/post perf artifacts for this micro-optimization slice (fixed and varying seed policy where applicable).

6. Predict/eval policy
- Verify hard-stop behavior is consistent across all user-facing predict/eval paths.
- Add clear error messages for offending variable(s).

7. Testing
- Add/extend tests for:
  - `none`/`range`/`fixed`
  - scalar vs vector bounds
  - invalid bounds (`a >= b`, training support violation)
  - eval support violation
  - finite bounds + non-fixed `bwtype` hard-stop (until activation tranche lands)
  - parity check: `(-Inf, Inf)` equals baseline behavior.

8. Documentation
- Update relevant `.Rd` files for new args and semantics:
  - bw constructors and front-end bw wrappers
  - estimator docs where eval hard-stop applies.
- Add note in package-level docs/changelog about bounded-kernel feature and caveats.
- Be explicit that current public support is `bwtype = "fixed"` only until a
  later nonfixed widening tranche is proven.
- Confirm no `options()`-based interface is required for this feature.

9. Porting to `np-npRmpi`
- After `np-master` stabilizes and tests pass:
  - port shared core changes
  - merge repo-specific R-layer docs/options carefully
  - run smoke benchmark/parity check in MPI path.

10. Conditional density/distribution `regtype` extension (`lc` -> `ll`)
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

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
