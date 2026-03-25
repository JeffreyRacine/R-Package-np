# Bounded Continuous Kernel Rollout Status (`np-npRmpi`)

## Status (2026-03-25)

Bounded continuous-kernel support in `npRmpi` is now substantially wider than
the older fixed-only MPI mirror.

Certified public `npRmpi` surface:

1. Bounded `fixed`
- enabled on the certified core public routes

2. Bounded `generalized_nn`
- enabled on:
  - `npudens*`
  - `npudist*`
  - `npreg*`
  - `npcdens*`
  - `npcdist*`
  - `npplreg*`
  - `npindex*`
  - `npscoef*`

3. Bounded `adaptive_nn`
- enabled on:
  - `npudens*`
  - `npreg*`
  - `npcdens*`
  - `npplreg*`
  - `npindex*`
  - `npscoef*`
- enabled on fit-only distribution routes with precomputed bandwidth objects:
  - `npudist()`
  - `npcdist()`
- intentionally still blocked at selector level:
  - `npudistbw()`
  - `npcdistbw()`

Release exception retained at closeout (2026-03-25):
- selector-level bounded `adaptive_nn` remains unsupported for
  `npudistbw()` / `npcdistbw()` in `npRmpi`
- keep the explicit public errors in place until a deeper dedicated
  MPI-selector campaign proves those routes safe

This file is now the canonical `npRmpi` bounded-kernel status note plus the
remaining MPI-aware backlog. The old statement “bounded public support is
`bwtype = \"fixed\"` only” is no longer true.

## CRAN Example Policy

- `npRmpi` is submitted to CRAN with heavy examples in `\\dontrun{}`.
- Do not treat full `\\dontrun{}` example execution as a release gate.
- Run those examples only when explicitly validating behavior after code changes.

## Done (MPI-Aware Rollout)

1. R-side plumbing/checks
- Added continuous bound args to constructors and persisted into bw objects:
  - `ckerbound/ckerlb/ckerub`
  - `cxkerbound/cxkerlb/cxkerub`
  - `cykerbound/cykerlb/cykerub`
- Added bound resolution helper (`none`/`range`/`fixed`) with strict validation.
- Added eval-time hard-stop helper for out-of-bounds prediction/evaluation.
- Applied eval hard-stop checks in:
  - `npudens`, `npudist`, `npreg`, `npcdens`, `npcdist`, `npksum`

2. Public bounded nonfixed activation
- Bounded `generalized_nn` is mirrored across the certified core plus
  semiparametric/index surface.
- Bounded `adaptive_nn` is mirrored on the certified non-distribution core plus
  semiparametric/index surface.
- Distribution-family bounded `adaptive_nn` fit-only routes are enabled when the
  user supplies precomputed bandwidth objects.

3. Route validation already completed for the widened kept surface
- certified widening tranches were validated in:
  - session/spawn
  - attach
  - profile/manual-broadcast
- final tarball-first `R CMD check --as-cran` gates remained green after the
  accepted keeps.

4. Native hot-path hardening
- Added positive-only denominator guards (`NZD_POS`) in bounded/native hot paths
  where denominator sign is known nonnegative:
  - kernel marginal denominator guards in convolution CV loops
  - ridge-adjusted `KWM[0][0]` correction terms
  - all-large-h CV leave-one-out denominator guard (`1 - hii`)

## Remaining (`npRmpi` Backlog)

1. Deferred selector gap
- Remaining bounded `adaptive_nn` public gap is narrow and explicit:
  - `npudistbw()`
  - `npcdistbw()`
- Keep those selectors blocked until a dedicated MPI selector-routing tranche
  proves them safe on:
  - session/spawn
  - attach
  - profile/manual-broadcast
- Treat this as an MPI selector/orchestration problem, not a broad bounded
  rollout problem.

2. C kernel normalization path audit
- Verify the bounded normalization path remains correct/active for all supported
  continuous kernels/orders in `src/jksum.c`:
  - Gaussian `2/4/6/8`
  - Epanechnikov `2/4/6/8`
  - rectangular
  - truncated Gaussian `2`
- Keep defaults exact when bounds are effectively unbounded.

3. C API plumbing review
- Re-audit that resolved bounds are passed into all native entry points that
  evaluate kernels directly:
  - `np_kernelsum`
  - `np_regression`
  - `np_density`
  - `np_density_conditional`
  - BW-selection entry points as needed for CV consistency

4. Convolution/CV completeness
- Fixed-bandwidth bounded convolution is already mirrored for public
  `npcdens cv.ls`.
- Remaining follow-up questions:
  - confirm bounded convolution semantics for the widened nonfixed public slice
    where convolution-based selectors are involved
  - determine whether any remaining selector-specific bounded convolution work
    is necessary beyond the current certified keep

5. Predict/eval policy audit
- Reconfirm hard-stop behavior is consistent across all user-facing
  predict/eval paths.
- Keep variable-specific error diagnostics crisp.

6. Documentation
- Update relevant `.Rd` files for the now-widened bounded semantics:
  - bw constructors and front-end bw wrappers
  - estimator docs where eval hard-stop applies
- Add a release-facing note summarizing the current `npRmpi` bounded surface,
  including the explicit deferred selector gap for bounded `adaptive_nn`
  `npudistbw()` / `npcdistbw()`.

7. Conditional density/distribution `regtype` extension (`lc` -> `ll`)
- Add `regtype = c("lc","ll")` to:
  - `conbandwidth` / `npcdens*`
  - `condbandwidth` / `npcdist*`
- Keep default `regtype = "lc"` unchanged.
- Thread `regtype` through R option lists and C indices (`headers.h` / `np.c`).
- Implement `ll` branch in conditional estimator core (`src/jksum.c`,
  `np_kernel_estimate_con_dens_dist_categorical`) rather than only wiring in R.
- Initial guardrails:
  - hard-stop `gradients=TRUE` with `regtype="ll"` until ll-gradient formulas
    are implemented
  - hard-stop ll bandwidth-selection paths (`cv.ml` / `cv.ls`) until ll CV
    objectives are validated/implemented
- Evaluate shape constraints for conditional CDF under ll (`npcdist`):
  - possible finite-sample non-monotonicity and out-of-`[0,1]` behavior; decide
    policy (raw vs optional projection)

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`,
   `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`,
   `bernstein.basis`, kernels, and bounds) must be propagated and used by the
   canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain
   semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant
   runtime overhead once canonical behavior exists.
