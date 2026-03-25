# Bounded Continuous Kernel Rollout Status (`np-master`)

## Status (2026-03-25)

Bounded continuous-kernel support in serial `np` is now publicly enabled for
all currently certified public families under:

- `bwtype = "fixed"`
- `bwtype = "generalized_nn"`
- `bwtype = "adaptive_nn"`

Certified serial public surface:

- core families:
  - `npudens*`
  - `npudist*`
  - `npreg*`
  - `npcdens*`
  - `npcdist*`
- semiparametric/index families:
  - `npplreg*`
  - `npindex*`
  - `npscoef*`

This file is no longer a “fixed-only activation” backlog. It is now the
canonical serial status note plus the remaining bounded-kernel backlog for work
that is still incomplete after the nonfixed widening tranches.

## CRAN Example Policy

- `np` is submitted to CRAN with heavy examples in `\\dontrun{}`.
- Do not treat full `\\dontrun{}` example execution as a release gate.
- Run those examples only when explicitly validating behavior after code changes.

## Done (Serial Public Rollout)

1. R-side plumbing/checks
- Added continuous bound args to constructors and persisted into bw objects:
  - `ckerbound/ckerlb/ckerub`
  - `cxkerbound/cxkerlb/cxkerub`
  - `cykerbound/cykerlb/cykerub`
- Added bound resolution helper (`none`/`range`/`fixed`) with strict validation:
  - `/Users/jracine/Development/np-master/R/util.R` (`npKernelBoundsResolve`)
- Added eval-time hard-stop helper for out-of-bounds prediction/evaluation:
  - `/Users/jracine/Development/np-master/R/util.R` (`npKernelBoundsCheckEval`)
- Applied eval hard-stop checks in:
  - `npudens`, `npudist`, `npreg`, `npcdens`, `npcdist`, `npksum`

2. Public bounded nonfixed activation
- Core serial public support is enabled for bounded `generalized_nn` on:
  - `npudens*`, `npudist*`, `npreg*`, `npcdens*`, `npcdist*`
- Core serial public support is enabled for bounded `adaptive_nn` on:
  - `npudens*`, `npudist*`, `npreg*`, `npcdens*`, `npcdist*`
- Serial semiparametric/index bounded nonfixed support is enabled for:
  - `npplreg*`
  - `npindex*`
  - `npscoef*`

3. Bounded convolution/CV repair status
- Fixed-bandwidth bounded convolution is implemented/used for public
  `npcdens cv.ls`.
- Serial bounded nonfixed public routes were widened only after proof bundles
  and tarball-first `R CMD check --as-cran` gates remained green.

4. Native hot-path hardening
- Added positive-only denominator guards (`NZD_POS`) in bounded/native hot paths
  where denominator sign is known nonnegative:
  - kernel marginal denominator guards in convolution CV loops
  - ridge-adjusted `KWM[0][0]` correction terms
  - all-large-h CV leave-one-out denominator guard (`1 - hii`)

## Remaining (Serial Backlog)

1. C kernel normalization path audit
- Verify the bounded normalization path remains correct/active for all supported
  continuous kernels/orders in `src/jksum.c`:
  - Gaussian `2/4/6/8`
  - Epanechnikov `2/4/6/8`
  - rectangular
  - truncated Gaussian `2`
- Keep defaults exact when bounds are effectively unbounded.

2. C API plumbing review
- Re-audit that resolved bounds are passed into all native entry points that
  evaluate kernels directly:
  - `np_kernelsum`
  - `np_regression`
  - `np_density`
  - `np_density_conditional`
  - BW-selection entry points as needed for CV consistency
- Reconfirm `src/np_init.c` registrations/signatures and `.Call(...)` sites.

3. Convolution/CV completeness
- Bounded nonfixed support is now public, but convolution remains the most
  specialized bounded sub-area.
- Remaining follow-up questions:
  - confirm bounded convolution semantics for the currently widened nonfixed
    public slice where convolution-based selectors are involved
  - decide whether additional analytic bounded-convolution kernels beyond
    Gaussian-2 are worth implementing for speed

4. Predict/eval policy audit
- Reconfirm hard-stop behavior is consistent across all user-facing
  predict/eval paths.
- Keep variable-specific error diagnostics crisp.

5. Documentation
- Update relevant `.Rd` files for the now-widened bounded semantics:
  - bw constructors and front-end bw wrappers
  - estimator docs where eval hard-stop applies
- Add a release-facing note summarizing the currently certified bounded surface.

6. Porting status
- `np-master` remains the source of truth for shared bounded-kernel behavior.
- Keep mirrored `npRmpi` status tracked separately in:
  - `/Users/jracine/Development/np-npRmpi/issue_notes/bounded_kernel_todo.md`

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
