# FFT Acceleration Plan for `np-master` (Tranche-Oriented)

Date: March 2, 2026  
Scope: `np-master` only (serial source of truth)  
Status: planning/specification only (no implementation in this document)

## 0. Development Isolation Policy

Implementation work does not begin on local `master`.

Execution policy for the FFT spike:

1. Create and use a separate throwaway worktree + branch for all FFT experimentation:
- branch: `codex/spike-fft-tranche1`
- worktree path: `/Users/jracine/Development/np-master-fft-spike`
2. Keep `/Users/jracine/Development/np-master` available for ongoing minor `master` work.
3. If spike fails, remove worktree and delete branch with no impact on `master`.

## 1. Objective and Rationale

Provide a practical FFT-accelerated path for low-dimensional, fixed-bandwidth workflows that are common in exploratory analysis, while preserving estimator intent and API semantics.

Primary goal:

- Deliver material speedups in 1D/2D continuous fixed-bandwidth cases for estimation and CV, with hard-stop eligibility checks when `fft=TRUE` is requested and conditions are not met.

Secondary goal:

- Use this as shared infrastructure for downstream methods that depend on repeated kernel sums, starting with `npksum` and then propagating to estimators and selectors.

Existence proof and precedent:

- `stats::density.default` in base R uses a binned + FFT convolution strategy (`C_BinDist` + FFT) with power-of-two grids and default `n = 512`. This supports conceptual feasibility for the core univariate density setting.

## 2. Non-Negotiable Design Principles

1. No silent fallback when user explicitly requests FFT.
- If `fft=TRUE` and path is ineligible, hard `stop()` with explicit reason.

2. `np-master` behavioral source of truth.
- Any later `np-npRmpi` work will be surgical parity porting, not direct copy.

3. Preserve existing defaults.
- `fft=FALSE` default everywhere in tranche-1.

4. Numerical integrity first.
- Acceptance is based on bandwidth argmin displacement and estimator parity, not objective value parity alone.

5. Clear auditability.
- Surface diagnostics for whether FFT was used and why not.

## 3. Scope by Tranche

## 3.1 Tranche-1 (Implement)

Eligibility-limited FFT acceleration for:

1. Unconditional density: `npudens`, `npudensbw` (`cv.ls`, `cv.ml`) in 1D/2D continuous fixed-bw cases.
2. Unconditional distribution: `npudist`, `npudistbw` (`cv.cdf`) in same continuous setting where convolution structure is usable.
3. Conditional density/distribution:
- `1y + 1x` (2D joint continuous block) only for tranche-1.
- Enforce shared-grid consistency for ratio operations.
4. Conditional quantile:
- `npqreg` for `1y + 1x` only, via conditional CDF backbone with explicit inversion validation.
5. Regression: `npreg`, `npregbw` where kernel sums map to eligible low-d continuous fixed-bw calls.
6. Shared primitive: `npksum` FFT backend under strict eligibility gates.
7. Immediate inherited beneficiaries where applicable:
- `npindex*`, `npscoef*`, `npregiv*` inner `npksum` calls when those calls meet FFT eligibility.

## 3.2 Tranche-1 Explicitly Out of Scope (Hard Stop when `fft=TRUE`)

1. Continuous dimension above 2 in FFT-smoothed block.
2. `bwtype != "fixed"` (`generalized_nn`, `adaptive_nn`).
3. Cases requiring unsupported kernels/boundary policies (until explicitly supported).
4. Workloads that violate memory/grid ceilings.
5. Conditional density with `1y + 2x` (3D joint requirement).

## 3.3 Tranche-2 Candidates (Defer)

1. Higher-dimensional FFT blocks.
2. Expanded kernel/boundary support matrix.
3. More aggressive semiparametric loop-level optimization.
4. Possible MPI-aware FFT route in `np-npRmpi`.

## 4. Integration Order (Dependency-First)

1. `npksum` backend expansion in C.
- Add FFT backend as an alternative execution route under existing kernel-sum architecture.

2. Estimation paths.
- Wire eligible estimation calls to FFT-capable kernel-sum backend.

3. Bandwidth objective paths.
- `cv.ls` first, then `cv.ml` with explicit LOO correction policy.

4. Semiparametric pass-through validation.
- Validate real gains for `npindex`/`npscoef` inner sums while documenting outer-loop limits.

## 5. API and Control Surface

Add `fft=FALSE` and optional `fft.control` to relevant `np*` and `np*bw` entry points.

`fft.control` proposed fields (subject to final lock):

1. `grid.max.per.axis` (upper bound; power-of-two target).
2. `grid.min.per.axis`.
3. `memory.max.bytes`.
4. `interpolation` (`linear` in tranche-1).
5. `verbose` (emit diagnostics).

Behavior contract:

1. `fft=FALSE`: existing behavior.
2. `fft=TRUE`: use FFT only if eligible; otherwise hard-stop with reason.
3. `fft="auto"` is not in tranche-1 (avoid silent route changes).
4. `strict` toggle is not exposed in tranche-1; hard-stop semantics are fixed when `fft=TRUE`.

## 6. Eligibility Gate Specification (Tranche-1)

A call is FFT-eligible only if all conditions hold:

1. `bwtype == "fixed"`.
2. Requested kernel is in tranche-1 supported kernel set.
3. Boundary mode is supported by tranche-1 boundary policy.
4. Continuous FFT block dimension is 1 or 2.
5. Requested outputs (gradients/stderr/etc.) are supported by FFT route for that method.
6. Categorical expansion (if any) stays below configured cell cap.
7. Dynamic grid sizing returns feasible grid under memory cap.

If any condition fails and `fft=TRUE`, stop with a message naming:

1. violated gate,
2. observed value,
3. threshold/requirement,
4. suggested action.

## 7. Kernel and Boundary Policy (Tranche-1)

Current intent:

1. Tranche-1 continuous kernels are restricted to order-2:
- Gaussian (order 2)
- Epanechnikov (order 2)
- Uniform
2. Higher-order continuous kernels (orders 4/6/8) are deferred to tranche-2.
3. Uniform remains tranche-1 only if anti-ringing parity tests pass at locked tolerances.

4. Bounded-kernel handling:
- Either explicitly implement correct boundary normalization on grid or hard-stop bounded modes in tranche-1.
- No partial/implicit treatment.

5. Mixed data:
- FFT only on continuous axes.
- Unordered/ordered effects remain multiplicative/sliced.

## 8. Dynamic Grid Policy (No Hardcoded 512)

`512` is a sensible default anchor, not a fixed rule.

Grid policy requirements:

1. Grid spacing must be tied to bandwidth magnitude and kernel support.
2. Grid size must be power-of-two per axis.
3. Domain must include explicit padding to avoid circular aliasing.
4. If required grid exceeds `memory.max.bytes` or `grid.max.per.axis`, hard-stop under `fft=TRUE`.
5. `grid.min.per.axis` default is `64`; values below `64` are permitted only with an explicit diagnostic warning.

Preflight must report:

1. selected grid shape,
2. estimated memory,
3. padding and interpolation mode.

## 9. CV Objective Treatment

## 9.1 `cv.ls`

1. Preserve decomposition structure (convolution term and leave-one-out term as defined in theory/canonical implementation).
2. In tranche-1, compute both `cv.ls` terms on the same FFT/grid approximation path (no mixed analytic+grid term computation).
3. Ensure both terms are computed on compatible discretization to avoid structural argmin drift.
4. Validate argmin stability against exact pairwise route.

## 9.2 `cv.ml`

1. Must be true LOO objective semantics, not pseudo-LOO.
2. Apply explicit self-contribution correction to obtain `f_{-i}(x_i)` before log evaluation.
3. Keep denominator floor/log handling consistent with existing `DBL_MIN` style guardrails.
4. Remaining approximation source should be only binning/interpolation, documented as such.
5. `n_eff<=1` in any active slice/profile is a hard stop under `fft=TRUE` for `cv.ml`.

## 9.3 `cv.cdf`

1. `cv.cdf` semantics are treated as first-class in tranche-1 (not inherited by analogy).
2. Preserve current estimator objective structure and explicit leave-one-out handling semantics:
- `F_{-i} = (n_eff/(n_eff-1))*F_i - G0/(n_eff-1)` in 1D, with product-form `G0` for multi-axis FFT blocks.
3. Apply the same selector-stability gates as `cv.ls`/`cv.ml` plus estimator-level parity checks on CDF values.
4. `n_eff<=1` in any active slice/profile is a hard stop under `fft=TRUE` for `cv.cdf`.

## 9.4 `npqreg` Inversion Policy

1. Invert conditional CDF using bracketed bisection (deterministic) in tranche-1.
2. Apply monotonicity enforcement on grid CDF (running-max within support) before inversion.
3. If inversion fails to converge within tolerance/iteration limits, hard-stop under `fft=TRUE` with a deterministic reason.

## 10. Conditional Density/Distribution Ratio Consistency

For cases such as `f(y|x)=f(y,x)/f(x)`:

1. Joint and marginal evaluations must be computed on compatible grids.
2. Ratio evaluation must avoid denominator-near-zero instability.
3. Domain trimming/guard policy must be explicit and deterministic.

This is mandatory even for tranche-1 `1y + 1x`.

## 11. Semiparametric Expectation Management

`npindex` and `npscoef` will not generally get uniform end-to-end speedup from FFT.

What can improve:

1. Inner `npksum` evaluations that are FFT-eligible.

What may still dominate:

1. Outer `optim()` multistart loops.
2. Backfitting/iteration logic.

Documentation and release notes must set this expectation clearly.

## 12. Build/Dependency Strategy

Tranche-1 decision direction:

1. Avoid external FFTW dependency in tranche-1 (portability/licensing/CRAN friction).
2. Use backend abstraction with tranche-1 baseline implementation via `stats::fft`/`stats::mvfft` (mixed-radix, no new dependency).
3. Defer custom internal C radix-2 FFT to tranche-2 unless tranche-1 performance gates fail.

Implication:

1. Lower operational risk for CRAN and reproducibility.

## 13. Diagnostics and User Transparency

Every FFT-eligible run should retain metadata in returned objects where feasible:

1. `used_fft` flag.
2. `fft_grid`.
3. `fft_interpolation`.
4. `fft_memory_estimate`.
5. `fft_reason` (if not used).
6. `fft_gate_failed` (if `fft=TRUE` and ineligible).

No hidden execution route changes.

Metadata placement contract:

1. Store `fft_*` metadata on both bandwidth objects (`np*bw`) and estimator objects (`np*`).
2. Estimator objects inherit `fft_*` metadata from supplied bw objects and append run-time fields (`used_fft`, realized grid, memory estimate).

## 14. Validation Protocol (Acceptance Gates)

A tranche does not pass unless all of the following pass on representative test batteries.

Numerical parity and selector stability:

1. Bandwidth argmin displacement:
- median relative error <= 5%
- 90th percentile <= 10%
2. Objective consistency:
- no systematic directional bias in selected bandwidths.
3. Downstream estimator parity at selected bandwidths:
- evaluate estimator at `h_fft` and `h_exact` on matched evaluation grids,
- require `L_inf` relative fitted-value error <= 2% in central support and <= 5% globally.
4. Sign-bias checks are run per simulation setting (not pooled across settings).

Performance:

1. Small and medium workload scales.
2. Report mean and median speedup.
3. Record crossover zones where FFT is slower.
4. Preserve raw artifacts in dated `/tmp` paths.

Scope matrix:

1. Methods: `npudens*`, `npudist*`, `npcdens*`, `npcdist*`, `npqreg*`, `npreg*`, `npksum`.
2. Kernels in tranche-1 support set.
3. Continuous dimensions 1 and 2.
4. Relevant mixed-data slices where applicable.

Locked validation sizes:

1. Small: `n in {50, 100}`
2. Medium: `n in {250, 500}`
3. Large: `n in {1000, 2500}`
4. Crossover documentation: include at least one very-large run (`n in {5000, 10000}` where feasible).

## 15. Additional Risks We Have Already Identified

1. Boundary normalization mismatch.
2. `bwscaling` representation mismatch.
3. NN bandwidth incompatibility.
4. Padding/aliasing mistakes.
5. Interpolation-induced objective roughness.
6. `cv.ml` LOO/log-floor instability.
7. `cv.ls` term inconsistency.
8. Conditional ratio instability.
9. Semiparametric overpromising.
10. Unsupported gradients/stderr requests.
11. Fast/fallback counters semantics drift.
12. Legacy/new C path asymmetry.
13. Global C state contamination.
14. Memory-thrash at 2D grid sizes in repeated CV loops.
15. NA/factor-level handling regressions.
16. Insufficient diagnostic visibility.

## 16. Pre-Mortem Checklist

## 16.1 Must Decide Before Code (Blocking)

1. Approximation acceptance contract:
- exact argmin displacement thresholds and parity gates.
2. FFT backend choice:
- self-contained internal implementation (no external FFTW) yes/no lock.
3. `cv.ml` semantics:
- true LOO correction formula and log-floor behavior lock.
4. Tranche-1 supported kernel list:
- exact inclusion/exclusion.
5. Boundary policy:
- supported boundary modes vs hard-stop list.
6. Dynamic grid policy:
- formula, padding rule, power-of-two strategy.
7. Memory policy:
- hard limits and fail behavior.
8. Interpolation policy:
- tranche-1 method (linear) and no hidden alternatives.
9. Eligibility error contract:
- hard-stop message schema and required fields.
10. Conditional ratio guard policy:
- denominator floor/trimming and diagnostics.
11. Output compatibility policy:
- whether FFT route supports gradients/stderr per method in tranche-1.
12. Validation battery and pass/fail thresholds:
- simulation design, locked `n` grids (small/medium/large/crossover), and required artifact outputs.
13. `cv.cdf` decision:
- explicit semantics lock and validation gates for distribution bandwidth selection.
14. Bandwidth-unit resolution (`bwscaling`):
- resolve all FFT grid construction in data units, with explicit unit conversion contract.
15. `cv.ml`/`cv.cdf` sparse-slice policy:
- explicit handling for `n_eff<=1`.
16. Eligibility gate ordering:
- explicit first-fail order for deterministic user messages.
17. CV optimizer loop architecture:
- per-call vs pre-binned workspace policy.
18. Backend call-stack contract:
- where binning/FFT/interpolation/LOO run (C vs R) under `stats::fft` backend.
19. `grid.min.per.axis` default:
- locked default and warning behavior for low values.
20. `cv.ls` integral-term computation policy:
- consistent approximation strategy (no mixed-mode ambiguity).
21. `npqreg` inversion policy:
- algorithm, monotonicity enforcement, convergence/failure contract.
22. `npudist` CDF recovery policy:
- cumulative integration domain, normalization, and fail guard.
23. `na.action` and FFT preflight NA policy:
- explicit contract for NA-free continuous blocks.
24. C workspace state policy:
- no global/static mutable state.
25. FFT metadata class-placement policy:
- bw and estimator inheritance/merge behavior.
## 16.2 Can Defer to Tranche-2 (Non-Blocking for Tranche-1)

1. Dimensions above 2 for FFT block.
2. Expanded kernel families and higher-order kernels beyond tranche-1 set.
3. Full bounded-kernel coverage if initially hard-stopped.
4. `fft="auto"` policy and heuristics.
5. Advanced interpolation variants.
6. Deeper semiparametric outer-loop optimization redesign.
7. MPI FFT strategy in `np-npRmpi`.
8. Additional UI/reporting niceties beyond required diagnostics.
9. Derivative-kernel FFT support (gradients/standard errors) implementation, but backend abstraction in tranche-1 must preserve a dispatch point for derivative kernels.

## 17. Critique-Driven Guardrails (From Current Review)

This plan explicitly incorporates critical points raised in review:

1. Validate argmin displacement, not only objective value drift.
2. Do not hardcode 512; use bandwidth-aware dynamic grids with hard caps.
3. Treat conditional ratio consistency as a first-class requirement.
4. Resolve dependency strategy up front (backend abstraction, no external FFT dependency).
5. Require true LOO semantics for `cv.ml`.
6. Set realistic performance expectations for semiparametric methods.

## 18. Immediate Next Step (Implementation Kickoff)

Decision-log annex now exists:

1. `FFT_ACCELERATION_DECISION_LOG_2026-03-02.md`

Next task (first implementation sprint):

1. Create isolated spike worktree/branch:
- `git worktree add /Users/jracine/Development/np-master-fft-spike -b codex/spike-fft-tranche1 master`
2. Implement tranche-1 eligibility preflight and diagnostics in `npksum` call path (no FFT compute yet):
- gate checks,
- deterministic hard-stop reasons,
- grid/memory preflight reporting.

Blocking decisions from Section 16.1 are now mapped in the companion annex and currently marked `PENDING LOCK v3`; implementation starts only after explicit lock sign-off.

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
