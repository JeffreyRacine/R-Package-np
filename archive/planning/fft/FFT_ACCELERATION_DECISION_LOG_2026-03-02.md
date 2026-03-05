# FFT Acceleration Decision Log Annex (Blocking Items)

Date: March 2, 2026  
Companion to: `FFT_ACCELERATION_PLAN_2026-03-02.md` (Section 16.1)

Purpose:

- Lock concrete tranche-1 policies for all blocking pre-code decisions.
- State direct testing consequences and user-visible behavior.

Decision status:

- `PENDING LOCK v3` (post-critique cycle 3 incorporation, March 2, 2026).
- Lock after final pre-code review confirms Decisions 1-25 below.

Interim supersession notice:

1. `FFT_ALIGNMENT_PLAN_V4_2026-03-02.md` Section 10 currently supersedes conflicting v3 items (notably sparse-profile selector handling, `cv.cdf` tranche status, selector reference-grid rule, and early vertical-slice gates).
2. Consolidate into a single lock document before implementation lock transition.
3. Execution-readiness authority now resides in `FFT_EXECUTION_READINESS_V5_2026-03-02.md`.

## Decision 1: Approximation Acceptance Contract

Chosen policy:

1. Primary acceptance metric is bandwidth selector stability (argmin), not objective value closeness.
2. Pass thresholds for selected bandwidths:
- median relative displacement <= 5%
- 90th percentile relative displacement <= 10%
3. For vector bandwidths, require both componentwise thresholds and vector relative norm threshold <= 10%.
4. Add downstream estimator parity gate:
- evaluate estimator at `h_fft` and `h_exact` on matched evaluation grids,
- require `L_inf` relative fitted-value error <= 2% in central support and <= 5% globally.
5. No systematic directional bias check is run per simulation setting (not pooled across settings).

Rationale:

1. Small objective drift can still produce materially different selected bandwidths.
2. Users care about fitted objects and inference, not CV scalar closeness in isolation.

Test implication:

1. Every selector benchmark stores both objective curves and chosen argmins.
2. Every selector benchmark stores downstream fitted-value parity for `h_fft` vs `h_exact`.
3. Bias sign tests are executed separately by DGP/method/kernel/dimension setting.

User-visible implication:

1. FFT path is documented as approximate but selector-stable and estimator-stable within explicit gates.

## Decision 2: FFT Backend Choice

Chosen policy:

1. No external FFTW dependency in tranche-1.
2. Implement backend abstraction.
3. Tranche-1 baseline backend uses `stats::fft` / `stats::mvfft` (mixed-radix, existing R infrastructure).
4. Custom internal C radix-2 FFT is deferred to tranche-2 unless tranche-1 performance gates fail.

Rationale:

1. Lowest maintenance and CRAN friction for tranche-1.
2. Provides immediate correctness path while preserving ability to optimize later.

Test implication:

1. Profile objective-evaluation overhead with `stats::fft` backend in small/medium regimes.
2. If tranche-1 speed gates are missed in representative medium cases, escalate custom C FFT before release lock.

User-visible implication:

1. No new external system dependency to enable `fft=TRUE`.

## Decision 3: `cv.ml` Semantics (True LOO + Kernel-Zero Table)

Chosen policy:

1. Implement true leave-one-out correction, not pseudo-LOO.
2. For each evaluation point in an FFT-eligible slice:
- `f_minus_i = (n_eff/(n_eff-1))*f_i - K0/((n_eff-1)*prod(h))`
3. Use `n_eff` per active slice/cell (not global `n` under slicing).
4. Apply floor after correction:
- `f_minus_i <- max(f_minus_i, DBL_MIN)` before log.
5. Maintain explicit tranche-1 `K0` table (order-2 kernels only):
- Gaussian (order 2): `K0 = 1/sqrt(2*pi)`
- Epanechnikov (order 2): `K0 = 3/4`
- Uniform: `K0 = 1/2`
6. Multivariate FFT blocks use product-kernel origin constants:
- `K0_d = prod_{j=1}^d K0_j`
- Example: 2D Gaussian product kernel uses `K0_2D = 1/(2*pi)`.
7. If any active slice/profile has `n_eff<=1`, `fft=TRUE` with `cv.ml` hard-stops with deterministic gate reason.
8. Higher-order kernels remain ineligible until analytic `K0` tables and kernel-sign LOO stress tests are added and passed.

Rationale:

1. Preserves current LOO objective semantics.
2. Makes correction explicit and auditable by kernel.

Test implication:

1. Compare corrected FFT-LOO `cv.ml` vs exact pairwise `cv.ml` for each supported tranche-1 kernel.
2. Add targeted stress tests for correction behavior near floor values.
3. Add explicit 2D Gaussian `K0_2D=1/(2*pi)` correction test.
4. Add sparse-slice tests where some profiles have `n_eff=1` and verify deterministic hard-stop.

User-visible implication:

1. `cv.ml` under FFT remains LOO ML by construction, with kernel-specific correction metadata available for audit.

## Decision 4: Tranche-1 Supported Kernel List

Chosen policy:

1. Tranche-1 continuous kernels:
- `gaussian` order 2
- `epanechnikov` order 2
- `uniform`
2. Tranche-1 FFT unsupported:
- `truncated gaussian`
- all higher-order continuous kernels (orders 4/6/8)
- any bounded-kernel mode (see Decision 5)

Rationale:

1. Lower-risk initial scope.
2. Avoids tranche-1 complexity from higher-order signed kernels.
3. Uniform is retained only with explicit ringing/aliasing validation gates.

Test implication:

1. Kernel-specific parity/performance matrix for exactly these three kernels.
2. Eligibility tests hard-stop unsupported kernels when `fft=TRUE`.
3. Uniform-specific anti-ringing parity tests are mandatory; failure defers uniform to tranche-2.

User-visible implication:

1. Clear supported-kernel contract for tranche-1.

## Decision 5: Boundary Policy

Chosen policy:

1. Tranche-1 FFT requires unbounded continuous kernel mode only:
- `ckerbound="none"`, `cxkerbound="none"`, `cykerbound="none"`.
2. Any finite bound mode with `fft=TRUE` hard-stops.

Rationale:

1. Correct bounded normalization on grids is nontrivial and high-risk for tranche-1.

Test implication:

1. Negative tests for bounded modes under `fft=TRUE`.
2. Deterministic stop-message tests.

User-visible implication:

1. Bounded-kernel users continue via standard path (`fft=FALSE`).

## Decision 6: Dynamic Grid Policy

Chosen policy:

1. No hardcoded fixed grid.
2. Per axis:
- `n_base_d = ceil((range_d + 2*pad_d)/delta_d)`
- `delta_d = h_d / m_k`
3. Tranche-1 kernel table (`m_k`, `s_k`):
- Gaussian (order 2): `m_k = 6`, `s_k = 4`
- Epanechnikov (order 2): `m_k = 4`, `s_k = 1`
- Uniform: `m_k = 4`, `s_k = 1`
4. Padding rule: `pad_d = s_k * h_d`.
5. FFT size per axis:
- `n_fft_d = 2^ceiling(log2(2 * max(n_base_d, grid.min.per.axis)))`
- bounded above by `grid.max.per.axis`.
6. Default `grid.min.per.axis = 64`; values below `64` require diagnostic warning.
7. Higher-order kernels cannot be enabled until kernel/order-specific `m_k` and `s_k` are re-derived and validated.

Rationale:

1. Links discretization to bandwidth and support.
2. Enforces linear-convolution safety via explicit zero-padding.

Test implication:

1. Grid-selection tests over range/bandwidth sweeps.
2. Aliasing regression tests with known convolution identities.
3. Uniform-kernel ringing tests across multiple `m_k` candidates; retain tranche-1 `m_k=4` only if parity gates pass.

User-visible implication:

1. Diagnostics report grid shape, spacing, and padding constants.

## Decision 7: Memory Policy

Chosen policy:

1. Default `memory.max.bytes = 268435456` (256 MB).
2. User override allowed via `fft.control`.
3. Workspace estimate:
- `bytes_est = c_work * prod(n_fft_d)`
- `c_work = 64` bytes/cell (tranche-1 planning estimate).
4. Hard-stop when `bytes_est > memory.max.bytes` and `fft=TRUE`.

Rationale:

1. More realistic default for interactive/shared environments.
2. Still allows power users to raise cap explicitly.

Test implication:

1. Threshold-boundary tests for preflight decisions.
2. Reproducibility tests for repeated calls and identical gate decisions.

User-visible implication:

1. Stop message reports estimated bytes, configured cap, and override hint.

## Decision 8: Interpolation Policy

Chosen policy:

1. Tranche-1 interpolation is fixed:
- linear (1D)
- bilinear (2D)
2. No hidden adaptive interpolation switching.

Rationale:

1. Predictable bias/smoothness behavior for optimization.

Test implication:

1. Interpolation consistency tests against analytic functions.
2. CV objective local smoothness checks around minima.

User-visible implication:

1. Interpolation mode is surfaced in diagnostics.

## Decision 9: Eligibility Error Contract

Chosen policy:

1. Hard-stop message schema is fixed for `fft=TRUE` ineligible calls:
- `gate`
- `observed`
- `required`
- `action`
2. First-fail deterministic ordering is fixed:
- Gate 1: `bwtype`
- Gate 2: kernel support
- Gate 3: boundary mode
- Gate 4: continuous FFT-block dimension
- Gate 5: output compatibility
- Gate 6: categorical expansion cap
- Gate 7: dynamic grid/memory feasibility

Rationale:

1. Improves debugging and reproducibility.

Test implication:

1. Golden-message tests per gate.
2. Deterministic ordering tests across platforms.
3. Multi-fail calls are tested to confirm exactly one primary first-fail gate reason.

User-visible implication:

1. Users can immediately identify required changes to become eligible.

## Decision 10: Conditional Ratio Guard Policy

Chosen policy:

1. For conditional ratio paths (`joint/marginal`): evaluate on grid-compatible support.
2. Denominator floor is fixed:
- `den = max(den, DBL_MIN)`
3. Emit diagnostic counter for floor hits.

Rationale:

1. Avoids ad hoc data-dependent floors and inconsistent bias.
2. Matches conservative existing numerical-guard style.

Test implication:

1. Stress tests in sparse tails and near-zero marginal regimes.
2. Bias check in central support region against exact route.

User-visible implication:

1. Stability in tails with transparent diagnostics when numerical guard activates.

## Decision 11: Output Compatibility Policy

Chosen policy:

1. Tranche-1 FFT enabled only for supported output contracts.
2. Under `fft=TRUE`, hard-stop for unsupported extras in tranche-1:
- regression/conditional gradients and standard errors where non-FFT moment pipelines are required.
- `npksum(return.kernel.weights=TRUE)`.
3. Base fitted-value paths are eligible.

Rationale:

1. Avoids partial correctness claims on advanced outputs.

Test implication:

1. Eligibility matrix tests by method and option combination.
2. Negative tests for unsupported-option hard-stops.

User-visible implication:

1. Clear tranche-1 support boundary for advanced outputs.

## Decision 12: Validation Battery and Artifact Contract

Chosen policy:

1. Required simulation families:
- Gaussian unimodal
- Gaussian mixtures (multimodal)
- heavy-tail (`t`)
- skewed continuous (lognormal/gamma)
2. Required sizes are locked:
- Small: `n in {50, 100}`
- Medium: `n in {250, 500}`
- Large: `n in {1000, 2500}`
- Crossover documentation: include at least one very-large run (`n in {5000, 10000}` where feasible).
3. Required dimensions: 1D and 2D continuous FFT block cases.
4. Required methods:
- `npudens*`, `npudist*`, `npcdens*` (`1y+1x`), `npcdist*` (`1y+1x`), `npqreg*` (`1y+1x`), `npreg*`, `npksum`.
5. Required artifacts under dated `/tmp` path:
- raw outputs
- selected bandwidths
- objective traces
- timing summaries
- parity summary table
- gate pass/fail report

Rationale:

1. Ensures selector and downstream estimator claims are evidence-backed.

Test implication:

1. No tranche pass without complete artifact set and threshold pass.

User-visible implication:

1. Release notes can state verified scope and known ineligible cases precisely.

## Decision 13: `cv.cdf` Semantics and Gates

Chosen policy:

1. `cv.cdf` is explicitly included in tranche-1 decision scope (not inferred from `cv.ls`/`cv.ml`).
2. Preserve existing objective semantics and explicit leave-one-out behavior for distribution bandwidth selection.
3. LOO correction uses cumulative-kernel self-contribution subtraction:
- `F_minus_i = (n_eff/(n_eff-1))*F_i - G0/(n_eff-1)` in 1D.
- For product-kernel FFT blocks, use `G0_d = prod_{j=1}^d G0_j`.
- For symmetric tranche-1 kernels, `G0_j = 1/2` on each continuous axis.
4. If any active slice/profile has `n_eff<=1`, `fft=TRUE` with `cv.cdf` hard-stops with deterministic gate reason.
5. Apply the same selector-stability and downstream parity thresholds as Decision 1.
6. If `cv.cdf` parity/selector gates fail in tranche-1 benchmarks, defer `npudistbw` FFT selector path to tranche-2 while keeping estimation FFT path eligible.

Rationale:

1. Distribution CV objective has distinct structure and must be validated directly.

Test implication:

1. Dedicated `cv.cdf` benchmark panel and artifact bundle are mandatory.
2. Add explicit LOO-correction formula parity tests against exact pairwise route for each tranche-1 kernel.
3. Add sparse-slice `n_eff<=1` deterministic hard-stop tests.

User-visible implication:

1. Clear pass/fail status for FFT-based distribution bandwidth selection in tranche-1 release notes.

## Decision 14: Bandwidth Unit Resolution (`bwscaling`)

Chosen policy:

1. FFT grid construction always uses bandwidths in data units.
2. Any stored/scaled bandwidth representation must be resolved to data units before computing `delta_d` and `pad_d`.
3. Unit-resolution helpers are centralized and reused by all FFT entry points (`npksum`, selectors, estimators).

Rationale:

1. Prevents silent over/under-smoothing from scale mismatch.

Test implication:

1. Add scale-factor invariance test (`sd=1` vs `sd=5`) that verifies expected 5x grid-spacing scaling in data units.

User-visible implication:

1. `fft.control` diagnostics report resolved data-unit bandwidths.

## Decision 15: CV Optimizer Architecture

Chosen policy:

1. Tranche-1 selectors use pre-binned workspace architecture by default.
2. Bin continuous data once on a reference grid before optimizer iterations.
3. Per bandwidth proposal, recompute kernel transform/convolution and objective using cached binned data.
4. Per-call re-binning is reserved for debug parity mode only.

Rationale:

1. Preserves FFT speed advantage inside optimization loops.

Test implication:

1. Benchmark objective-evaluation throughput with and without pre-binning and retain artifacts.

User-visible implication:

1. CV speedups are expected in low-dimensional fixed-bw settings, not only in one-off evaluations.

## Decision 16: Backend Call-Stack Contract (`stats::fft` Baseline)

Chosen policy:

1. Binning, interpolation, LOO corrections, and gate checks remain in C (`.Call` path).
2. FFT transform calls use `stats::fft`/`stats::mvfft` through a narrow R-level wrapper.
3. Interface contract is explicit: C builds vectors -> R wrapper transforms -> C resumes correction/objective work.
4. Keep backend dispatch pluggable for future internal C FFT.

Rationale:

1. Minimizes tranche-1 refactor risk while preserving swap-in path for faster backend later.

Test implication:

1. Add microbench for boundary-crossing overhead in objective loops.

User-visible implication:

1. No new external dependencies; backend can evolve without API changes.

## Decision 17: `cv.ls` Integral-Term Computation

Chosen policy:

1. In tranche-1, both `cv.ls` components are computed on the same FFT/grid approximation route.
2. No mixed analytic + grid term computation in tranche-1.
3. Closed-form analytic alternatives are deferred to tranche-2.

Rationale:

1. Avoids structural asymmetry in approximation bias across objective terms.

Test implication:

1. Objective-surface parity tests around minima are mandatory for each supported kernel.

User-visible implication:

1. `cv.ls` approximation behavior is consistent and documented.

## Decision 18: `npqreg` CDF Inversion Contract

Chosen policy:

1. Use bracketed bisection inversion in tranche-1.
2. Enforce monotonicity on grid CDF via running-max before inversion.
3. Lock tolerance and iteration caps in implementation constants and expose in diagnostics.
4. If inversion fails to converge, hard-stop under `fft=TRUE` with deterministic reason.

Rationale:

1. Prevents silent quantile mis-location from non-monotone numerical CDF artifacts.

Test implication:

1. Add heavy-tail and multimodal conditional tests with inversion convergence checks.

User-visible implication:

1. Quantile failures are explicit and reproducible under `fft=TRUE`.

## Decision 19: `npudist` CDF Recovery from FFT Density

Chosen policy:

1. Recover CDF by cumulative integration over the unpadded data domain only.
2. Exclude zero-padding cells from CDF accumulation.
3. Normalize CDF at right data boundary; if normalization deviates beyond tolerance, hard-stop under `fft=TRUE`.

Rationale:

1. Prevents tail distortion from padded-domain accumulation errors.

Test implication:

1. Add analytic-reference tests for CDF normalization and boundary behavior.

User-visible implication:

1. Distribution outputs include normalization diagnostics when `verbose=TRUE`.

## Decision 20: `strict` Control Surface

Chosen policy:

1. Remove `fft.control$strict` from tranche-1 API surface.
2. Hard-stop semantics are fixed whenever `fft=TRUE` and any eligibility/contract check fails.

Rationale:

1. Avoids contradictory soft-fail semantics in tranche-1.

Test implication:

1. Ensure no code path silently degrades from `fft=TRUE` to non-FFT.

User-visible implication:

1. `fft=TRUE` behavior is deterministic and explicit.

## Decision 21: NA Handling Contract

Chosen policy:

1. FFT path requires NA-free continuous FFT-block inputs after existing `na.action` handling.
2. Preflight validates NA-free requirement before binning.
3. Any NA detected at FFT preflight under `fft=TRUE` hard-stops with deterministic reason.

Rationale:

1. Prevents NA poisoning or undefined C behavior in grid construction.

Test implication:

1. Add NA-focused negative tests to confirm stop behavior matches existing `na.action` contract.

User-visible implication:

1. Failures identify the NA gate and required remediation.

## Decision 22: C Workspace State Policy

Chosen policy:

1. No global/static mutable C state in FFT backend.
2. Workspace is per-call allocated and freed via R-compatible memory management.
3. Cross-call correctness must not depend on prior call sizes or settings.

Rationale:

1. Prevents cross-call contamination and hidden state bugs.

Test implication:

1. Add repeated-call tests with varying `n`/`h` orderings and verify reproducible outputs.

User-visible implication:

1. Deterministic behavior across sessions and call sequences.

## Decision 23: FFT Metadata Placement and Inheritance

Chosen policy:

1. FFT metadata (`fft_*`) is attached to both bandwidth objects and estimator objects.
2. Estimator objects inherit bw metadata when a bw object is supplied, then append run-time realization fields.
3. Merge policy is deterministic: estimator run-time values take precedence over inherited defaults.

Rationale:

1. Preserves provenance across bandwidth-selection and estimation stages.

Test implication:

1. Add object-structure tests for metadata presence and deterministic merge behavior.

User-visible implication:

1. Users can audit FFT provenance on both bw and fit objects.

## Decision 24: Development Branch/Worktree Isolation

Chosen policy:

1. FFT implementation starts in throwaway worktree branch:
- branch: `codex/spike-fft-tranche1`
- worktree: `/Users/jracine/Development/np-master-fft-spike`
2. Local `master` remains available for unrelated minor work.
3. Spike may be deleted entirely with no merge if feasibility fails.

Rationale:

1. De-risks experimentation without blocking routine maintenance.

Test implication:

1. All tranche artifacts and commits reference spike branch provenance.

User-visible implication:

1. FFT experiment can be abandoned cleanly with no `master` contamination.

## Decision 25: Lock Transition and Reopen Conditions

Chosen policy:

1. Transition from `PENDING LOCK v3` to `LOCKED` requires explicit sign-off after reading this annex revision.
2. Reopen tranche-1 policy decisions only if one of the following occurs:
- validation battery fails any hard gate,
- `stats::fft` backend misses tranche-1 speed gates in representative medium workloads,
- a user-facing correctness issue is observed in supported scope.

Rationale:

1. Preserves rigor while allowing controlled iteration if evidence contradicts assumptions.

Test implication:

1. Lock status change is recorded as a dated doc revision with no implicit state change.

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
