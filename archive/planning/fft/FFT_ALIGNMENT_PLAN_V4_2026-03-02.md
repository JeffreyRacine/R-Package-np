# FFT Alignment Plan v4 (Two Introspective Rounds)

Date: March 2, 2026
Scope: `np-master` planning only (no implementation in this document)
Status: pre-code plan, iteration v4

Authority note:

1. This v4 memo is the active planning authority for FFT tranche sequencing.
2. A supersession table in Section 10 explicitly resolves conflicts with prior planning docs.
3. Superseded by `FFT_EXECUTION_READINESS_V5_2026-03-02.md` for execution-readiness lock decisions.

## 0. Intent

This memo is a self-critical refinement pass over prior FFT planning artifacts.

Goal:

1. Produce a plan I fully endorse technically.
2. Surface where I changed position and why.
3. Avoid architectural dead ends (including future derivative support).

Development isolation policy (retained):

1. FFT work starts in throwaway worktree branch:
- `git worktree add /Users/jracine/Development/np-master-fft-spike -b codex/spike-fft-tranche1 master`
2. `master` stays free for unrelated maintenance.

## 1. Ground Truth Re-Audit (Code-Level)

This pass was tied to current `np-master` internals, not abstract formulas.

Key observations:

1. CV objectives are implemented by dedicated C routines in `src/kernelcv.c` and `src/jksum.c`:
- density ML/LS wrappers: `np_cv_func_density_categorical_ml`, `np_cv_func_density_categorical_ls`
- conditional density ML/LS wrappers: `np_cv_func_con_density_categorical_ml`, `np_cv_func_con_density_categorical_ls_npksum`
- distribution LS wrappers: `cv_func_distribution_categorical_ls`, `cv_func_con_distribution_categorical_ls`
2. `cv.cdf` is not a trivial density-style constant subtraction; current code computes leave-one-out via explicit weighted sums with indicator terms.
3. `bwscaling` is real and active in estimator/CV paths; conversion between scale factors and raw bandwidths is explicit in R bw objects (`bw`, `bandwidth`, `sfactor`) and in C setup (`int_LARGE_SF`, `kernel_bandwidth_mean`).
4. `npksum` already has operator semantics (`OP_NORMAL`, `OP_DERIVATIVE`, `OP_INTEGRAL`) and kernel-index offsets (`OP_CFUN_OFFSETS`) in C; this is the right abstraction anchor for derivative/integral future support.
5. Existing code already has large-bandwidth shortcuts and kernel-at-zero machinery (`np_cont_largeh_k0`), including higher-order kernels and rectangular kernel, which confirms kernel-zero constants are kernel/order specific.
6. `kbandwidth` explicitly rejects `bwscaling=TRUE` for `npksum`, while bandwidth selector objects in other paths carry both scaled and raw representations.

## 2. Position-by-Position Stance Review

Each item is explicitly accepted, modified, or rejected.

1. Hard-stop semantics for `fft=TRUE`: ACCEPT.
- Reason: reproducibility and benchmarking integrity.
- Policy: no silent global fallback.

2. One/two-dimensional fixed-bandwidth tranche entry: ACCEPT.
- Reason: most realistic speedup zone with manageable risk.
- Policy: `q in {1,2}` and `bwtype == "fixed"` for FFT-smoothed block.

3. Uniform in tranche-1: MODIFY.
- Reason: keep only with mandatory anti-ringing gate.
- Policy: uniform is conditionally enabled by test outcome, not presumed safe.

4. Higher-order kernels in tranche-1: REJECT.
- Reason: high risk for LOO and grid-support calibration.
- Policy: defer orders 4/6/8.

5. `bwscaling` handling as a risk note only: REJECT.
- Reason: this is a silent-wrong-answer risk.
- Policy: hard decision with explicit unit-resolution contract.

6. Single-constant `K0` logic for all CV: REJECT.
- Reason: valid only for specific objectives/contexts; conditional and distribution CV are not equivalent simplifications.
- Policy: objective-specific LOO formulations mirrored from current C semantics.

7. `n_eff<=1` unconditional hard stop in all CV contexts: MODIFY.
- Reason: too blunt for optimizer flow and profile-sliced contexts.
- Policy: explicit invalid-objective candidate handling inside selectors; hard-stop only when call-level preconditions are structurally violated.

8. Pre-binned CV architecture with bandwidth-dependent grid formula unchanged: REJECT as internally inconsistent.
- Reason: fixed cached bins and `delta = h/m` cannot both be primary policy.
- Policy: reference-grid policy tied to optimizer search bounds with deterministic regrid rules.

9. `stats::fft` backend as tranche-1 default: ACCEPT WITH PERFORMANCE GATE.
- Reason: best portability/maintenance starting point.
- Policy: keep only if overhead gates pass; else escalate internal C FFT sooner.

10. `cv.ls` mixed analytic/grid terms: REJECT.
- Reason: asymmetric approximation bias can move argmins.
- Policy: same approximation family for both objective terms in tranche-1.

11. `cv.cdf` included by analogy: REJECT.
- Reason: objective structure differs materially.
- Policy: separate decision track and potential deferral.

12. `npqreg` inversion loosely specified: REJECT.
- Reason: non-monotone CDF artifacts can silently break quantiles.
- Policy: deterministic inversion + monotonicity enforcement + fail contract.

13. Keep `strict` toggle in tranche-1: REJECT.
- Reason: conflicts with hard-stop principle and multiplies edge behavior.
- Policy: remove `strict` for tranche-1.

14. Metadata only on one object class: REJECT.
- Reason: provenance loss across bw/fit workflow.
- Policy: attach to both bw and fit objects with deterministic merge.

15. Ignore derivative future now: REJECT.
- Reason: high retrofit risk.
- Policy: enforce operator-dispatch-ready architecture now, even if derivative kernels are tranche-2.

## 3. Round 1 Update (Ambitious Draft)

Round-1 draft (intentionally broad):

1. Tranche-1 includes estimation and selector FFT for:
- `npudens*`, `npudist*`, `npcdens* (1y+1x)`, `npcdist* (1y+1x)`, `npqreg* (1y+1x)`, `npreg*`, `npksum`.
2. CV acceleration includes `cv.ls`, `cv.ml`, and `cv.cdf` for eligible cases.
3. Backend baseline is `stats::fft`/`stats::mvfft` through a C<->R wrapper.
4. Pre-binned objective acceleration used in selector loops.
5. Kernel set: Gaussian2/Epanechnikov2/Uniform (conditional on anti-ringing tests).
6. No bounded kernels in tranche-1.
7. Hard-stop global eligibility failures.

## 4. Round 1 Critique (Why Round 1 Is Not Yet Safe)

Critical issues discovered in self-review:

1. `cv.cdf` formula lock is premature.
- Existing implementation uses structured leave-one-out weighted sums with indicator terms, not a simple constant correction.

2. Conditional ML/LS selectors are materially more complex than unconditional selectors.
- Current C code computes leave-one-out ratios from joint and marginal weighted sums with profile/kernel interactions.

3. Pre-binning and bandwidth-tied grid policy conflict.
- A selector loop cannot be both fixed-grid cached and strictly `delta = h/m` without regrid rules.

4. `n_eff` policy was over-simplified.
- Per-profile degeneracy should not automatically kill whole optimization unless structurally unavoidable.

5. `stats::fft` overhead risk is unresolved for tight objective loops.
- Architecture must include a measurable boundary-crossing budget.

Conclusion of Round-1 critique:

- Round-1 scope is too wide for a first coding sprint if the objective is low-regret execution.

## 5. Round 2 Update (Plan I Endorse)

Round-2 narrows and sequences work to preserve momentum while keeping the long-run architecture intact.

### 5.1 Tranche Structure

Tranche 1A (implementation start target):

1. FFT infrastructure in `npksum` for 1D/2D continuous fixed-bandwidth blocks (`OP_NORMAL` only).
2. Estimation-path acceleration (no selector acceleration yet) for eligible calls:
- `npudens`, `npudist`, `npreg`, `npcdens (1y+1x)`, `npcdist (1y+1x)`, `npqreg (1y+1x)`.
3. Kernel set: Gaussian2, Epanechnikov2.
4. Uniform enabled only after anti-ringing gate passes.
5. Hard-stop on global ineligibility under `fft=TRUE`.

Tranche 1B (after 1A parity/speed pass):

1. Selector acceleration for unconditional density:
- `npudensbw` with `cv.ls`, `cv.ml`.
2. `npudistbw` (`cv.cdf`) only if objective-equivalence derivation + benchmark gate passes.
3. Conditional selectors (`npcdensbw`, `npcdistbw`) deferred unless dedicated equivalence proof succeeds.

Tranche 2:

1. Higher-order kernels.
2. Bounded-kernel FFT support.
3. Conditional selector FFT.
4. Derivative/stderr FFT compute paths.
5. MPI adaptation (`np-npRmpi`).

### 5.2 CV Semantics Contract (Revised)

1. Never replace existing C objective definitions with simplified formulas unless analytically proven equivalent for that objective.
2. For each selector we mirror existing algebra from current C routines:
- `cv.ml` unconditional: leave-one-out density likelihood objective as currently computed.
- `cv.ml` conditional: ratio-of-leave-one-out weighted sums objective as currently computed.
- `cv.ls` density: preserve convolution + LOO structure.
- `cv.cdf`: preserve indicator-vs-LOO-CDF LS structure.
3. Any FFT-specific closed form (e.g., constant self-contribution) is objective-specific and must be validated against exact pairwise baseline.
4. Selector infeasible-candidate policy for sparse profiles:
- if any profile/slice yields `n_eff<=1` for a bandwidth proposal, return `+Inf` for that objective evaluation;
- do not silently drop that profile from the objective sum;
- record diagnostic counters for infeasible-candidate hits.

### 5.3 Grid/Caching Contract (Resolved)

1. Selector FFT uses a reference grid cache per optimization run.
2. Reference grid lower bound policy when user does not supply explicit optimizer bounds:
- compute initialization bandwidth `h_init` using the existing selector initialization path (normal-reference or method-specific initializer);
- set finest reference bandwidth to `h_ref_min = h_init/10`;
- derive `delta_ref` from `h_ref_min` and locked cells-per-bandwidth constants;
- compute once before optimizer entry.
3. Candidate `h` checks:
- if `h/delta` stays above minimum cells-per-bandwidth threshold, reuse cache.
- otherwise trigger deterministic regrid-and-rebin (counted and logged) or reject candidate per locked policy.
4. Preflight computes memory for worst-case active grid state under this policy.

### 5.4 `bwscaling` Unit Contract

1. FFT compute always operates on data-unit bandwidths.
2. For selector internals, bandwidth conversion follows the same path currently used by C (`kernel_bandwidth_mean` pathway).
3. For estimator calls using bw objects, resolved data-unit bandwidths are explicit in diagnostics (`fft_bandwidth_data_units`).

### 5.5 Backend Contract

1. Backend interface remains in C and is dispatch-based:
- `NP_FFT_BACKEND_NONE`
- `NP_FFT_BACKEND_STATSFFT`
- future `NP_FFT_BACKEND_C`.
2. Tranche-1 uses `stats::fft`/`stats::mvfft` backend.
3. Boundary-crossing overhead is benchmarked explicitly; fail gate triggers backend escalation.
4. Quantitative gate for selector loops:
- R<->C boundary-crossing overhead must be < 20% of total objective-evaluation time at `n=500` in the 1D Gaussian2 selector benchmark;
- otherwise escalate backend optimization before widening selector scope.

### 5.6 Operator Extensibility Contract (Derivative-Ready)

1. FFT core API accepts operator tags aligned with existing `npksum` semantics:
- `OP_NORMAL`, `OP_INTEGRAL`, `OP_DERIVATIVE`.
2. Tranche-1 implements `OP_NORMAL` only.
3. Tranche-1 must include dispatch slots/tests proving unsupported operators fail deterministically (not silently).
4. This prevents architectural dead ends for gradients/stderr in tranche-2.

### 5.7 Robustness Contracts

1. No global/static mutable FFT workspace state in C.
2. NA checks enforce FFT-block cleanliness after existing `na.action` behavior.
3. Metadata attached to both bw and fit objects.
4. Deterministic gate ordering retained.
5. Validation sizes locked (`n={50,100}`, `{250,500}`, `{1000,2500}`, and one very-large run).

## 6. Round 2 Critique (Residual Risks Still Worth Challenging)

Even after Round-2 tightening, three meaningful risks remain:

1. `stats::fft` call-overhead uncertainty in optimizer-heavy CV loops.
- Mitigation: explicit microbench gate before committing selector acceleration scope.

2. Uniform-kernel ringing risk.
- Mitigation: treat uniform as probationary; auto-defer if parity gate fails.

3. `cv.cdf` feasibility risk.
- Mitigation: keep `npudistbw` selector FFT as conditional deliverable in tranche-1B, not guaranteed tranche-1A.
4. `cv.cdf` equivalence proof risk.
- Mitigation: require a short written derivation artifact mapping FFT terms to existing C objective terms before enabling `npudistbw` selector FFT.

These are acceptable residuals because they are now explicit go/no-go gates, not hidden assumptions.

## 7. Must-Decide-Before-Code (v4 Lock List)

These are the only blockers I consider truly mandatory before coding:

1. Approve tranche sequencing (`1A estimation-first`, `1B selector-follow-on`).
2. Approve CV semantics rule: objective-equivalence-to-current-C is mandatory; no simplified replacement without proof.
3. Approve grid/caching rule (reference grid + deterministic regrid criteria).
4. Approve `bwscaling` data-unit contract and diagnostic fields.
5. Approve backend policy (`stats::fft` first, benchmark-gated escalation).
6. Approve operator-extensibility requirement for derivative future-proofing.
7. Approve conditional deferral policy for selector FFT (`cv.cdf` and conditional selectors contingent on proof/gates).

## 8. First Coding Sprint Plan (After Lock)

1. Create spike worktree branch.
2. Add FFT gate/preflight scaffolding in `npksum` with metadata output only (no compute).
3. Add backend dispatch skeleton and `OP_*` contract checks.
4. Implement 1D Gaussian2 estimation FFT for unconditional density as vertical slice.
5. Minimal vertical-slice exit gate (must all pass before any scope expansion):
- parity: `L_inf` relative error <= 2% on held-out Gaussian unimodal data for `n in {100,500}` (1D);
- speed: measured speedup >= 2x at `n=500` in the locked benchmark script;
- safety: zero regressions on existing `np-master` test suite.
6. Validate and archive parity/performance artifacts under dated `/tmp` paths.
7. Expand to 2D and additional tranche-1 kernels only after Step 5 gate passes.

This sequence maximizes early signal with minimal sunk cost.

## 9. Why I Am Fully Aligned With This Plan

1. It keeps the end-state ambition intact (broad estimator family, derivative-ready architecture).
2. It removes the mathematically least-certain commitments from the first coding sprint.
3. It encodes uncertainty as explicit gates and deferrals, not hidden assumptions.
4. It preserves your practical objective: meaningful low-dimensional speedups for common workflows.

## 10. Supersession Table (v4 vs Prior Planning Docs)

This table resolves policy conflicts with:

1. `FFT_ACCELERATION_DECISION_LOG_2026-03-02.md` (PENDING LOCK v3)
2. `FFT_ACCELERATION_PLAN_2026-03-02.md`

Superseded positions:

1. Sparse-profile selector behavior:
- prior v3 language: hard-stop on `n_eff<=1` in `cv.ml`/`cv.cdf`;
- v4 replacement: infeasible-candidate handling returns `+Inf` per bandwidth proposal with diagnostics.

2. `cv.cdf` tranche status:
- prior v3 language: tranche-1 first-class inclusion;
- v4 replacement: tranche-1B conditional deliverable requiring objective-equivalence derivation + benchmark gates.

3. Selector grid policy:
- prior v3 language: dynamic grid from candidate `h` without explicit selector reference-grid lower-bound rule;
- v4 replacement: reference-grid cached policy with `h_ref_min = h_init/10` when no explicit bounds are supplied.

4. Backend performance gate specificity:
- prior v3 language: qualitative overhead concern;
- v4 replacement: quantitative selector-loop overhead gate (<20% at `n=500`, 1D Gaussian2 benchmark).

5. Vertical-slice go/no-go criterion:
- prior docs: generic “validate artifacts” wording;
- v4 replacement: explicit three-part minimal exit gate (parity, speed, safety) before any scope expansion.

Precedence rule:

1. Where conflicts exist, this v4 memo governs until a consolidated decision log revision is published.
