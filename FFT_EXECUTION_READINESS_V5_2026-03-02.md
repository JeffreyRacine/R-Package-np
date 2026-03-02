# FFT Execution Readiness v5 (Critical-of-Critique Pass)

Date: March 2, 2026
Scope: `np-master` planning only
Status: near-execution planning contract

Consolidation status:

1. This is the single canonical FFT planning document in repo root.
2. Superseded drafts are archived in `archive/planning/fft/`.

Execution strategy lock:

1. Execute FFT work in `np-master` only through tranche-1 completion gates.
2. `np-npRmpi` work is explicitly deferred until `np-master` tranche gates pass and feasibility is confirmed.

## 1. Purpose

This document intentionally critiques the prior critique cycle rather than accepting it wholesale.

Objective:

1. Lock only the recommendations that are statistically and engineering-sound.
2. Modify or reject recommendations that are too brittle, underspecified, or machine-dependent.
3. Define milestone gates that determine go/no-go progression.

## 2. Critique of Prior Critique (Accept / Modify / Reject)

### 2.1 Document Reconciliation

Decision: ACCEPT.

Rationale:

1. Multiple planning artifacts with conflicting policies are implementation-risk multipliers.

Action:

1. v5 is the execution-readiness authority.
2. Prior docs are reference context only and archived under `archive/planning/fft/`.

### 2.2 `n_eff<=1` Candidate Policy in Selectors

Prior suggestion:

1. Return `+Inf` objective for any candidate with a profile/slice where `n_eff<=1`.

Decision: MODIFY.

Why modify:

1. Unconditional `+Inf` is safe but can over-constrain search in mixed-data cases if rare sparse profiles appear.
2. The existing C objective behavior differs by selector/method; we should match method semantics, not impose a universal new rule.

Locked policy:

1. For tranche-1 selector FFT, if `n_eff<=1` occurs in any required LOO term:
- default objective return is `DBL_MAX` (equivalent optimizer penalty),
- diagnostics record invalid-profile count and share.
2. If invalid-profile share exceeds `fft.control$invalid.profile.max.share` (default `0` in tranche-1), hard-stop the selector call under `fft=TRUE`.
3. No silent profile dropping.

### 2.3 Reference-Grid Lower Bound (`h_init/10`)

Prior suggestion:

1. Use `h_ref_min = h_init/10` when user bounds absent.

Decision: MODIFY.

Why modify:

1. Single factor `10` is arbitrary and can be too coarse/fine across kernels, dimension, and scale.
2. This can create avoidable memory pressure or under-resolution.

Locked policy:

1. If user lower bounds exist, use them.
2. Else compute
- `h_ref_min_d = max(h_floor_d, h_init_d / c_ref)`
- defaults: `c_ref = 8`, `h_floor_d = max(machine_floor, 0.01 * sd_d * n^(-1/5))`.
3. Preflight must verify resulting grid passes memory cap; if not, increase `h_ref_min_d` deterministically and report adjustment.

### 2.4 Vertical Slice Speed Gate (`>=2x at n=500`)

Prior suggestion:

1. Require speedup `>=2x` at `n=500`.

Decision: MODIFY.

Why modify:

1. Single-point speed gate is hardware- and implementation-path dependent.
2. It can reject good implementations that still deliver strong crossover performance.

Locked policy:

1. Vertical-slice performance gate requires:
- no regression worse than 15% at `n=100`,
- speedup >= 1.5x at one of `n in {500, 1000}`,
- documented crossover estimate `n*` where FFT overtakes exact route.
2. Keep parity and test-regression gates mandatory.

### 2.5 Backend Boundary-Crossing Overhead Gate (`<20% at n=500`)

Prior suggestion:

1. Require R<->C overhead <20% at `n=500`.

Decision: MODIFY.

Why modify:

1. Fixed 20% threshold at one `n` is too brittle for route changes and CPU variance.

Locked policy:

1. Measure overhead share at `n={250,500,1000}`.
2. Gate passes if:
- median overhead share <= 30%, and
- total objective speed gate from 2.4 passes.
3. If gate fails, escalate backend optimization before widening selector scope.

### 2.6 `cv.cdf` Written Derivation Requirement

Decision: ACCEPT.

Rationale:

1. Benchmark agreement is insufficient without term-level objective correspondence.

Action:

1. Require a derivation memo mapping FFT-computed terms to existing C objective terms before enabling `npudistbw` FFT selector path.

## 3. Final Direction I Endorse

### 3.1 Scope and Sequencing

Tranche 1A (execute first):

1. `npksum` FFT infrastructure for 1D/2D continuous fixed-bandwidth blocks (`OP_NORMAL` only).
2. Estimation-path FFT acceleration (no selector FFT yet):
- `npudens`, `npudist`, `npreg`, `npcdens (1y+1x)`, `npcdist (1y+1x)`, `npqreg (1y+1x)`.
3. Kernel set:
- Gaussian2, Epanechnikov2.
- Uniform only if anti-ringing gate passes.

Tranche 1B (only after 1A gates pass):

1. Selector FFT for `npudensbw` (`cv.ls`, `cv.ml`).
2. `npudistbw` (`cv.cdf`) only after derivation memo + benchmark gates.
3. Conditional selector FFT remains deferred unless objective-equivalence proof succeeds.

### 3.2 Non-Negotiables

1. `fft=TRUE` never silently falls back globally.
2. `bwtype == "fixed"` and continuous FFT block dimension in `{1,2}`.
3. No bounded-kernel FFT in tranche-1.
4. `bwscaling` resolved to data units before grid construction.
5. No global/static mutable FFT workspace state.
6. Metadata on both bw and fit objects.
7. Operator-ready backend interface includes `OP_NORMAL`, `OP_INTEGRAL`, `OP_DERIVATIVE` tags even if only `OP_NORMAL` implemented now.

## 4. Milestones and Success Gates

### M0: Governance Lock (Planning Closure)

Deliverables:

1. This v5 file approved.
2. Consolidated policy map in prior docs marked non-authoritative for conflicts.

Success gate:

1. Zero unresolved policy conflicts across active planning artifacts.
2. `np-master`-only strategy acknowledged and locked for tranche-1 execution.

### M1: Preflight + Dispatch Skeleton

Deliverables:

1. Eligibility checks and deterministic first-fail ordering in `npksum` FFT path.
2. Backend dispatch skeleton and operator contract checks.
3. Grid-sizing scaffold using `nextn`-based axis sizing under memory caps.
4. Diagnostic metadata scaffold.

Success gate:

1. `fft=FALSE` behavior unchanged.
2. `fft=TRUE` ineligible paths always deterministic and informative.
3. Existing tests pass.

### M2: Vertical Slice (1D Gaussian2, Unconditional Density Estimation)

Deliverables:

1. Working FFT estimation path for `npudens` in 1D Gaussian2.
2. Benchmark and parity harness artifacts.

Success gate:

1. Parity:
- `L_inf` relative error <= 2% central support, <= 5% global at `n={100,500}`.
2. Performance:
- no slowdown worse than 15% at `n=100`,
- speedup >= 1.5x at one of `n in {500,1000}`,
- crossover `n*` estimated and reported.
3. Safety:
- zero regression on test suite.

### M3: 2D + Kernel Expansion Within 1A

Deliverables:

1. 2D continuous block estimation support.
2. Epanechnikov2 support.
3. Uniform support candidate.

Success gate:

1. Same parity/performance gates as M2.
2. Uniform anti-ringing test pass; else uniform deferred automatically.

### M4: Selector Enablement (`npudensbw` only)

Deliverables:

1. FFT selector path for `cv.ls` and `cv.ml` unconditional density.
2. Objective-equivalence notes vs existing C routines.

Success gate:

1. Argmin displacement gates:
- median <= 5%, p90 <= 10%.
2. Downstream fit parity gate passes at selected bandwidths.
3. Overhead gate:
- median boundary-crossing share <= 30% across `n={250,500,1000}`.

### M5: Conditional `cv.cdf` Decision Gate

Deliverables:

1. Written derivation memo for `cv.cdf` term mapping.
2. Prototype benchmarks for `npudistbw` selector FFT path.

Success gate:

1. Derivation accepted.
2. Benchmarks satisfy selector and downstream parity thresholds.
3. If either fails, defer `npudistbw` selector FFT to tranche-2 (estimation path may remain enabled).

## 5. Measurement and Artifacts

Required artifact bundle per milestone:

1. Script and exact invocation.
2. Data size/dimension/kernel/bwtype settings.
3. Timing summaries (mean/median, crossover estimate).
4. Parity summaries and threshold results.
5. Gate pass/fail report.
6. Dated `/tmp` storage paths.

Timing reliability protocol under background load:

1. Assume concurrent machine load (including user simulations) can perturb wall-clock timing.
2. For each timing claim, run at least 5 repetitions and report median and MAD (or IQR).
3. Do not use a single-run timing to pass/fail a milestone gate.
4. If variability is high (MAD/median > 15%), increase repetitions and report uncertainty explicitly.
5. Prefer crossover-zone conclusions from trend across `n` rather than one-size point estimates.

## 5.1 Checkpoint and Rollback Protocol

1. Create one local git checkpoint commit at the end of each milestone (`M1`, `M2`, ...).
2. Commit message format:
- `checkpoint: fft Mx <short-description>`
3. Record checkpoint hash and artifact directory in milestone notes.
4. Rollback path:
- `git checkout <checkpoint-hash>` for exact-state restore, or
- `git revert <checkpoint-hash>` for preserving linear history.
5. Never use destructive reset for milestone rollback in this workflow.

## 6. Immediate Next Planning Step (Before Code)

1. Confirm v5 as lock candidate.
2. Reflect v5 supersession status into prior decision log and plan headers.
3. Then create worktree and start M1 only.

## 7. Pre-Execution Investigation: Base R FFT (Local Environment)

Environment checked:

1. R version: `4.5.2` (`aarch64-apple-darwin20`).
2. `stats::fft` and `stats::mvfft` are `.Call` wrappers (`C_fft`, `C_mvfft`).
3. `stats::nextn` uses smooth-factor targets by default (`factors = c(2,3,5)`).

Functional checks:

1. `fft()` output is complex; inverse reconstruction requires division by transform length.
2. `mvfft()` reconstruction behavior is consistent with column-wise FFT semantics.

Microbench observations (local, indicative):

1. Repeated FFT spectral-convolution proxy scales as expected with `n log n`.
2. Caching the data FFT (selector-style inner loop) reduced per-eval cost materially versus recomputing data FFT each candidate (roughly `1.37x` to `1.75x` faster in tested sizes).
3. For non-smooth lengths, `nextn` padding can improve speed substantially (example tested `n=1537`, `nextn=1600` gave ~2x per-eval speed improvement).

Planning implications:

1. Keep pre-binned/cached selector architecture.
2. Replace any power-of-two-only assumption with:
- `n_fft_d = nextn(n_base_d)` (default factors `2,3,5`) subject to memory cap.
3. Retain backend-overhead milestone gates; local findings are positive but not a substitute for milestone validation on target workloads.
