# npcdens Fast-Path TODO

Date: 2026-02-19
Repo: /Users/jracine/Development/np-master
Primary file: /Users/jracine/Development/np-master/src/jksum.c

## Completed in this pass
- Added fast second-term collapse for `gate_x_all_large_fixed && BANDWIDTH_den==BW_FIXED` using `ky_jk` aggregates (`total/row/col/diag`) to avoid `kx_ij/kx_ik` work in this regime.
- Added fast joint-term (`jmean`) path that computes Y-only leave-one-out term when X is all-large-fixed, instead of full XY kernel build.
- Verified numerical identity for fixed-point fast-path benchmark (exact match for `fval`, `xbw`, `ybw`, and eval counters; elapsed differs only).

## Benchmark artifacts
- `/tmp/npcdens_point_prepost_summary_20260219_2.csv`
- `/tmp/npcdens_point_prepost_parity_20260219_2.csv`
- `/tmp/npcdens_prepost_summary_20260219_2.csv`
- `/tmp/npcdens_prepost_parity_20260219_2.csv`
- `/tmp/npcdens_point_prepost_summary_20260219_5.csv`
- `/tmp/npcdens_point_prepost_parity_20260219_5.csv`
- `/tmp/npcdens_prepost_summary_20260219_5.csv`
- `/tmp/npcdens_prepost_parity_20260219_5.csv`

## Current observations
- Fast-focused fixed-point run showed modest mean/median speedup after joint-term reduction.
- Full optimizer run also improved modestly, but those runs had `num.feval.fast = 0` and therefore do not directly isolate fast-branch benefit.

## Additional implementation in latest pass
- Removed heap allocation for fast joint operator (`y_joint_operator`) by using a stack buffer.
- Hoisted repeated `np_gate_ctx_set/clear` for X and Y out of inner loops in the second-term block.
- Removed dead `if(BANDWIDTH_den != BW_FIXED)` branches inside the fixed-bandwidth fast branch.

## Latest measured deltas
- Fast-focused fixed-point benchmark (10 paired reps): mean elapsed change improved to about `+0.37%` faster post, with exact parity on `fval`, `xbw`, `ybw`, and eval counters.
- Full optimizer benchmark (`npcdens(y~x)`, `n=1000`, 3 fixed + 3 varying):
  - fixed mean about `+5.71%`, median about `+12.44%`,
  - varying mean about `+2.78%`, median about `+2.05%`,
  - parity remained tight on objective (`~3e-12`) with expected optimizer-path differences in selected bandwidth/eval counts.

## Next painstaking opportunities
1. Workspace reuse
- Reuse `row_sum`, `col_sum`, `diag_sum`, and fast joint operator buffers across objective evaluations to reduce allocator overhead.

2. Fast-joint gating
- Evaluate adding a dedicated Y-normal fast gate for the fast joint (`jmean`) path to avoid unnecessary Y kernel work when Y bandwidth is effectively all-large.
- Keep parity checks strict because this path affects the first objective term directly.

3. Symmetry + block reductions (guarded)
- Consider symmetry reduction for `ky_jk` only under conditions that guarantee exact symmetry (bounded kernels may violate assumptions).
- Add explicit guards and fallback to existing path when conditions are not met.

4. Measurement discipline
- Continue fixed-point fast-hit timing (direct branch signal) plus fixed-seed and varying-seed optimizer timing.
- Always report mean/median percent change and parity (`max abs diff` for objective and selected bandwidths).

## Rejected variant in this pass
- Tree-enabled indexing for fast joint Y-only path (`int_TREE_Y`/`kdt_extern_Y`) did not improve this workload and was reverted.
