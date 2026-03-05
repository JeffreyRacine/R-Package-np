# Issue Notes Index (`np-npRmpi`)

Last refresh: 2026-03-01

## Active / Keep In Place
- `issue_notes/bounded_kernel_todo.md` (active TODO)
- `issue_notes/issue_3_local_linear_high_order.md` (open)
- `issue_notes/issue_44_npindex_bvcov_recycling.md` (open)
- `issue_notes/issue_runmode_check_completion_hang.md` (open)
- `issue_notes/NATIVE_REGISTRATION_POLICY_20260301.md` (active policy record)

## Archived (Resolved / Historical)
Archive directory:
- `issue_notes/archive/resolved_20260301/`

Archived notes:
- `forensic_audit_autodispatch_plot_2026-02-20.md`
- `issue_10_npcdens_gradients_abort.md`
- `issue_13_npudens_avar_categorical.md`
- `issue_16_npregbw_rule_of_thumb.md`
- `issue_18_npregiv_multidim_instruments.md`
- `issue_26_formula_eval_explodepipe.md`
- `issue_4_npcdensbw_cvls_cvml_segfault.md`
- `issue_4_npcdensbw_segfault_cvml.md`
- `issue_50_np_pairs_plot.md`
- `issue_51_npregiv_exogenous_w.md`
- `issue_5_npqreg_tau_validation.md`
- `issue_6_npindexbw_bwtype.md`
- `issue_7_uocquantile_factor_subset.md`

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
