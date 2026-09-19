# Native source map

This is a navigation guide, not a replacement for the implementation, public
method documentation, or the ownership contracts beside declarations. Symbol
names below are search anchors, not a promise that every option visits that
function. Search objectives, fitted values, derivatives and uncertainty have
distinct owners even when they share kernel rows or solves.

## Start from the boundary

1. Find the public method in `../R/`, then its registered `C_np_*` call.
2. Check the signature and arity in [np_init.c](np_init.c).
3. Follow the entry in [np.c](np.c) or its smaller implementation file below.
4. Follow preparation, accumulation, solve, result assembly and cleanup together.
   A helper's return value alone does not describe its caller's error policy.

Examples: `C_np_regression` reaches
`kernel_estimate_regression_categorical_tree_np`; the
`C_np_regression_lp_apply_conditional` family connects to
`np_regression_lp_apply_matrix`, `np_regression_lp_hat_matrix` and
`np_regression_lp_sigtest_iid`. `C_np_kernelsum` reaches the weighted-sum
machinery. These are separate consumer contracts, not interchangeable estimates.

## Responsibilities and stable anchors

| Area | Files and anchors | Boundary to retain |
| --- | --- | --- |
| R/native boundary, search preparation, progress | [np.c](np.c), [np_init.c](np_init.c), [headers.h](headers.h); `C_np_regression`, `C_np_density`, `C_np_density_conditional`, `C_np_kernelsum` | Argument layouts, output requests and cleanup belong to the entry/lifecycle owner; registration is not the computation. |
| Prepared search state | [np.c](np.c); `np_regression_prepared_context_eval`, `np_density_prepared_context_destroy`, `C_np_density_conditional_prepared_destroy` | Follow prepare/refresh/evaluate/destroy as one lifecycle; do not assume state is local or reentrant. |
| Kernel sums and estimator orchestration | [jksum.c](jksum.c); `kernel_weighted_sum_np_ctx`, `kernel_weighted_sum_np_route`, `kernel_estimate_regression_categorical_tree_np` | Route selection, row traversal and final normalization are coupled; inspect the caller's operator and bandwidth topology. |
| LP objectives versus fit/apply | [jksum.c](jksum.c); `np_kernel_estimate_regression_categorical_ls_aic`, `np_conditional_density_cvls_lp_stream_ctx`, `np_regression_lp_apply_matrix` | Do not transplant an objective's delete-one or normalization rules into fitting/inference. |
| Basis conditioning | [jksum_lp_basis.h](jksum_lp_basis.h), [jksum_lp_basis.c](jksum_lp_basis.c); `np_lp_conditioned_basis_prepare` | The header specifies caller-owned matrices and column-space preservation. Public degree/basis preparation also involves np.c and jksum.c. |
| Resident LP accumulation | [jksum_lp_row.h](jksum_lp_row.h), [jksum_lp_row.c](jksum_lp_row.c); `NPLPDenseRowContext`, `np_lp_accumulate_dense_resident_row` | Adds to caller-owned moments/RHS. It does not own solve policy or MPI reductions. |
| Solve and rank policy | [jksum_lp_solve.h](jksum_lp_solve.h), [jksum_lp_solve.c](jksum_lp_solve.c); `NPLPSolveWorkspace`, `np_lp_solve_workspace_solve_response_ranked` | Source versus destructive LAPACK buffers, retained factors, reserve/clear and failure statuses are documented in the header. |
| Block dimensions | [jksum_block_plan.h](jksum_block_plan.h), [jksum_block_plan.c](jksum_block_plan.c); `NPJksumBlockPlan` | Plans describe checked dimensions; the caller allocates storage and performs communication. |
| Kernel/operator metadata | [kernel_registry.h](kernel_registry.h), [kernel_registry.c](kernel_registry.c), [kernel.c](kernel.c); `NPContinuousKernelRoute` | Family/order, operator and support metadata must remain consistent through the row consumer. |
| Continuous and beta rows | [continuous_kernel_row.h](continuous_kernel_row.h), [continuous_kernel_row.c](continuous_kernel_row.c), [continuous_kernel_row_gradient.c](continuous_kernel_row_gradient.c), [beta_kernel.c](beta_kernel.c), [beta_bandwidth.c](beta_bandwidth.c), [beta_scaled_row.h](beta_scaled_row.h) | Prepared beta allocation scopes and scaled-row/derivative contracts are explicit in the headers; do not substitute response-independent algebra at endpoints without proof. |
| Categorical product tiles | [categorical_profile_tile.h](categorical_profile_tile.h), [categorical_profile_tile.c](categorical_profile_tile.c), [categorical_profile_tile_r.c](categorical_profile_tile_r.c); `NPCategoricalProfileKernelSpec` | Native kernel codes/operators and the R caller determine units; similar-looking regression and density weights need not share normalization. |
| Geometry, support, NN radii | [tree.h](tree.h), [tree.c](tree.c), [tree_capability.h](tree_capability.h), [tree_capability.c](tree_capability.c), [kernelb.c](kernelb.c), [nn_radius_error.h](nn_radius_error.h) | Eligibility, traversal, bandwidth geometry and zero-radius reporting are distinct responsibilities. |
| Bounded Gaussian siblings | [jksum_gaussian_fixed.h](jksum_gaussian_fixed.h), [jksum_gaussian_fixed.c](jksum_gaussian_fixed.c), [jksum_gaussian_density.h](jksum_gaussian_density.h), [jksum_gaussian_density.c](jksum_gaussian_density.c) | Capability-selected arithmetic helpers rejoin their estimator owner; they do not define a second statistical method. |
| Hat/apply and batch projections | [reghat_fast.c](reghat_fast.c), [lp_batch_project.c](lp_batch_project.c), [npscoef_batch_solve.c](npscoef_batch_solve.c); `C_np_reghat_lp_matrix_fast`, `C_np_lp_batch_project_ranked`, `C_np_npscoef_batch_project` | R orchestrators choose matrix/apply/response-batch layouts. A projection is not automatically a reusable studentized statistic. |
| Residual and contrast inference | [regression_contract.h](regression_contract.h), [regression_residual.h](regression_residual.h), [regression_contrast.h](regression_contrast.h), [regression_inference_reuse.h](regression_inference_reuse.h), [conditional_lp_pair_se.c](conditional_lp_pair_se.c), [conditional_ann_direct.h](conditional_ann_direct.h), [conditional_kernel_moments.h](conditional_kernel_moments.h) | Residual preparation, paired contrasts and requested uncertainty remain distinct from mean-only work. Follow the enclosing owner in jksum.c. |
| Other inference/helpers | [ann_se.c](ann_se.c), [copula_se.c](copula_se.c), [quantile.c](quantile.c), [statmods.c](statmods.c), [entropy_gaussian_integrand.c](entropy_gaussian_integrand.c), [entropy_gaussian_summation.c](entropy_gaussian_summation.c) | Read the registered caller and estimator-specific target before reusing these routines. |
| Supporting infrastructure | [np_native_safety.h](np_native_safety.h), [kernelcv.c](kernelcv.c), [mat_vec.c](mat_vec.c), [nr.c](nr.c), [hash.c](hash.c), [bspline.c](bspline.c) | Includes checked sizes, CV callbacks, allocation/numerical utilities and basis support; filenames alone do not establish storage ownership. |

## Ownership and change checklist

- A pointer-bearing struct may be a borrowed view, an owning workspace, or
  process-retained state. Check its initializer, consumers and clear/destroy
  functions; a shallow copy does not transfer ownership.
- R-managed temporary storage, protected R results, manually allocated buffers
  and retained external state have different lifetimes. Check the actual
  cleanup/unwind path; this map does not certify every error path.
- Preserve array order, strides, tree permutations, fitted-row identity,
  deletion rules and family-specific kernel normalization.
- A change to a shared helper or selector needs its reverse caller inventory,
  including fitting, prediction, uncertainty and bootstrap consumers.
- Keep portable and accelerated siblings under the same arithmetic contract.
  The hot-engine alignment comment near `NP_HOT_ALIGN` in jksum.c records why
  moving function bodies or changing translation units is not cosmetic.
- Build inputs live in [Makevars.in](Makevars.in), [Makevars.win](Makevars.win)
  and `../configure`. This guide does not authorize changing those inputs.
- Update this map when an owner moves. Update declaration comments when
  ownership changes. Prefer stable symbols to line numbers and do not copy
  every detailed header contract here.

