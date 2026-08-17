np_prepared_migration_root <- function() {
  dirname(npRmpi_test_source_path("DESCRIPTION"))
}

np_prepared_migration_text <- function(root, paths) {
  paste(unlist(lapply(
    file.path(root, paths),
    readLines,
    warn = FALSE
  )), collapse = "\n")
}

test_that("prepared-objective migration closure cannot restore retired owners", {
  root <- np_prepared_migration_root()
  native <- np_prepared_migration_text(root, c("src/np.c", "src/np_init.c"))
  r_code <- np_prepared_migration_text(root, c(
    "R/np.density.bw.R",
    "R/np.distribution.bw.R",
    "R/np.regression.bw.R",
    "R/np.condensity.bw.R",
    "R/np.condistribution.bw.R"
  ))

  retired <- c(
    "np_density_nomad_native_eval_once",
    "np_distribution_nomad_native_eval_once",
    "np_regression_nomad_native_eval_once",
    "np_distribution_conditional_nomad_native_eval_once",
    "NPRegressionNomadShadowCtx",
    "NPConditionalDensityNomadShadowCtx",
    "C_np_density_nomad_native_fixed_eval",
    "C_np_distribution_nomad_native_fixed_eval",
    "C_np_regression_nomad_shadow_prepare",
    "C_np_regression_nomad_shadow_eval",
    "C_np_regression_nomad_shadow_native_search",
    "C_np_regression_nomad_shadow_clear",
    "C_np_density_conditional_nomad_shadow_prepare",
    "C_np_density_conditional_nomad_shadow_eval",
    "C_np_density_conditional_nomad_shadow_native_search",
    "C_np_density_conditional_nomad_shadow_fixed_native_search",
    "C_np_density_conditional_nomad_shadow_clear"
  )
  for (symbol in retired) {
    expect_false(grepl(symbol, native, fixed = TRUE), info = symbol)
    expect_false(grepl(symbol, r_code, fixed = TRUE), info = symbol)
  }

  retired_r <- c(
    ".npregbw_nomad_shadow_",
    ".npcdensbw_nomad_shadow_",
    "npRmpiNomadShadowSearchConditionalDistribution"
  )
  for (symbol in retired_r)
    expect_false(grepl(symbol, r_code, fixed = TRUE), info = symbol)
})
