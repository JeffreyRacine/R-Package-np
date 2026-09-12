# Shared contract for beta uncertainty: fixed/GNN retain their contribution
# oracle; finite-bound ANN publishes unavailable SEs, not a fixed-radius formula.
expect_beta_se_fit <- function(bwtype, code) {
  state <- new.env(parent=emptyenv())
  state$warnings <- character()
  fit <- withCallingHandlers(force(code), warning=function(w) {
    state$warnings <- c(state$warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  if (identical(bwtype, "adaptive_nn")) {
    expect_length(state$warnings, 1L)
    expect_match(state$warnings,
      "ANN uncertainty with finite kernel bounds is not yet implemented",
      fixed=TRUE)
  } else {
    expect_length(state$warnings, 0L)
  }
  fit
}

expect_beta_se <- function(fit, bwtype, expected, tolerance) {
  if (identical(bwtype, "adaptive_nn")) expected[] <- NA_real_
  expect_equal(se(fit), expected, tolerance=tolerance)
}
