test_that("conditional-density raw probes preserve resident exploration", {
  probe <- function(type, method) {
    ns <- asNamespace("np")
    old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE,
      npRmpi.autodispatch.disable = TRUE)
    on.exit(options(old), add = TRUE)
    x <- data.frame(x = c(rep(0, 20L), 1:4))
    y <- data.frame(y = sin((1:24) / 3))
    bw <- get("npcdensbw", ns)(xdat = x, ydat = y, bws = c(22, 22),
      bwtype = type, bwmethod = method, regtype = "lp", degree = 1L,
      nomad = FALSE, bandwidth.compute = FALSE)
    prep <- get(".npcdensbw_prepared_prepare_args", ns)(xdat = x, ydat = y,
      bws = bw, invalid.penalty = "baseline", degree.search = TRUE)
    names(prep)[names(prep) == "penalty_mode"] <- "penalty.mode"
    names(prep)[names(prep) == "penalty_multiplier"] <- "penalty.multiplier"
    prepare <- get("npPreparedObjectivePrepareConditionalDensity", ns)
    evaluate <- get("npPreparedObjectiveEvalConditionalDensity", ns)
    raw <- get("npPreparedObjectiveEvalConditionalDensityRaw", ns)
    destroy <- get("npPreparedObjectiveDestroyConditionalDensity", ns)
    capture <- function(expr) tryCatch(expr, error = function(e) conditionMessage(e))
    stopifnot(isTRUE(do.call(prepare, prep)))
    observed <- tryCatch({
      before <- evaluate(c(5, 22), 2L)
      invalid <- raw(c(5, 22), 2L)
      valid <- raw(c(22, 22), 1L)
      list(before = before, invalid = invalid, valid = valid,
        after = evaluate(c(5, 22), 2L),
        malformed = capture(raw(numeric(), 2L)),
        degree = capture(raw(c(22, 22), integer())))
    }, finally = destroy())
    observed$inactive <- capture(raw(c(22, 22), 1L))
    observed
  }
  for (type in c("generalized_nn", "adaptive_nn")) for (method in c("cv.ml", "cv.ls")) {
    observed <- probe(type, method)
    expect_identical(observed$before, observed$after)
    expect_true(abs(observed$before[[1L]]) < .Machine$double.xmax)
    expect_identical(observed$invalid[[1L]], -.Machine$double.xmax)
    expect_true(abs(observed$valid[[1L]]) < .Machine$double.xmax)
    expect_length(observed$invalid, 4L)
    expect_identical(observed$invalid[[2L]], 1)
    expect_identical(observed$valid[[2L]], 1)
    expect_match(observed$malformed, "bandwidth vector of unexpected length", fixed = TRUE)
    expect_match(observed$degree, "degree vector of unexpected length", fixed = TRUE)
    expect_match(observed$inactive, "state is not active", fixed = TRUE)
  }
})
