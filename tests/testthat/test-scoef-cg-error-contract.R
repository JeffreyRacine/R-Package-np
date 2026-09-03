test_that("CG diagnostics advise only native optim nonfinite-parameter failures", {
  original <- tryCatch(stats::optim(c(Inf, 0), function(x) sum(x*x), method = "CG"),
                       error = identity)
  caught <- tryCatch(.npscoef_cg_diagnostic(
    stats::optim(c(Inf, 0), function(x) sum(x*x), method = "CG"), "CG"),
    error = identity)
  expect_s3_class(caught, "np_cg_numerical_error")
  expect_match(conditionMessage(caught), 'Try optim.method = "Nelder-Mead"', fixed = TRUE)
  expect_identical(conditionMessage(caught$parent), conditionMessage(original))
  expect_identical(conditionCall(caught), conditionCall(caught$parent))

  for (message in c("unrelated objective failure",
                   gettext("non-finite value supplied by optim", domain = "R-stats"))) {
    for (original.call in list(NULL, quote(objective(point)))) {
      sentinel <- errorCondition(message, class = "objective_test_error", call = original.call)
      caught <- tryCatch(.npscoef_cg_diagnostic(
        stats::optim(c(1, 2), function(x) stop(sentinel), method = "CG"), "CG"),
        error = identity)
      expect_identical(caught, sentinel)
    }
  }
  sentinel <- errorCondition("non-finite value supplied by optim")
  for (method in c("Nelder-Mead", "BFGS", "CG"))
    expect_identical(tryCatch(.npscoef_cg_diagnostic(stop(sentinel), method),
                             error = identity), sentinel)
  interrupt <- structure(list(message = "interrupt", call = NULL),
                         class = c("interrupt", "condition"))
  expect_identical(tryCatch(.npscoef_cg_diagnostic(stop(interrupt), "CG"),
                           interrupt = identity), interrupt)
})

test_that("CG diagnostic setup leaves successful optim results and RNG unchanged", {
  objective <- function(x) sum((x-c(.2, -.4))^2)
  for (method in c("Nelder-Mead", "BFGS", "CG")) {
    set.seed(9303)
    before <- .Random.seed
    original <- stats::optim(c(1, 2), objective, method = method)
    wrapped <- .npscoef_cg_diagnostic(stats::optim(c(1, 2), objective, method = method),
                                     method)
    expect_identical(wrapped, original)
    expect_identical(.Random.seed, before)
  }
  expect_warning(expect_identical(.npscoef_cg_diagnostic({
    warning("unchanged warning"); 7L
  }, "CG"), 7L), "unchanged warning", fixed = TRUE)
})

test_that("every smooth-coefficient optim call shares the diagnostic boundary", {
  counts <- c(optim = 0L, wrapped = 0L)
  visit <- function(node, inside = FALSE) {
    if (missing(node))
      return(invisible(NULL))
    if (!is.call(node) && !is.pairlist(node))
      return(invisible(NULL))
    if (is.call(node)) {
      if (identical(node[[1L]], as.name(".npscoef_cg_diagnostic")))
        inside <- TRUE
      if (identical(node[[1L]], as.name("optim"))) {
        counts[["optim"]] <<- counts[["optim"]] + 1L
        if (inside)
          counts[["wrapped"]] <<- counts[["wrapped"]] + 1L
      }
    }
    for (part in as.list(node))
      visit(part, inside)
  }
  visit(body(npscoefbw.scbandwidth))
  expect_gte(counts[["optim"]], 3L)
  expect_identical(counts[["wrapped"]], counts[["optim"]])
})
