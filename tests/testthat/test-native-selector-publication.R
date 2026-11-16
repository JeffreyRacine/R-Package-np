test_that("selector finish preserves public names, calls and owners", {
  local <- new.env(parent = baseenv())
  caller <- new.env(parent = baseenv())
  call <- quote(npregbw(xdat = .np_native_xdat,
                       ydat = .np_native_ydat, bws = .np_native_bws))
  environment(call) <- local
  x <- list(call = call, xnames = ".np_native_xdat",
            ynames = ".np_native_ydat", bw = .5)
  expressions <- list(xdat = quote(x), ydat = quote(`response name`), bws = .5)
  answer <- .np_bw_selector_finish(x, expressions, value.env = local,
                                   caller_env = caller)
  expected <- quote(npregbw(xdat = x, ydat = `response name`, bws = .5))
  environment(expected) <- caller
  expect_identical(answer$call, expected)
  expect_identical(answer$xnames, "x")
  expect_identical(answer$ynames, "response name")
  expect_identical(environment(x$call), local)
  expect_identical(.np_bw_selector_finish(answer, expressions), answer)
  expect_null(.np_bw_selector_finish(NULL, expressions))
  expect_identical(.np_bw_selector_finish(1, expressions), 1)
})

test_that("native selector certificate does not infer ownership from names", {
  call <- quote(npregbw.default(xdat = .np_native_xdat, ydat = .np_native_ydat))
  expect_null(.npRmpi_bw_native_selector_plan(call, "npregbw.default",
    npregbw.default, environment(), parent.frame()))
  expect_null(.npRmpi_bw_native_selector_plan(call, "npregbw.NULL",
    npregbw.NULL, environment(), parent.frame()))
})

test_that("named native regression fits use an unchanged canonical bandwidth", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI runtime not available")
  on.exit(close_mpi_slaves(), add = TRUE)
  ns <- asNamespace("npRmpi")
  old <- options(l1.test.native.valid = 0L)
  on.exit(options(old), add = TRUE)
  trace(".npRmpi_lease_resolve", exit = quote({
    answer <- returnValue()
    if (identical(consumer, "npreg") && !is.null(answer))
      options(l1.test.native.valid = getOption("l1.test.native.valid", 0L) + 1L)
  }), where = ns, print = FALSE)
  on.exit(untrace(".npRmpi_lease_resolve", where = ns), add = TRUE)
  i <- seq_len(24L); x <- data.frame(x = sin(i*.7)); y <- cos(i*.3)
  fit <- npreg(txdat = x, tydat = y, bws = .5,
               bandwidth.compute = FALSE, se = FALSE)
  expect_gt(getOption("l1.test.native.valid"), 0L)
  expect_identical(fit$bws$ynames, "y")
  expect_identical(fit$bws$call$xdat, quote(x))
  expect_identical(fit$bws$call$ydat, quote(y))
  options(l1.test.native.valid = 0L)
  other <- npreg(txdat = x, tydat = y + 0, bws = .5,
                 bandwidth.compute = FALSE, se = FALSE)
  expect_identical(getOption("l1.test.native.valid"), 0L)
  expect_identical(fitted(other), fitted(fit))
})
