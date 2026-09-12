test_that("conditional categorical demand keeps one existing endpoint pair", {
  pkg <- "np"
  helper <- getFromNamespace("npConditionalCategoricalFirstDifferences", pkg)
  hat <- getFromNamespace("npcdenshat", pkg)
  x <- data.frame(
    x = c(-0.2, 0.1, 0.3),
    u = factor(c("b", "a", "c")),
    o = ordered(c("low", "middle", "high"),
                levels = c("low", "middle", "high")))
  y <- data.frame(y = seq_len(nrow(x)))
  bws <- list(xndim = 3L, ixcon = c(TRUE, FALSE, FALSE),
              ixuno = c(FALSE, TRUE, FALSE), ixord = c(FALSE, FALSE, TRUE),
              xdati = list(icon = c(TRUE, FALSE, FALSE)))
  calls <- list()
  local_mocked_bindings(
    .np_conditional_higher_hat = function(hat.args, cdf, base.rows,
                                          allow.external, return.norm) {
      expect_false(allow.external)
      expect_false(return.norm)
      calls[[length(calls) + 1L]] <<- hat.args$exdat
      as.double(hat.args$exdat$x + 2 * as.integer(hat.args$exdat$u) +
                  3 * as.integer(hat.args$exdat$o))
    }, .package = pkg)
  args <- list(hat.fun = hat, bws = bws, txdat = x, tydat = y,
               exdat = x, eydat = y, where = "demand test",
               allow.external = FALSE, .np.defer.empty.rows = TRUE)
  full <- do.call(helper, args)
  full.calls <- calls
  expect_length(full.calls, 4L)
  for (j in 2:3) {
    calls <- list()
    one <- do.call(helper, c(args, list(gradient.target = j)))
    expect_identical(dim(one), dim(full))
    expect_identical(one[, j], full[, j])
    expect_true(all(is.na(one[, -j, drop = FALSE])))
    expect_identical(calls, full.calls[(2L * j - 3L):(2L * j - 2L)])
  }
  expect_error(do.call(helper, c(args, list(gradient.target = 1L))),
               "categorical")
  expect_error(do.call(helper, c(args, list(gradient.target = 2.5))),
               "invalid conditional gradient coordinate")
})

test_that("private conditional demand preserves native target and rejects uncertainty", {
  pkg <- "np"
  eval.raw <- getFromNamespace(".np_conditional_eval_selected", pkg)
  bw.fun <- getFromNamespace("npcdensbw", pkg)
  with.local <- if (pkg == "npRmpi")
    getFromNamespace(".npRmpi_with_local_regression", pkg) else
      function(expr) force(expr)
  eval.fun <- function(...) with.local(eval.raw(...))
  i <- seq_len(36L)
  x <- data.frame(x = (i - 18) / 18,
                  f = factor(rep(c("a", "b"), length.out = length(i))),
                  z = sin(i * 0.7))
  y <- data.frame(y = sin(i * 0.31) + x$x)
  bws <- with.local(bw.fun(xdat = x, ydat = y,
                 bws = c(0.75, 0.8, 0.2, 0.8),
                 bandwidth.compute = FALSE, regtype = "ll"))
  args <- list(bws = bws, xdat = x, ydat = y,
               exdat = x[c(4L, 11L, 24L), , drop = FALSE],
               eydat = y[c(4L, 11L, 24L), , drop = FALSE],
               gradients = TRUE, se = FALSE)
  full <- do.call(eval.fun, args)
  for (j in seq_len(ncol(x))) {
    one <- do.call(eval.fun, c(args, list(gradient.target = j)))
    expect_identical(one$congrad[, j], full$congrad[, j])
    expect_identical(one$condens, full$condens)
    expect_identical(attr(one, ".np.empty.rows", exact = TRUE),
                     attr(full, ".np.empty.rows", exact = TRUE))
    expect_identical(attr(one, ".np.empty.base.rows", exact = TRUE),
                     attr(full, ".np.empty.base.rows", exact = TRUE))
  }
  bad <- args
  bad$se <- TRUE
  expect_error(do.call(eval.fun, c(bad, list(gradient.target = 1L))),
               "gradients = TRUE.*se = FALSE.*proper = FALSE")
  bad <- args
  bad$gradients <- FALSE
  expect_error(do.call(eval.fun, c(bad, list(gradient.target = 1L))),
               "gradients = TRUE.*se = FALSE.*proper = FALSE")
  bad <- args
  bad$proper <- TRUE
  expect_error(do.call(eval.fun, c(bad, list(gradient.target = 1L))),
               "gradients = TRUE.*se = FALSE.*proper = FALSE")
  expect_error(do.call(eval.fun,
    c(args, list(categorical.effects = FALSE, gradient.target = 2L))),
    "categorical.effects = TRUE")
})

test_that("smooth bootstrap forwards its physical continuous gradient demand", {
  pkg <- "np"
  helper <- getFromNamespace(".np_plot_conditional_gradient_smooth_boot", pkg)
  x <- data.frame(u = factor(c("a", "b", "a")), x = c(-.2, .1, .3),
                  o = ordered(c("low", "middle", "high")))
  y <- data.frame(y = 1:3)
  bws <- list(type = "fixed", bandwidth = list(x = c(.2, .8, .2), y = .8),
              xdati = list(icon = c(FALSE, TRUE, FALSE)),
              ydati = list(icon = TRUE), cxkertype = "gaussian", cxkerorder = 2L,
              cykertype = "gaussian", cykerorder = 2L)
  calls <- 0L
  target <- 2L
  local_mocked_bindings(
    .np_plot_conditional_pilot_prepare = function(xdat, ydat, ...)
      list(x=xdat, y=ydat),
    .np_plot_pilot_draw_side = function(side, idx) side[idx, , drop = FALSE],
    .np_plot_conditional_pilot_reference = function(pilot, exdat, ...)
      rep(0, nrow(exdat)),
    npConditionalRegEngineSpec = function(...) list(reg.engine="lp"),
    npConditionalGradientOrder = function(...) 2L,
    .np_plot_conditional_eval = function(exdat, gradients, gradient.order,
      lp.first.se.demand, cat.se.demand, se, gradient.target, ...) {
      calls <<- calls + 1L
      expect_true(gradients)
      expect_identical(gradient.order, 2L)
      expect_false(lp.first.se.demand)
      expect_false(cat.se.demand)
      expect_false(se)
      expect_identical(gradient.target, target)
      list(congrad = matrix(rep(c(10, 20, 30), each = nrow(exdat)),
                            nrow = nrow(exdat)))
    }, .package = pkg)
  withr::local_options(np.messages = FALSE, np.plot.progress = FALSE)
  args <- list(xdat = x, ydat = y, exdat = x[1:2, , drop = FALSE],
    eydat = y[1:2, , drop = FALSE], bws = bws, cdf = FALSE,
    gradient.index = 2L, gradient.order = 2L, plot.errors.boot.method = "inid",
    plot.errors.boot.blocklen = 2L, plot.errors.boot.num = 1L, progress.label = NULL)
  out <- do.call(helper, args)
  expect_identical(calls, 2L)
  expect_identical(out$t0, c(20, 20))
  expect_identical(out$t, matrix(20, nrow = 1L, ncol = 2L))
  args$gradient.index <- target <- 1L
  out <- do.call(helper, args)
  expect_identical(out$t0, c(10, 10))
  expect_identical(out$t, matrix(10, nrow = 1L, ncol = 2L))
  expect_identical(calls, 4L)
  args$gradient.index <- 0L
  expect_error(do.call(helper, args), "invalid conditional gradient coordinate")
  expect_identical(calls, 4L)
})
