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
