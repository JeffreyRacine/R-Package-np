t4_capture_warning <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}

test_that("external LP empty rows are explicit across matrix, apply and prediction", {
  old <- options(np.messages = FALSE, np.tree = FALSE, warn = 0)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 64L))
  y <- sin(x$x)
  ex <- data.frame(x = c(0, 10, -.25))
  bw <- npregbw(xdat = x, ydat = y, bws = .3, bandwidth.compute = FALSE,
                regtype = "lp", degree = 2L, ckertype = "epanechnikov")
  for (tree in list(FALSE, TRUE, "auto")) {
    options(np.tree = tree)
    captured <- t4_capture_warning(npreghat(bw, txdat = x, exdat = ex))
    H <- captured$value
    expect_length(captured$warnings, 1L)
    expect_match(captured$warnings, "all computed kernel weights are zero")
    expect_match(captured$warnings, "1 external evaluation.*row")
    expect_true(all(is.na(H[2, ])))
    expect_true(all(is.finite(H[c(1, 3), ])))
    expect_null(attr(H, ".np.empty.rows", exact = TRUE))
    expect_true(is.na(attr(H, "ridge.used")[2]))
    supported <- npreghat(bw, txdat = x, exdat = ex[c(1, 3), , drop = FALSE])
    expect_identical(unname(H[c(1, 3), ]), unname(supported[, ]))
    for (rhs in list(y, cbind(first = y, second = y^2))) {
      a <- t4_capture_warning(npreghat(bw, txdat = x, exdat = ex,
                                       y = rhs, output = "apply"))
      expect_length(a$warnings, 1L)
      expected <- H %*% rhs
      if (is.null(dim(rhs))) expected <- as.vector(expected)
      expect_equal(a$value, expected, tolerance = 1e-10)
      expect_null(attr(a$value, ".np.empty.rows", exact = TRUE))
    }
    A <- t4_capture_warning(npreghat(bw, txdat = x, exdat = ex,
                                     y = y, output = "constraint"))
    expect_length(A$warnings, 1L)
    expect_equal(A$value[, ], t(H[, ]) * y, tolerance = 1e-10)
    P <- t4_capture_warning(predict(supported, newdata = ex))
    expect_length(P$warnings, 1L)
    expect_equal(P$value[, ], H[, ], tolerance = 1e-10)
  }
  options(np.tree = FALSE, warn = 2)
  expect_error(npreghat(bw, txdat = x, exdat = ex), "all computed kernel weights")
  options(warn = 0)
  expect_true(all(is.finite(npreghat(bw, txdat = x, exdat = ex[1, , drop = FALSE]))))
})

test_that("strict consumers and nonfinite or NN failures do not become partial results", {
  old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 64L))
  y <- sin(x$x)
  ex <- data.frame(x = c(0, 10))
  bw <- npregbw(xdat = x, ydat = y, bws = .3, bandwidth.compute = FALSE,
                regtype = "lp", degree = 2L, ckertype = "epanechnikov")
  expect_error(.npreghat_complete(bw, txdat = x, exdat = ex), "LP hat helper failed")
  expect_error(.npreghat_exact_lp_apply_from_regression_core(
    bw, x, matrix(y, ncol = 1), exdat = ex, degree = 2L, s = 0L),
    "LP apply helper failed")
  bad.y <- y
  bad.y[32] <- Inf
  expect_error(npreghat(bw, txdat = x, exdat = ex[1, , drop = FALSE],
                       y = bad.y, output = "apply"), "LP apply helper failed")
  tied <- data.frame(x = c(rep(0, 32), seq_len(32)))
  nn <- npregbw(xdat = tied, ydat = y, bws = 2, bandwidth.compute = FALSE,
                regtype = "lp", degree = 2L, bwtype = "generalized_nn")
  expect_error(npreghat(nn, txdat = tied, exdat = data.frame(x = 0)),
               "radius|radii")
  # Positive but tiny weights are not thresholded; exact computed zeros are
  # reported without pretending that their deeper cause is known.
  gaussian <- npregbw(xdat = x, ydat = y, bws = .3, bandwidth.compute = FALSE,
                      regtype = "lp", degree = 2L)
  far <- t4_capture_warning(npreghat(gaussian, txdat = x,
                                     exdat = data.frame(x = 1e4)))
  expect_length(far$warnings, 1L)
  expect_false(grepl("outside.*support|ties|underflow", far$warnings))
  expect_true(all(is.na(far$value)))
})

test_that("generalized NN keeps positive radii but classifies an empty product row", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 64L),
    f = factor(rep(c("a", "b"), 32L), levels = c("a", "b", "c")))
  y <- sin(x$x)
  bw <- npregbw(xdat = x, ydat = y, bws = c(5, 0),
    bandwidth.compute = FALSE, regtype = "lp", degree = 2L,
    bwtype = "generalized_nn")
  ex <- data.frame(x = c(0, 0), f = factor(c("a", "c"), levels = levels(x$f)))
  h <- t4_capture_warning(npreghat(bw, txdat = x, exdat = ex))
  a <- t4_capture_warning(npreghat(bw, txdat = x, exdat = ex,
                                  y = y, output = "apply"))
  expect_length(h$warnings, 1L)
  expect_length(a$warnings, 1L)
  expect_true(all(is.finite(h$value[1L, ])))
  expect_true(all(is.na(h$value[2L, ])))
  expect_equal(a$value, as.vector(h$value %*% y), tolerance = 1e-10)
})

test_that("worker-local empty-row completion is safe with warnings promoted to errors", {
  skip_on_cran()
  owned <- !.npRmpi_has_active_slave_pool(comm = 1L)
  if (owned) {
    if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  remote <- getFromNamespace("mpi.remote.exec", "npRmpi")
  before <- remote(Sys.getpid())
  replies <- remote(local({
    old <- options(np.messages = FALSE, warn = 2)
    on.exit(options(old))
    x <- data.frame(x = seq(-1, 1, length.out = 32L))
    bw <- npregbw(xdat = x, ydat = sin(x$x), bws = .3,
      bandwidth.compute = FALSE, regtype = "lp", degree = 2L,
      ckertype = "epanechnikov")
    H <- npreghat(bw, txdat = x, exdat = data.frame(x = c(0, 10)))
    all(is.finite(H[1L, ])) && all(is.na(H[2L, ])) &&
      is.null(attr(H, ".np.empty.rows", exact = TRUE))
  }))
  expect_true(length(replies) > 0L)
  expect_true(all(vapply(replies, isTRUE, logical(1))))
  expect_identical(remote(Sys.getpid()), before)
})
