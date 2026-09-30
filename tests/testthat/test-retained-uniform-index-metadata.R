# A small condition/metadata guard; broad pre/post coverage is external.
test_that("retained index uniform order is not a fresh user request", {
  old <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x1 = seq(-1, 1, length.out = 24L),
                  x2 = sin(seq_len(24L)))
  y <- sin(x$x1 + .3 * x$x2)
  expect_warning(b <- npindexbw(xdat = x, ydat = y,
    bws = c(1, .3, .6), ckertype = "uniform", ckerorder = 4L,
    bandwidth.compute = FALSE),
    "ignoring kernel order specified with uniform kernel type", fixed = TRUE)
  expect_identical(as.integer(b$ckerorder), 4L)
  expect_warning(fit <- npindex(b, se = FALSE), NA)
  expect_warning(H <- npindexhat(b, txdat = x, exdat = x,
                                output = "matrix"), NA)
  expect_warning(got <- npindexhat(b, txdat = x, exdat = x, y = y,
                                  output = "apply"), NA)
  index <- x$x1 + .3 * x$x2
  W <- abs(outer(index, index, "-")) < .6
  reference <- colSums(W * y) / colSums(W)
  expect_equal(as.numeric(fitted(fit)), reference, tolerance = 1e-12)
  expect_equal(as.numeric(got), reference, tolerance = 1e-12)
  expect_equal(as.numeric(H %*% y), reference, tolerance = 1e-12)
  context <- get(".np_retained_uniform_constructor",
                 asNamespace(getNamespaceName(environment(npindex))))
  expect_warning(context(function(...) warning("unrelated", call. = FALSE),
                         ckertype = "uniform"), "unrelated", fixed = TRUE)
  expect_error(context(function(...) stop("unrelated-error"),
                       ckertype = "uniform"), "unrelated-error", fixed = TRUE)
})
