test_that("explicit ridges are rejected only by owners that cannot honor them", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    spawn_mpi_slaves()
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  set.seed(7439)
  x <- data.frame(x = runif(24), z = runif(24))
  y <- x$x + x$z^2
  b <- npregbw(xdat = x, ydat = y, bws = c(.5, .5),
    regtype = "lp", degree = c(2L, 2L), bandwidth.compute = FALSE)
  e <- x[1:3, ]
  for (out in c("matrix", "apply", "constraint")) {
    for (s in list(c(0L, 0L), c(1L, 0L))) {
      expect_error(npreghat(b, exdat = e, y = y, output = out,
        s = s, ridge = .1), "nonzero 'ridge' is not supported", fixed = TRUE)
      a <- npreghat(b, exdat = e, y = y, output = out, s = s)
      z <- npreghat(b, exdat = e, y = y, output = out, s = s, ridge = 0)
      expect_equal(as.vector(a), as.vector(z), tolerance = 0)
    }
  }
  h <- npreghat(b, exdat = e, s = c(1L, 1L), ridge = .1)
  expect_true(all(attr(h, "ridge.used") >= .1))
  expect_equal(as.vector(predict(h, exdat = e)), as.vector(h), tolerance = 1e-12)
  expect_gt(max(abs(h - npreghat(b, exdat = e, s = c(1L, 1L)))), .01)
  expect_error(npreghat(b, y = y, output = "apply",
                       leave.one.out = TRUE, ridge = .1), "nonzero 'ridge'")
  for (bad in list(NA_real_, Inf, -1, c(0, 1), "0"))
    expect_error(npreghat(b, ridge = bad), "finite nonnegative")
})
