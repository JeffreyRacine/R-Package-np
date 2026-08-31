test_that("R LP hat fallback uses the pristine intercept ridge anchor", {
  old <- options(matprod = "internal", np.tree = FALSE, np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  x <- data.frame(x = c(0, 1, 2, 3))
  y <- 1 + 2 * x$x
  evaluation <- data.frame(x = 0)

  for (bernstein in c(FALSE, TRUE)) {
    bw <- npregbw(
      xdat = x, ydat = y, bws = 0.1, bandwidth.compute = FALSE,
      bwtype = "fixed", regtype = "lp", degree = 1L,
      degree.select = "manual", basis = "glp",
      bernstein.basis = bernstein, ckertype = "epanechnikov"
    )

    options(matprod = "internal")
    H.r <- npreghat(
      bw, txdat = x, exdat = evaluation, output = "matrix", s = 0L
    )
    D.r <- npreghat(
      bw, txdat = x, exdat = evaluation, output = "matrix", s = 1L
    )
    options(matprod = "default")
    H.native <- npreghat(
      bw, txdat = x, exdat = evaluation, output = "matrix", s = 0L
    )
    D.native <- npreghat(
      bw, txdat = x, exdat = evaluation, output = "matrix", s = 1L
    )

    expect_equal(sum(H.r), 1.0, tolerance = 5e-15)
    expect_equal(sum(D.r), 0.0, tolerance = 5e-15)
    expect_equal(unclass(H.r), unclass(H.native), tolerance = 5e-15,
                 ignore_attr = TRUE)
    expect_equal(unclass(D.r), unclass(D.native), tolerance = 5e-15,
                 ignore_attr = TRUE)
  }
})

test_that("R LP hat ridge correction is coherent for NN bandwidths", {
  old <- options(matprod = "internal", np.tree = FALSE, np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  x <- data.frame(x = c(0, 0, 2, 3))
  y <- 1 + 2 * x$x
  evaluation <- data.frame(x = 0.1)

  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    bw <- npregbw(
      xdat = x, ydat = y, bws = 2, bandwidth.compute = FALSE,
      bwtype = bwtype, regtype = "lp", degree = 1L,
      degree.select = "manual", basis = "glp",
      bernstein.basis = FALSE, ckertype = "epanechnikov"
    )
    H <- npreghat(
      bw, txdat = x, exdat = evaluation, output = "matrix", s = 0L
    )
    expect_equal(sum(H), 1.0, tolerance = 5e-15)
    expect_equal(drop(H %*% rep.int(3.25, nrow(x))), 3.25,
                 tolerance = 2e-14)
  }
})
