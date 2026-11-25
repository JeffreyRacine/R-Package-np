test_that("fixed uniform CDF primitive, process and range agree with signed integration", {
  x <- c(-.9, -.4, -.15, .22, .57, .82)
  g <- data.frame(x = c(-2, -.38, .29, 2))
  state <- list(type = "fixed", x = data.frame(x), n = length(x),
                q = 1L, s = 0L, numeric = TRUE, h = .25, category = list(),
                max.bytes = 256*1024^2, metadata = list(kernel = "uniform2"))
  primitive <- function(t) integrate(function(u) (9 - 15*u^2)/8, -1, max(-1, min(1, t)))$value
  rows <- .np_udist_bc_rows(state, g)
  t <- outer(g$x, x, "-") / state$h
  J <- matrix(vapply(t, primitive, 0.0), nrow(g))
  H <- (1 + pmin(pmax(t, -1), 1))/2
  expect_equal(rows$corrected, J, tolerance = 1e-13)
  expect_equal(rows$original, H, tolerance = 1e-13)
  expect_equal(.np_udist_bc_rows(state, g, FALSE)$corrected, H, tolerance = 1e-13)
  expect_equal(rows$center, rowMeans(J), tolerance = 1e-13)
  factor <- sweep(J, 1L, rowMeans(J), "-") / length(x)
  expect_equal(tcrossprod(rows$factor), tcrossprod(factor), tolerance = 1e-13)
  expect_equal(unname(.np_udist_marginal_range(1L, "uniform2")),
               vapply(c(-1, 1)*sqrt(3/5), primitive, 0.0), tolerance = 1e-13)
  for (law in c("inid", "fixed", "geom"))
    expect_error(.np_udist_joint_inference(state, g, law = law), "uniform CDF")
  old <- options(np.largeh = TRUE)
  on.exit(options(old))
  state$h <- diff(range(c(x, g$x)))
  expect_error(.np_udist_bc_rows(state, g), "large-h")
  options(np.largeh = FALSE)
  expect_equal(.np_udist_bc_rows(state, g)$center,
               rowMeans(matrix(vapply(outer(g$x, x, "-") / state$h, primitive, 0.0), nrow(g))),
               tolerance = 1e-13)
})
