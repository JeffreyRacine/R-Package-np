test_that("partial-linear public and cached-plot owners share rank validation", {
  ns <- asNamespace("np")
  set.seed(683)
  n <- 48L
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  Z <- data.frame(z = runif(n))
  y <- 1 + X$x1 - .5*X$x2 + sin(Z$z)
  bw <- npplregbw(xdat = X, ydat = y, zdat = Z,
                 bws = matrix(.3, 3L, 1L), bandwidth.compute = FALSE)
  a <- npplreg(bws = bw, txdat = X, tydat = y, tzdat = Z, se = FALSE)
  b <- npplreg(bws = bw, txdat = X, tydat = y, tzdat = Z, se = TRUE)
  expect_identical(a$mean, b$mean)
  expect_identical(coef(a), coef(b))
  expect_true(all(is.finite(vcov(b))))
  common <- get(".np_plot_plreg_apply_common", ns)(bw, X, Z)
  expect_length(common[["formation.error"]], 2L)
  expect_true(all(is.finite(common[["formation.error"]])))
  X$x2 <- X$x1
  bw <- npplregbw(xdat = X, ydat = y, zdat = Z,
                 bws = matrix(.3, 3L, 1L), bandwidth.compute = FALSE)
  for (s in c(FALSE, TRUE))
    expect_error(npplreg(bws = bw, txdat = X, tydat = y, tzdat = Z, se = s),
                 "rank deficient after smoothing", fixed = TRUE)
  expect_error(get(".np_plot_plreg_local_fit", ns)(bw, X, y, Z, X, Z),
               "rank deficient after smoothing", fixed = TRUE)
  common <- get(".np_plot_plreg_apply_common", ns)(bw, X, Z)
  state <- get(".np_plot_plreg_apply_eval_state", ns)(common, X, Z)
  expect_error(get(".np_plot_plreg_apply_from_state", ns)(state, y),
               "rank deficient after smoothing", fixed = TRUE)
})
