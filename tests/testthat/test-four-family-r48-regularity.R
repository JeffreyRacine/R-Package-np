test_that("working zero-variance disagreement is typed without hiding invariant failures", {
  A <- rbind(c(1, -1), c(0, 0))
  condition <- tryCatch(.np_plot_four_calibrate(A, c(1, 1), .05, 9999L,
    "working"), error = identity)
  expect_s3_class(condition, "np_plot_joint_unavailable")
  expect_match(conditionMessage(condition), "disagree about zero variance")
  expect_warning(value <- .np_plot_joint_try(.np_plot_four_calibrate(A,
    c(1, 1), .05, 9999L, "working")), "retaining any pointwise/Bonferroni")
  expect_null(value)

  exact <- tryCatch(.np_plot_four_calibrate(A, c(1, 1), .05, 9999L,
    "exact", exact.diagonal = TRUE), error = identity)
  expect_s3_class(exact, "error")
  expect_false(inherits(exact, "np_plot_joint_unavailable"))
  expect_error(.np_plot_joint_try(.np_plot_four_calibrate(A,
    c(1, 1), .05, 9999L, "exact", exact.diagonal = TRUE)), "disagree")
  expect_error(.np_plot_joint_try(.np_plot_four_calibrate(matrix(NA_real_, 2, 2),
    c(1, 1), .05, 9999L, "invalid")), "invalid four-family influence")
  expect_error(.np_plot_joint_try(.np_plot_four_calibrate(A[1, , drop=FALSE],
    2, .05, 9999L, "exact", exact.diagonal = TRUE)), "marginal variance")
  expect_error(.np_plot_lsq_joint(list(tau=.5, reg.bws=list()),
    data.frame(x=1:2), data.frame(x=1), .5, NA_real_),
    class="np_plot_joint_unavailable")
})

test_that("regular and algebraic-zero Gaussian rows retain exact calibration", {
  A <- rbind(c(1, -1), c(0, 0))
  actual <- .np_with_seed(71, .np_plot_four_calibrate(A, c(sqrt(2), 0),
    .05, 9999L, "exact", exact.diagonal = TRUE))
  oracle <- .np_with_seed(71, .np_plot_gaussian_calibrate(diag(c(1, 0)),
    alpha=.05, draws=9999L))
  expect_identical(actual$critical, oracle$critical)
  expect_equal(actual$grid.size, 2L)
  expect_error(.np_plot_four_calibrate(matrix(1e308, 2, 4), c(1, 1),
    .05, 9999L, "overflow"), class="np_plot_joint_unavailable")
  expect_error(.np_plot_four_calibrate(A * 1e-250, c(0, 0),
    .05, 9999L, "underflow"), class="np_plot_joint_unavailable")
})

test_that("projection-induced zero Jacobians cannot certify class certainty", {
  map <- function(raw, k, proper=TRUE)
    .np_plot_conmode_delta_row(raw, k, proper, sqrt(.Machine$double.eps))
  expect_error(map(c(-.02, .42, .6), 1L),
    class="np_plot_joint_unavailable")
  expect_error(map(c(1.08, -.04, -.04), 1L),
    class="np_plot_joint_unavailable")
  expect_error(map(c(1.08, -.04, -.04), 2L),
    class="np_plot_joint_unavailable")
  # The remaining active two-class face has a nonzero tangent map.
  active <- map(c(-.02, .42, .6), 2L)
  expect_equal(active$point, .41)
  expect_equal(active$jacobian, c(0, .5, -.5))
  raw <- map(c(-.02, .42, .6), 1L, FALSE)
  expect_identical(raw$point, -.02)
  expect_equal(raw$jacobian, c(1, 0, 0))
  interior <- map(c(.2, .3, .5), 2L)
  expect_equal(interior$point, .3)
  expect_equal(interior$jacobian, c(-.3, .7, -.3))
  # Ordinary zero response keeps its separate raw-factor certificate path.
  constant <- map(c(0, .4, .6), 1L)
  expect_equal(constant$point, 0)
  expect_equal(constant$jacobian, c(0, 0, 0))
  expect_identical(constant$zeros, 1L)
  expect_error(map(c(-1e-9, .4, .600000001), 1L),
    class="np_plot_joint_unavailable")
  expect_error(map(c(1+1e-9, -1e-9, 0), 1L),
    class="np_plot_joint_unavailable")
  expect_error(.np_plot_conmode_delta_row(c(-.1, .4, .7), 1L, TRUE, .2),
    class="np_plot_joint_unavailable")
  # A zero raw class still needs the caller's exact zero-response certificate;
  # another clipped class alone must not remove that certificate path.
  zero <- map(c(1+1e-9, -1e-9, 0), 3L)
  expect_identical(zero$zeros, 3L)
  expect_equal(zero$jacobian, c(0, 0, 0))
})
