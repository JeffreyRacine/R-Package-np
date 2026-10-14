# One focused regression contract: exact same-sample ratio cancellation in
# new copula simultaneous requests. The matrix fixture supplies raw marginal
# kernels independently; it does not change package namespace bindings.
copula_constant_ratio_fixture <- function(weights, point = NULL, se = NULL) {
  pkg <- unname(getNamespaceName(environment(npcopula)))
  owner <- getFromNamespace(".np_plot_copula_influence", pkg)
  scope <- new.env(parent = environment(owner))
  scope$.npcopula_marginal_bw <- function(bws, data, j, ...) list(index = j)
  scope$.np_plot_density_matrix <- function(kbw, x, ex, operator)
    weights[[kbw$index]][as.integer(ex[[1L]]), , drop = FALSE]
  environment(owner) <- scope
  n <- ncol(weights[[1L]]); m <- nrow(weights[[1L]]); p <- length(weights)
  data <- as.data.frame(setNames(rep(list(seq_len(n)), p), paste0("x", seq_len(p))))
  ex <- as.data.frame(setNames(rep(list(seq_len(m)), p), names(data)))
  means <- lapply(weights, rowMeans)
  if (is.null(point)) point <- rowMeans(Reduce(`*`, weights))/Reduce(`*`, means)
  if (is.null(se)) se <- rep(.3, m)
  object <- list(bws = list(type = "fixed", xnames = names(data)),
    target = "density", copula = point)
  owner(object, data, ex, se)
}

test_that("copula ratio cancellation uses exact raw donor identity row by row", {
  withr::local_options(np.messages = FALSE, np.plot.progress = FALSE)
  u <- c(1, 2, 3, 4); v <- c(2, 4, 1, 3)
  left <- rbind(u, rep(-2, 4), rep(-2, 4), u, u)
  right <- rbind(rep(3, 4), v, rep(3, 4), v, c(3, 3, 3, 3 + 1e-6))
  got <- copula_constant_ratio_fixture(list(left, right))
  expect_identical(got$A[1:3, , drop = FALSE], matrix(0, 3, 4))
  expect_identical(got$simultaneous.se, c(0, 0, 0, .3, .3))
  expect_identical(got$estimate[1:3], rep(1, 3))
  # A coincidentally unit ratio is not a constant random ratio. Raw quotient
  # differentiation supplies the independent common-donor reference at row 4.
  a <- mean(u); b <- mean(v); joint <- mean(u * v); ratio <- joint/(a*b)
  oracle <- (u*v-joint)/(4*a*b) -
    ratio*((u-a)/(4*a)+(v-b)/(4*b))
  expect_equal(got$A[4, ], oracle, tolerance = 1e-14)
  expect_gt(sum(got$A[4, ]^2), 0)
  # No numerical-nearness predicate: the nearly constant fifth row stays live.
  expect_gt(sum(got$A[5, ]^2), 0)
  wrong.point <- got$estimate
  wrong.point[1] <- 1 + 1e-5
  expect_error(copula_constant_ratio_fixture(list(left, right), point = wrong.point),
    "constant copula ratio does not reproduce the plotted target", fixed = TRUE)
  expect_error(copula_constant_ratio_fixture(
    list(matrix(u, 1), matrix(0, 1, 4)), point = 1),
    "zero or nonfinite marginal denominator", fixed = TRUE)
})

test_that("copula ratio cancellation counts every marginal in higher dimensions", {
  withr::local_options(np.messages = FALSE, np.plot.progress = FALSE)
  u <- c(1, 2, 3, 4); v <- c(2, 4, 1, 3)
  got <- copula_constant_ratio_fixture(list(rbind(u, u),
    rbind(v, rep(-2, 4)), rbind(rep(3, 4), rep(3, 4))))
  expect_gt(sum(got$A[1, ]^2), 0)
  expect_identical(got$A[2, ], rep(0, 4))
  expect_identical(got$simultaneous.se, c(.3, 0))
  two <- copula_constant_ratio_fixture(list(matrix(u, 1), matrix(v, 1)))
  expect_equal(got$A[1, ], as.vector(two$A), tolerance = 1e-14)
  expect_null(two$simultaneous.se)
})

test_that("copula cancellation excludes every marginal at or below the public floor", {
  withr::local_options(np.messages = FALSE, np.plot.progress = FALSE)
  eps <- .Machine$double.eps
  margins <- eps * c(1 - eps, 1, 1 + eps, -(1 - eps), -1, -(1 + eps))
  point <- abs(margins) / pmax(abs(margins), eps)
  expected.scale <- c(.3, .3, 0, .3, .3, 0)
  # Even near-floor public points roundoff-close to one cannot certify a zero.
  constant <- matrix(rep(margins, each = 4), ncol = 4, byrow = TRUE)
  ordinary <- matrix(rep(c(1, 2, 3, 4), 6), ncol = 4, byrow = TRUE)
  got <- copula_constant_ratio_fixture(list(ordinary, constant), point = point)
  expect_identical(got$simultaneous.se, expected.scale)
  # A varying mean exactly at the floor is a one-sided kink, not an identity.
  varying <- outer(margins, c(0, 1, 1, 2))
  got <- copula_constant_ratio_fixture(list(varying, matrix(3, 6, 4)), point = point)
  expect_identical(got$simultaneous.se, expected.scale)
  expect_error(copula_constant_ratio_fixture(
    list(matrix(c(1, 2, 3, 4), 1), matrix(eps/8, 1, 4)), point = 1/8),
    "copula influence does not reproduce the plotted target", fixed = TRUE)
})

test_that("copula cancellation does not hide nonfinite intermediate arithmetic", {
  withr::local_options(np.messages = FALSE, np.plot.progress = FALSE)
  # Diagnostic-only representability boundary: finite inputs and nonzero mean,
  # but division by that mean overflows. It must never become a certified zero.
  huge <- .Machine$double.xmax
  expect_error(copula_constant_ratio_fixture(
    list(matrix(c(huge, -huge, 1, 1), 1), matrix(1, 1, 4)), point = 1),
    "invalid four-family influence factor", fixed = TRUE)
})

test_that("copula cancellation leaves NN and distribution owners outside its cone", {
  withr::local_options(np.messages = FALSE, np.plot.progress = FALSE)
  pkg <- unname(getNamespaceName(environment(npcopula)))
  owner <- getFromNamespace(".np_plot_copula_influence", pkg)
  scope <- new.env(parent = environment(owner))
  environment(owner) <- scope
  data <- data.frame(x = 1:4, y = 4:1)
  ex <- data[1, , drop = FALSE]
  sentinel <- list(A = matrix(c(-.1, -.05, .05, .1), 1), estimate = .25,
    simultaneous.se = .2, exact.diagonal = TRUE)
  calls <- character()
  scope$.np_plot_copula_ratio_influence <- function(...) {
    calls <<- c(calls, "NN density"); sentinel
  }
  scope$.np_plot_ann_bounded_distribution_factor <- function(...) {
    calls <<- c(calls, "ANN distribution"); sentinel
  }
  scope$.np_plot_four_budget <- function(...) stop("entered fixed density allocation")
  for (type in c("generalized_nn", "adaptive_nn")) {
    object <- list(bws = list(type = type, icon = c(TRUE, TRUE)),
      target = "density", copula = .25)
    expect_identical(owner(object, data, ex, .2), sentinel)
  }
  object <- list(bws = list(type = "adaptive_nn", icon = c(TRUE, TRUE)),
    target = "distribution", copula = .25)
  expect_identical(owner(object, data, ex, .2), sentinel)
  expect_identical(calls, c("NN density", "NN density", "ANN distribution"))
  rm(".np_plot_four_budget", envir = scope)
  scope$.np_plot_density_matrix <- function(...) stop("entered density kernel owner")
  scope$.np_plot_density_unconditional <- function(...) sentinel
  for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
    object <- list(bws = list(type = type, icon = c(FALSE, FALSE)),
      target = "distribution", copula = .25)
    got <- owner(object, data, ex, .2)
    expect_identical(got$A, sentinel$A)
    expect_identical(got$estimate, sentinel$estimate)
    expect_false(got$exact.diagonal)
    expect_null(got$simultaneous.se)
  }
})

test_that("constant categorical copula margins retain marginal bands and allow joint bands", {
  pkg <- unname(getNamespaceName(environment(npcopula)))
  if (identical(pkg, "npRmpi")) {
    skip_on_cran()
    if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  withr::local_options(np.messages = FALSE, np.plot.progress = FALSE)
  i <- seq_len(96L)
  data <- data.frame(x = sin(i*sqrt(2)) + .013*i/length(i),
    category = ordered(c("low", "middle", "high")[1L + i %% 3L],
      levels = c("low", "middle", "high")))
  bw <- npudensbw(dat = data, bws = c(.48, 1), bandwidth.compute = FALSE,
    ckertype = "gaussian", okertype = "racineliyan")
  object <- npcopula(bws = bw, data = data,
    u = data.frame(u1 = c(.15, .3, .5, .7, .85), u2 = c(.25, .5, .75, .85, .95)),
    target = "density", se = TRUE, n.quasi.inv = 100L)
  before <- object
  method <- getS3method("plot", "npcopula", envir = asNamespace(pkg))
  request <- function(band) method(object, output = "data",
    errors = "asymptotic", band = band)
  pointwise <- request("pointwise")
  bonferroni <- request("bonferroni")
  joint <- request("simultaneous")
  all.bands <- request("all")
  expect_identical(object, before)
  expect_identical(joint$lower, joint$center)
  expect_identical(joint$upper, joint$center)
  expect_identical(all.bands$simultaneous.lower, all.bands$center)
  expect_identical(all.bands$simultaneous.upper, all.bands$center)
  expect_identical(all.bands$pointwise.lower, pointwise$lower)
  expect_identical(all.bands$pointwise.upper, pointwise$upper)
  expect_identical(all.bands$bonferroni.lower, bonferroni$lower)
  expect_identical(all.bands$bonferroni.upper, bonferroni$upper)
  expect_identical(attr(joint, "band.calibration")$critical, 0)
})

