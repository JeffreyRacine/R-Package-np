npplreg_audit_ready <- function() {
  if (exists("spawn_mpi_slaves", mode = "function"))
    return(spawn_mpi_slaves())
  TRUE
}

npplreg_audit_ns <- function()
  getNamespaceName(environment(npplreg))

npplreg_audit_fixture <- function(n = 18L) {
  z <- seq(-1, 1, length.out = n)
  x <- 0.4 * z + seq(-0.25, 0.25, length.out = n)
  y <- 1.25 * x + sin(2 * z) + 0.05 * cos(seq_len(n))
  list(xdat = data.frame(x = x), zdat = data.frame(z = z), ydat = y)
}

npplreg_audit_bw <- function(dat, bw = 0.55, bwtype = "fixed") {
  npplregbw(
    xdat = dat$xdat,
    zdat = dat$zdat,
    ydat = dat$ydat,
    bws = matrix(rep(bw, 2L), nrow = 2L),
    bwtype = bwtype,
    bandwidth.compute = FALSE
  )
}

test_that("npplreg direct default route preserves vector response and explicit-bw equivalence", {
  if (!npplreg_audit_ready()) skip("Could not initialize npplreg audit context")

  dat <- npplreg_audit_fixture()
  bw.mat <- matrix(rep(0.55, 2L), nrow = 2L)
  direct <- npplreg(
    bws = bw.mat,
    txdat = dat$xdat,
    tydat = dat$ydat,
    tzdat = dat$zdat,
    bandwidth.compute = FALSE
  )
  bw <- npplregbw(
    xdat = dat$xdat,
    ydat = dat$ydat,
    zdat = dat$zdat,
    bws = bw.mat,
    bandwidth.compute = FALSE
  )
  explicit <- npplreg(bws = bw)
  auto <- npplreg(
    txdat = dat$xdat,
    tydat = dat$ydat,
    tzdat = dat$zdat,
    nmulti = 1L,
    tol = 0.5,
    ftol = 0.5
  )

  expect_s3_class(direct, "plregression")
  expect_s3_class(auto, "plregression")
  expect_equal(as.numeric(direct$mean), as.numeric(explicit$mean), tolerance = 1e-10)
  expect_equal(as.numeric(coef(direct)), as.numeric(coef(explicit)), tolerance = 1e-10)
  expect_equal(vcov(direct), vcov(explicit), tolerance = 1e-10)
})

test_that("npplreg residualized rank deficiency has a single informative condition", {
  if (!npplreg_audit_ready()) skip("Could not initialize npplreg audit context")

  z <- factor(rep(letters[1:4], each = 4L))
  x <- as.double(z)
  y <- 2.0 + 0.7 * x + rep(c(-0.03, 0.01, 0.02, 0.00), 4L)
  dat <- list(xdat = data.frame(x = x), zdat = data.frame(z = z), ydat = y)
  bw <- npplregbw(
    xdat = dat$xdat,
    zdat = dat$zdat,
    ydat = dat$ydat,
    bws = matrix(0, nrow = 2L, ncol = 1L),
    bandwidth.compute = FALSE
  )

  expect_error(
    npplreg(bws = bw),
    "residualized linear regressors are rank deficient after smoothing on z",
    fixed = TRUE
  )

  ns <- npplreg_audit_ns()
  common <- get(".np_plot_plreg_apply_common", envir = asNamespace(ns), inherits = FALSE)(
    bws = bw,
    txdat = dat$xdat,
    tzdat = dat$zdat
  )
  state <- get(".np_plot_plreg_apply_eval_state", envir = asNamespace(ns), inherits = FALSE)(
    common = common,
    exdat = dat$xdat,
    ezdat = dat$zdat
  )
  expect_error(
    get(".np_plot_plreg_apply_from_state", envir = asNamespace(ns), inherits = FALSE)(
      state = state,
      y = dat$ydat
    ),
    "residualized linear regressors are rank deficient after smoothing on z",
    fixed = TRUE
  )
})

test_that("npplreg evaluation rows retain alignment and external GOF is explicit", {
  if (!npplreg_audit_ready()) skip("Could not initialize npplreg audit context")

  dat <- npplreg_audit_fixture()
  bw <- npplreg_audit_bw(dat)
  ex <- data.frame(x = c(-0.8, NA, -0.1, 0.4, 0.8))
  ez <- data.frame(z = c(-0.9, -0.4, NA, 0.35, 0.9))
  ey <- 1.25 * ex$x + sin(2 * ez$z)

  fit <- npplreg(bws = bw, exdat = ex, ezdat = ez)
  expect_length(fit$mean, nrow(ex))
  expect_equal(fit$nobs, nrow(ex), tolerance = 0)
  expect_true(all(is.na(fit$mean[c(2L, 3L)])))
  expect_true(all(is.finite(fit$mean[-c(2L, 3L)])))
  expect_true(all(is.na(c(fit$R2, fit$MSE, fit$MAE, fit$MAPE, fit$CORR, fit$SIGN))))

  train <- npplreg(bws = bw)
  pred <- predict(train, exdat = ex, ezdat = ez, se.fit = TRUE)
  expect_length(pred$fit, nrow(ex))
  expect_length(pred$se.fit, nrow(ex))
  expect_true(all(is.na(pred$fit[c(2L, 3L)])))
  expect_true(all(is.na(pred$se.fit[c(2L, 3L)])))
  expect_true(all(is.finite(pred$fit[-c(2L, 3L)])))
  expect_true(all(is.finite(pred$se.fit[-c(2L, 3L)])))

  fit.ey <- npplreg(bws = bw, exdat = ex, ezdat = ez, eydat = ey)
  expect_length(fit.ey$mean, nrow(ex))
  expect_true(is.finite(fit.ey$R2))
  expect_true(is.finite(fit.ey$MSE))
})

test_that("npplreg factor response and factor linear regressor routes remain usable", {
  if (!npplreg_audit_ready()) skip("Could not initialize npplreg audit context")

  dat <- npplreg_audit_fixture()
  factor.y <- dat
  factor.y$ydat <- factor(ifelse(dat$ydat > median(dat$ydat), "hi", "lo"))
  fit.y <- npplreg(bws = npplreg_audit_bw(factor.y))
  expect_s3_class(fit.y, "plregression")
  expect_true(all(is.finite(fit.y$mean)))
  expect_true(all(is.finite(coef(fit.y))))

  factor.x <- dat
  factor.x$xdat <- data.frame(x = factor(ifelse(dat$xdat$x > median(dat$xdat$x), "high", "low")))
  fit.x <- npplreg(bws = npplreg_audit_bw(factor.x))
  expect_s3_class(fit.x, "plregression")
  expect_true(all(is.finite(fit.x$mean)))
  expect_true(all(is.finite(coef(fit.x))))
})

test_that("npplreg plot-center local fit stays parity-clean across bandwidth types", {
  if (!npplreg_audit_ready()) skip("Could not initialize npplreg audit context")

  dat <- npplreg_audit_fixture()
  ex <- data.frame(x = c(-0.8, -0.25, 0.4, 0.85))
  ez <- data.frame(z = c(-0.9, -0.2, 0.35, 0.95))
  local_fit <- get(".np_plot_plreg_local_fit",
                   envir = asNamespace(npplreg_audit_ns()),
                   inherits = FALSE)

  for (bwtype in c("fixed", "adaptive_nn", "generalized_nn")) {
    bw.value <- if (identical(bwtype, "fixed")) 0.55 else 5
    bw <- npplreg_audit_bw(dat, bw = bw.value, bwtype = bwtype)
    fit <- npplreg(bws = bw, exdat = ex, ezdat = ez)
    plotfit <- local_fit(
      bws = bw,
      xdat = dat$xdat,
      ydat = dat$ydat,
      zdat = dat$zdat,
      exdat = ex,
      ezdat = ez
    )

    expect_equal(as.numeric(plotfit$mean), as.numeric(fit$mean), tolerance = 1e-10)
    expect_equal(as.numeric(coef(plotfit)), as.numeric(coef(fit)), tolerance = 1e-10)
    expect_equal(vcov(plotfit), vcov(fit), tolerance = 1e-10)
  }
})
