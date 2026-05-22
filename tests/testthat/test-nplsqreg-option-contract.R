test_that("nplsqreg honors estimator bandwidth options separately from pilot options", {
  options(np.messages = FALSE)
  set.seed(20260521)
  dat <- data.frame(
    y = sin(seq(0, 2 * pi, length.out = 30L)) + rnorm(30L, sd = 0.1),
    x = seq(0, 1, length.out = 30L),
    u = factor(rep(letters[1:3], length.out = 30L)),
    o = ordered(rep(1:4, length.out = 30L))
  )
  scale0 <- rep(1, nrow(dat))

  fit.ll <- nplsqreg(y ~ x, data = dat, regtype = "ll",
                     scale = scale0, nmulti = 1L,
                     optim.control = list(maxit = 2L))
  expect_identical(fit.ll$bws$reg.bws$regtype, "ll")

  fit.lp <- nplsqreg(y ~ x, data = dat, regtype = "lp", degree = 1L,
                     scale = scale0, nmulti = 1L,
                     optim.control = list(maxit = 2L))
  expect_identical(fit.lp$bws$reg.bws$regtype, "lp")
  expect_identical(as.integer(fit.lp$bws$reg.bws$degree), 1L)

  fit.nn <- nplsqreg(bws = 4L, txdat = dat["x"], tydat = dat$y,
                     bwtype = "generalized_nn", scale = scale0,
                     bandwidth.compute = FALSE)
  expect_identical(fit.nn$bws$reg.bws$type, "generalized_nn")

  fit.ker <- nplsqreg(
    y ~ x + u + o,
    data = dat,
    cxkertype = "epanechnikov",
    cxkerorder = 4L,
    uxkertype = "liracine",
    oxkertype = "wangvanryzin",
    scale = scale0,
    nmulti = 1L,
    optim.control = list(maxit = 2L)
  )
  expect_identical(fit.ker$bws$reg.bws$ckertype, "epanechnikov")
  expect_identical(as.integer(fit.ker$bws$reg.bws$ckerorder), 4L)
  expect_identical(fit.ker$bws$reg.bws$ukertype, "liracine")
  expect_identical(fit.ker$bws$reg.bws$okertype, "wangvanryzin")

  expect_error(
    nplsqreg(y ~ x, data = dat, cykertype = "epanechnikov",
             scale = scale0, nmulti = 1L,
             optim.control = list(maxit = 2L)),
    "response-side"
  )
  expect_error(
    nplsqreg(y ~ x, data = dat, total_nonsense = TRUE,
             scale = scale0, nmulti = 1L,
             optim.control = list(maxit = 2L)),
    "unused nplsqregbw argument"
  )
})

test_that("nplsqreg formula route honors native exdat precedence", {
  options(np.messages = FALSE)
  set.seed(20260521)
  dat <- data.frame(
    y = sin(seq(0, 2 * pi, length.out = 30L)) + rnorm(30L, sd = 0.1),
    x = seq(0, 1, length.out = 30L),
    u = factor(rep(letters[1:3], length.out = 30L))
  )
  scale0 <- rep(1, nrow(dat))
  newdata <- dat[seq_len(5L), c("x", "u"), drop = FALSE]
  exdat <- dat[seq_len(2L), c("x", "u"), drop = FALSE]

  fit <- nplsqreg(y ~ x + u, data = dat, newdata = newdata,
                  exdat = exdat, scale = scale0,
                  nmulti = 1L, optim.control = list(maxit = 2L))
  expect_identical(fit$nobs, 2L)
  expect_false(fit$trainiseval)
})

test_that("nplsqreg vector tau applies estimator options to each full-search tau", {
  options(np.messages = FALSE)
  set.seed(20260521)
  dat <- data.frame(
    y = cos(seq(0, 2 * pi, length.out = 24L)) + rnorm(24L, sd = 0.1),
    x = seq(0, 1, length.out = 24L)
  )

  fit <- nplsqreg(
    y ~ x,
    data = dat,
    tau = c(0.25, 0.5),
    tau.search = "full",
    regtype = "ll",
    scale = rep(1, nrow(dat)),
    nmulti = 1L,
    optim.control = list(maxit = 2L)
  )

  expect_identical(fit$bws$tau.search, "full")
  expect_identical(fit$bws$fit.order, seq_along(fit$tau))
  expect_true(all(vapply(
    fit$bws$tau.bws,
    function(bw) identical(bw$reg.bws$regtype, "ll"),
    logical(1)
  )))

  out.fit <- capture.output(summary(fit))
  out.bw <- capture.output(summary(fit$bws))
  expect_true(any(grepl("Location-Scale Quantile Regression Data", out.fit, fixed = TRUE)))
  expect_true(any(grepl("Exp. Var. Bandwidth", out.fit, fixed = TRUE)))
  expect_true(any(grepl("Search Method: powell", out.fit, fixed = TRUE)))
  expect_true(any(grepl("Location-Scale Quantile Regression Bandwidth Data", out.bw, fixed = TRUE)))
  expect_true(any(grepl("Exp. Var. Bandwidth", out.bw, fixed = TRUE)))
  expect_true(any(grepl("Bandwidth Selection Method: Check-Loss Cross-Validation", out.bw, fixed = TRUE)))
})

test_that("nplsqreg plot uses quantile plot contract for multiple slices", {
  options(np.messages = FALSE)
  set.seed(20260522)
  dat <- data.frame(
    y = sin(seq(0, 2, length.out = 32L)) + rnorm(32L, sd = 0.1),
    x = seq(0, 1, length.out = 32L),
    z = ordered(rep(1:4, length.out = 32L))
  )
  scale0 <- rep(1, nrow(dat))

  fit <- nplsqreg(y ~ x + z, data = dat, tau = c(0.25, 0.5),
                  tau.search = "full", scale = scale0, nmulti = 1L,
                  optim.control = list(maxit = 2L))

  out <- plot(fit, output = "data", perspective = FALSE, neval = 6L)
  expect_identical(names(out), c("cd1", "cd2"))
  expect_s3_class(out$cd1, "lsqregression")
  expect_s3_class(out$cd2, "lsqregression")
  expect_true(is.matrix(out$cd1$quantile))
  expect_identical(colnames(out$cd1$quantile), c("tau=0.25", "tau=0.50"))

  gout <- plot(fit, output = "data", perspective = FALSE, neval = 6L,
               gradients = TRUE)
  expect_identical(names(gout), c("cd1", "cd2"))
  expect_true(length(dim(gout$cd1$quantgrad)) == 3L)

  aout <- plot(fit, output = "data", perspective = FALSE, neval = 6L,
               errors = "asymptotic")
  expect_true(all(c("quanterr", "bias") %in% names(aout$cd1)))

  bout <- plot(fit, output = "data", perspective = FALSE, neval = 5L,
               errors = "bootstrap", B = 3L)
  expect_true(all(c("quanterr", "bias") %in% names(bout$cd1)))

  dat2 <- data.frame(
    y = dat$y,
    x1 = seq(0, 1, length.out = nrow(dat)),
    x2 = cos(seq(0, 2, length.out = nrow(dat)))
  )
  fit2 <- nplsqreg(y ~ x1 + x2, data = dat2, tau = 0.5,
                   scale = scale0, nmulti = 1L,
                   optim.control = list(maxit = 2L))
  out2 <- plot(fit2, output = "data", perspective = FALSE, neval = 6L)
  expect_identical(names(out2), c("cd1", "cd2"))
  expect_s3_class(out2$cd1, "lsqregression")
  expect_s3_class(out2$cd2, "lsqregression")

  expect_error(
    plot(fit$tau.fits[[1L]], tau = c(0.25, 0.5), output = "data",
         perspective = FALSE),
    "must have been fitted"
  )
})
