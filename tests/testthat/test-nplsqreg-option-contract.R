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

test_that("predict.lsqregression handles newdata, exdat, se.fit, and vector tau", {
  options(np.messages = FALSE)
  set.seed(20260523)
  dat <- data.frame(
    y = sin(seq(0, 2 * pi, length.out = 30L)) + rnorm(30L, sd = 0.1),
    x = seq(0, 1, length.out = 30L),
    u = factor(rep(letters[1:3], length.out = 30L))
  )
  scale0 <- rep(1, nrow(dat))
  newdata <- dat[seq_len(5L), c("x", "u"), drop = FALSE]
  exdat <- dat[seq_len(2L), c("x", "u"), drop = FALSE]

  fit <- nplsqreg(y ~ x + u, data = dat, scale = scale0,
                  nmulti = 1L, optim.control = list(maxit = 2L))
  p.new <- predict(fit, newdata = newdata)
  ref.new <- fitted(nplsqreg(bws = fit$bws, exdat = newdata))
  expect_equal(p.new, ref.new)

  p.both <- predict(fit, newdata = newdata, exdat = exdat)
  ref.exdat <- predict(fit, exdat = exdat)
  expect_length(p.both, nrow(exdat))
  expect_equal(p.both, ref.exdat)

  p.se <- predict(fit, newdata = newdata, se.fit = TRUE)
  expect_true(all(c("fit", "se.fit", "df", "residual.scale") %in% names(p.se)))
  expect_length(p.se$fit, nrow(newdata))
  expect_length(p.se$se.fit, nrow(newdata))

  fit.vec <- nplsqreg(y ~ x, data = dat, tau = c(0.25, 0.5),
                      scale = scale0, nmulti = 1L,
                      optim.control = list(maxit = 2L))
  nd.x <- dat[seq_len(4L), "x", drop = FALSE]
  pv <- predict(fit.vec, newdata = nd.x, se.fit = TRUE)
  expect_true(is.matrix(pv$fit))
  expect_true(is.matrix(pv$se.fit))
  expect_identical(dim(pv$fit), c(nrow(nd.x), length(fit.vec$tau)))
  expect_identical(colnames(pv$fit), c("tau=0.25", "tau=0.50"))
  expect_identical(colnames(pv$se.fit), c("tau=0.25", "tau=0.50"))
})

test_that("nplsqreg formula subset and na.action match explicit mixed-data path", {
  options(np.messages = FALSE)
  set.seed(20260523)
  dat <- data.frame(
    y = sin(seq(0, 2 * pi, length.out = 36L)) + rnorm(36L, sd = 0.1),
    x = seq(0, 1, length.out = 36L),
    o = ordered(rep(1:4, length.out = 36L))
  )
  dat$y[c(4L, 21L)] <- NA_real_
  dat$x[9L] <- NA_real_
  keep <- dat$x > 0.15
  mf <- model.frame(y ~ x + o, data = dat, subset = keep,
                    na.action = na.omit)
  xdat <- mf[, c("x", "o"), drop = FALSE]
  ydat <- model.response(mf)
  scale.sub <- rep(1, nrow(mf))

  bw.form <- nplsqregbw(y ~ x + o, data = dat, subset = keep,
                        na.action = na.omit, scale = scale.sub,
                        nmulti = 1L, optim.control = list(maxit = 2L))
  bw.exp <- nplsqregbw(xdat = xdat, ydat = ydat, scale = scale.sub,
                       nmulti = 1L, optim.control = list(maxit = 2L))
  expect_equal(bw.form$bw, bw.exp$bw)
  expect_identical(bw.form$nobs, nrow(mf))
  expect_identical(bw.form$xnames, c("x", "o"))
  expect_identical(bw.form$ynames, "y")

  fit.form <- nplsqreg(y ~ x + o, data = dat, subset = keep,
                       na.action = na.omit, scale = scale.sub,
                       nmulti = 1L, optim.control = list(maxit = 2L))
  expect_identical(fit.form$ntrain, nrow(mf))
  expect_identical(fit.form$xnames, c("x", "o"))
  expect_identical(fit.form$ynames, "y")
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

test_that("nplsqreg NOMAD vector tau summary reports per-tau degrees", {
  options(np.messages = FALSE)
  set.seed(20260524)
  dat <- data.frame(
    y = sin(seq(0, 2 * pi, length.out = 24L)) + rnorm(24L, sd = 0.1),
    x1 = seq(0, 1, length.out = 24L),
    x2 = cos(seq(0, 2, length.out = 24L))
  )

  fit <- nplsqreg(
    y ~ x1 + x2,
    data = dat,
    tau = c(0.25, 0.5),
    tau.search = "full",
    nomad = TRUE,
    scale = rep(1, nrow(dat)),
    nmulti = 1L,
    nomad.nmulti = 0L,
    degree.max = 1L,
    degree.restarts = 0L,
    optim.control = list(maxit = 1L)
  )

  out.fit <- capture.output(summary(fit))
  out.bw <- capture.output(summary(fit$bws))
  expect_true(any(grepl("Continuous LP Degree(s):", out.fit, fixed = TRUE)))
  expect_true(any(grepl("Continuous LP Degree(s):", out.bw, fixed = TRUE)))
  expect_true(any(grepl("tau=0.25", out.fit, fixed = TRUE)))
  expect_true(any(grepl("tau=0.50", out.fit, fixed = TRUE)))
  expect_true(any(grepl("x1", out.fit, fixed = TRUE)))
  expect_true(any(grepl("x2", out.fit, fixed = TRUE)))
})

test_that("nplsqreg residuals accessor exposes requested residuals", {
  options(np.messages = FALSE)
  set.seed(20260522)
  dat <- data.frame(
    y = cos(seq(0, 2 * pi, length.out = 24L)) + rnorm(24L, sd = 0.1),
    x = seq(0, 1, length.out = 24L)
  )
  scale0 <- rep(1, nrow(dat))

  fit.no <- nplsqreg(y ~ x, data = dat, tau = 0.5, scale = scale0,
                     nmulti = 1L, optim.control = list(maxit = 2L))
  expect_error(residuals(fit.no), "refit with residuals=TRUE")

  fit.scalar <- nplsqreg(y ~ x, data = dat, tau = 0.5, scale = scale0,
                         residuals = TRUE, nmulti = 1L,
                         optim.control = list(maxit = 2L))
  expect_length(residuals(fit.scalar), nrow(dat))
  expect_equal(residuals(fit.scalar), fit.scalar$resid)
  expect_equal(residuals(fit.scalar), dat$y - fitted(fit.scalar))

  eval.dat <- dat[seq_len(5L), "x", drop = FALSE]
  fit.eval <- nplsqreg(y ~ x, data = dat, newdata = eval.dat, tau = 0.5,
                       scale = scale0, residuals = TRUE, nmulti = 1L,
                       optim.control = list(maxit = 2L))
  expect_length(fitted(fit.eval), nrow(eval.dat))
  expect_length(residuals(fit.eval), nrow(dat))

  fit.vector <- nplsqreg(y ~ x, data = dat, tau = c(0.25, 0.5),
                         scale = scale0, residuals = TRUE, nmulti = 1L,
                         optim.control = list(maxit = 2L))
  r <- residuals(fit.vector)
  expect_true(is.matrix(r))
  expect_identical(dim(r), c(nrow(dat), length(fit.vector$tau)))
  expect_identical(colnames(r), c("tau=0.25", "tau=0.50"))
  expect_equal(r[, 1L], fit.vector$tau.fits[[1L]]$resid)
  expect_equal(r[, 2L], fit.vector$tau.fits[[2L]]$resid)
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

  dat.lp <- data.frame(x = seq(0, 1, length.out = 40L))
  dat.lp$y <- dat.lp$x + dat.lp$x^2 + rnorm(nrow(dat.lp), sd = 0.01)
  fit.lp <- nplsqreg(y ~ x, data = dat.lp, tau = c(0.25, 0.5),
                     scale = rep(1, nrow(dat.lp)), regtype = "lp",
                     degree = 1L, nmulti = 1L,
                     optim.control = list(maxit = 2L))
  q.lp <- plot(fit.lp, output = "data", perspective = FALSE, neval = 41L)$cd1
  g.lp <- plot(fit.lp, output = "data", perspective = FALSE, neval = 41L,
               gradients = TRUE)$cd1
  x.eval <- q.lp$xeval[[1L]]
  q.mat <- fitted(q.lp)
  g.arr <- gradients(g.lp)
  fd.center <- function(x, y) {
    n <- length(y)
    out <- numeric(n)
    out[1L] <- (y[2L] - y[1L]) / (x[2L] - x[1L])
    out[n] <- (y[n] - y[n - 1L]) / (x[n] - x[n - 1L])
    out[2:(n - 1L)] <- (y[3:n] - y[1:(n - 2L)]) / (x[3:n] - x[1:(n - 2L)])
    out
  }
  for (j in seq_len(ncol(q.mat))) {
    expect_lt(max(abs(g.arr[, 1L, j] - fd.center(x.eval, q.mat[, j]))), 0.05)
  }

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
