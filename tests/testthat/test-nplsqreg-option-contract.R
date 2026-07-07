test_that("session nplsqreg honors estimator bandwidth options separately from pilot options", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "options(np.messages = FALSE)",
      "set.seed(20260521)",
      "dat <- data.frame(y = sin(seq(0, 2 * pi, length.out = 30L)) + rnorm(30L, sd = 0.1), x = seq(0, 1, length.out = 30L), u = factor(rep(letters[1:3], length.out = 30L)), o = ordered(rep(1:4, length.out = 30L)))",
      "scale0 <- rep(1, nrow(dat))",
      "fit.ll <- nplsqreg(y ~ x, data = dat, regtype = 'll', scale = scale0, nmulti = 1L, optim.control = list(maxit = 2L))",
      "stopifnot(identical(fit.ll$bws$reg.bws$regtype, 'll'))",
      "fit.lp <- nplsqreg(y ~ x, data = dat, regtype = 'lp', degree = 1L, scale = scale0, nmulti = 1L, optim.control = list(maxit = 2L))",
      "stopifnot(identical(fit.lp$bws$reg.bws$regtype, 'lp'))",
      "stopifnot(identical(as.integer(fit.lp$bws$reg.bws$degree), 1L))",
      "fit.nn <- nplsqreg(bws = 4L, txdat = dat['x'], tydat = dat$y, bwtype = 'generalized_nn', scale = scale0, bandwidth.compute = FALSE)",
      "stopifnot(identical(fit.nn$bws$reg.bws$type, 'generalized_nn'))",
      "fit.ker <- nplsqreg(y ~ x + u + o, data = dat, cxkertype = 'epanechnikov', cxkerorder = 4L, uxkertype = 'liracine', oxkertype = 'wangvanryzin', scale = scale0, nmulti = 1L, optim.control = list(maxit = 2L))",
      "stopifnot(identical(fit.ker$bws$reg.bws$ckertype, 'epanechnikov'))",
      "stopifnot(identical(as.integer(fit.ker$bws$reg.bws$ckerorder), 4L))",
      "stopifnot(identical(fit.ker$bws$reg.bws$ukertype, 'liracine'))",
      "stopifnot(identical(fit.ker$bws$reg.bws$okertype, 'wangvanryzin'))",
      "stopifnot(inherits(try(nplsqreg(y ~ x, data = dat, cykertype = 'epanechnikov', scale = scale0, nmulti = 1L, optim.control = list(maxit = 2L)), silent = TRUE), 'try-error'))",
      "stopifnot(inherits(try(nplsqreg(y ~ x, data = dat, total_nonsense = TRUE, scale = scale0, nmulti = 1L, optim.control = list(maxit = 2L)), silent = TRUE), 'try-error'))",
      "nd <- dat[seq_len(5L), c('x', 'u'), drop = FALSE]",
      "ed <- dat[seq_len(2L), c('x', 'u'), drop = FALSE]",
      "fit.eval <- nplsqreg(y ~ x + u, data = dat, newdata = nd, exdat = ed, scale = scale0, nmulti = 1L, optim.control = list(maxit = 2L))",
      "stopifnot(identical(fit.eval$nobs, 2L))",
      "stopifnot(!isTRUE(fit.eval$trainiseval))",
      "fit.pred <- nplsqreg(y ~ x + u, data = dat, scale = scale0, nmulti = 1L, optim.control = list(maxit = 2L))",
      "p.new <- predict(fit.pred, newdata = nd)",
      "ref.new <- fitted(nplsqreg(bws = fit.pred$bws, exdat = nd))",
      "stopifnot(isTRUE(all.equal(p.new, ref.new)))",
      "p.both <- predict(fit.pred, newdata = nd, exdat = ed)",
      "ref.ed <- predict(fit.pred, exdat = ed)",
      "stopifnot(identical(length(p.both), nrow(ed)))",
      "stopifnot(isTRUE(all.equal(p.both, ref.ed)))",
      "p.se <- predict(fit.pred, newdata = nd, se.fit = TRUE)",
      "stopifnot(all(c('fit', 'se.fit', 'df', 'residual.scale') %in% names(p.se)))",
      "stopifnot(identical(length(p.se$fit), nrow(nd)))",
      "stopifnot(identical(length(p.se$se.fit), nrow(nd)))",
      "fit.vec <- nplsqreg(y ~ x, data = dat, tau = c(0.25, 0.5), tau.search = 'full', regtype = 'll', scale = scale0, nmulti = 1L, optim.control = list(maxit = 2L))",
      "pv <- predict(fit.vec, newdata = dat[seq_len(4L), 'x', drop = FALSE], se.fit = TRUE)",
      "stopifnot(is.matrix(pv$fit))",
      "stopifnot(is.matrix(pv$se.fit))",
      "stopifnot(identical(dim(pv$fit), c(4L, length(fit.vec$tau))))",
      "stopifnot(identical(colnames(pv$fit), c('tau=0.25', 'tau=0.50')))",
      "stopifnot(identical(colnames(pv$se.fit), c('tau=0.25', 'tau=0.50')))",
      "stopifnot(identical(fit.vec$bws$tau.search, 'full'))",
      "stopifnot(identical(fit.vec$bws$fit.order, seq_along(fit.vec$tau)))",
      "stopifnot(all(vapply(fit.vec$bws$tau.bws, function(bw) identical(bw$reg.bws$regtype, 'll'), logical(1))))",
      "out.fit <- capture.output(summary(fit.vec))",
      "out.bw <- capture.output(summary(fit.vec$bws))",
      "stopifnot(any(grepl('Location-Scale Quantile Regression Data', out.fit, fixed = TRUE)))",
      "stopifnot(any(grepl('Exp. Var. Bandwidth', out.fit, fixed = TRUE)))",
      "stopifnot(any(grepl('Search Method: powell', out.fit, fixed = TRUE)))",
      "stopifnot(any(grepl('Location-Scale Quantile Regression Bandwidth Data', out.bw, fixed = TRUE)))",
      "stopifnot(any(grepl('Exp. Var. Bandwidth', out.bw, fixed = TRUE)))",
      "stopifnot(any(grepl('Bandwidth Selection Method: Check-Loss Cross-Validation', out.bw, fixed = TRUE)))",
      "fit.no.resid <- nplsqreg(y ~ x, data = dat, tau = 0.5, scale = scale0, nmulti = 1L, optim.control = list(maxit = 2L))",
      "stopifnot(inherits(try(residuals(fit.no.resid), silent = TRUE), 'try-error'))",
      "fit.scalar.resid <- nplsqreg(y ~ x, data = dat, tau = 0.5, scale = scale0, residuals = TRUE, nmulti = 1L, optim.control = list(maxit = 2L))",
      "stopifnot(identical(length(residuals(fit.scalar.resid)), nrow(dat)))",
      "stopifnot(isTRUE(all.equal(residuals(fit.scalar.resid), fit.scalar.resid$resid)))",
      "stopifnot(isTRUE(all.equal(residuals(fit.scalar.resid), dat$y - fitted(fit.scalar.resid))))",
      "eval.dat <- dat[seq_len(5L), 'x', drop = FALSE]",
      "fit.eval.resid <- nplsqreg(y ~ x, data = dat, newdata = eval.dat, tau = 0.5, scale = scale0, residuals = TRUE, nmulti = 1L, optim.control = list(maxit = 2L))",
      "stopifnot(identical(length(fitted(fit.eval.resid)), nrow(eval.dat)))",
      "stopifnot(identical(length(residuals(fit.eval.resid)), nrow(dat)))",
      "fit.vector.resid <- nplsqreg(y ~ x, data = dat, tau = c(0.25, 0.5), scale = scale0, residuals = TRUE, nmulti = 1L, optim.control = list(maxit = 2L))",
      "rv <- residuals(fit.vector.resid)",
      "stopifnot(is.matrix(rv))",
      "stopifnot(identical(dim(rv), c(nrow(dat), length(fit.vector.resid$tau))))",
      "stopifnot(identical(colnames(rv), c('tau=0.25', 'tau=0.50')))",
      "stopifnot(isTRUE(all.equal(rv[, 1L], fit.vector.resid$tau.fits[[1L]]$resid)))",
      "stopifnot(isTRUE(all.equal(rv[, 2L], fit.vector.resid$tau.fits[[2L]]$resid)))",
      "dat.plot <- data.frame(y = sin(seq(0, 2, length.out = 32L)) + rnorm(32L, sd = 0.1), x = seq(0, 1, length.out = 32L), z = ordered(rep(1:4, length.out = 32L)))",
      "fit.plot <- nplsqreg(y ~ x + z, data = dat.plot, tau = c(0.25, 0.5), tau.search = 'full', scale = rep(1, nrow(dat.plot)), nmulti = 1L, optim.control = list(maxit = 2L))",
      "pout <- plot(fit.plot, output = 'data', perspective = FALSE, neval = 6L)",
      "stopifnot(identical(names(pout), c('cd1', 'cd2')))",
      "stopifnot(inherits(pout$cd1, 'lsqregression'))",
      "stopifnot(is.matrix(pout$cd1$quantile))",
      "stopifnot(identical(colnames(pout$cd1$quantile), c('tau=0.25', 'tau=0.50')))",
      "gout <- plot(fit.plot, output = 'data', perspective = FALSE, neval = 6L, gradients = TRUE)",
      "stopifnot(identical(names(gout), c('cd1', 'cd2')))",
      "dat.lp <- data.frame(x = seq(0, 1, length.out = 40L))",
      "dat.lp$y <- dat.lp$x + dat.lp$x^2 + rnorm(nrow(dat.lp), sd = 0.01)",
      "fit.lp <- nplsqreg(y ~ x, data = dat.lp, tau = c(0.25, 0.5), scale = rep(1, nrow(dat.lp)), regtype = 'lp', degree = 1L, nmulti = 1L, optim.control = list(maxit = 2L))",
      "q.lp <- plot(fit.lp, output = 'data', perspective = FALSE, neval = 41L)$cd1",
      "g.lp <- plot(fit.lp, output = 'data', perspective = FALSE, neval = 41L, gradients = TRUE)$cd1",
      "x.eval <- q.lp$xeval[[1L]]",
      "q.mat <- fitted(q.lp)",
      "g.arr <- gradients(g.lp)",
      "fd.center <- function(x, y) { n <- length(y); out <- numeric(n); out[1L] <- (y[2L] - y[1L]) / (x[2L] - x[1L]); out[n] <- (y[n] - y[n - 1L]) / (x[n] - x[n - 1L]); out[2:(n - 1L)] <- (y[3:n] - y[1:(n - 2L)]) / (x[3:n] - x[1:(n - 2L)]); out }",
      "for (j in seq_len(ncol(q.mat))) stopifnot(max(abs(g.arr[, 1L, j] - fd.center(x.eval, q.mat[, j]))) < 0.05)",
      "aout <- plot(fit.plot, output = 'data', perspective = FALSE, neval = 6L, errors = 'asymptotic')",
      "stopifnot(all(c('quanterr', 'bias') %in% names(aout$cd1)))",
      "bout <- plot(fit.plot, output = 'data', perspective = FALSE, neval = 5L, errors = 'bootstrap', B = 3L)",
      "stopifnot(all(c('quanterr', 'bias') %in% names(bout$cd1)))",
      "dat2 <- data.frame(y = dat.plot$y, x1 = seq(0, 1, length.out = nrow(dat.plot)), x2 = cos(seq(0, 2, length.out = nrow(dat.plot))))",
      "fit2 <- nplsqreg(y ~ x1 + x2, data = dat2, tau = 0.5, scale = rep(1, nrow(dat2)), nmulti = 1L, optim.control = list(maxit = 2L))",
      "out2 <- plot(fit2, output = 'data', perspective = FALSE, neval = 6L)",
      "stopifnot(identical(names(out2), c('cd1', 'cd2')))",
      "stopifnot(inherits(out2$cd1, 'lsqregression'))",
      "stopifnot(inherits(out2$cd2, 'lsqregression'))",
      "stopifnot(inherits(try(plot(fit.plot$tau.fits[[1L]], tau = c(0.25, 0.5), output = 'data', perspective = FALSE), silent = TRUE), 'try-error'))",
      "cat('NPLSQREG_OPTION_CONTRACT_OK\\n')"
    ),
    timeout = 90L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NPLSQREG_OPTION_CONTRACT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session nplsqreg formula subset and na.action match explicit fit path", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "options(np.messages = FALSE)",
      "set.seed(20260523)",
      "dat <- data.frame(y = sin(seq(0, 2 * pi, length.out = 30L)) + rnorm(30L, sd = 0.1), x = seq(0, 1, length.out = 30L))",
      "dat$y[c(4L, 17L)] <- NA_real_",
      "dat$x[9L] <- NA_real_",
      "keep <- dat$x > 0.15",
      "mf <- model.frame(y ~ x, data = dat, subset = keep, na.action = na.omit)",
      "xdat <- mf[, 'x', drop = FALSE]",
      "ydat <- model.response(mf)",
      "scale.sub <- rep(1, nrow(mf))",
      "fit.form <- nplsqreg(y ~ x, data = dat, subset = keep, na.action = na.omit, scale = scale.sub, nmulti = 1L, optim.control = list(maxit = 2L))",
      "fit.exp <- nplsqreg(bws = fit.form$bws, txdat = xdat, tydat = ydat)",
      "stopifnot(identical(fit.form$ntrain, nrow(mf)))",
      "stopifnot(identical(fit.form$xnames, 'x'))",
      "stopifnot(identical(fit.form$ynames, 'y'))",
      "stopifnot(identical(length(fitted(fit.form)), fit.form$nobs))",
      "stopifnot(identical(length(fitted(fit.exp)), nrow(mf)))",
      "stopifnot(identical(as.integer(which(is.na(fitted(fit.form)))), as.integer(fit.form$omit)))",
      "stopifnot(isTRUE(all.equal(as.numeric(stats::na.omit(fitted(fit.form))), as.numeric(fitted(fit.exp)))))",
      "stopifnot(inherits(fit.form$bws$formula, 'formula'))",
      "cat('NPLSQREG_FORMULA_SUBSET_NA_OK\\n')"
    ),
    timeout = 90L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NPLSQREG_FORMULA_SUBSET_NA_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("session nplsqreg NOMAD vector tau summary reports per-tau degrees", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "options(np.messages = FALSE)",
      "set.seed(20260524)",
      "dat <- data.frame(y = sin(seq(0, 2 * pi, length.out = 24L)) + rnorm(24L, sd = 0.1), x1 = seq(0, 1, length.out = 24L), x2 = cos(seq(0, 2, length.out = 24L)))",
      "fit <- nplsqreg(y ~ x1 + x2, data = dat, tau = c(0.25, 0.5), tau.search = 'full', nomad = TRUE, scale = rep(1, nrow(dat)), nmulti = 1L, nomad.nmulti = 0L, degree.max = 1L, degree.restarts = 0L, optim.control = list(maxit = 1L))",
      "out.fit <- capture.output(summary(fit))",
      "out.bw <- capture.output(summary(fit$bws))",
      "stopifnot(any(grepl('Continuous LP Degree(s):', out.fit, fixed = TRUE)))",
      "stopifnot(any(grepl('Continuous LP Degree(s):', out.bw, fixed = TRUE)))",
      "stopifnot(any(grepl('tau=0.25', out.fit, fixed = TRUE)))",
      "stopifnot(any(grepl('tau=0.50', out.fit, fixed = TRUE)))",
      "stopifnot(any(grepl('x1', out.fit, fixed = TRUE)))",
      "stopifnot(any(grepl('x2', out.fit, fixed = TRUE)))",
      "cat('NPLSQREG_NOMAD_DEGREE_SUMMARY_OK\\n')"
    ),
    timeout = 90L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NPLSQREG_NOMAD_DEGREE_SUMMARY_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
