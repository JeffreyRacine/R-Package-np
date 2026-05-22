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
      "fit.vec <- nplsqreg(y ~ x, data = dat, tau = c(0.25, 0.5), tau.search = 'full', regtype = 'll', scale = scale0, nmulti = 1L, optim.control = list(maxit = 2L))",
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
      "dat.plot <- data.frame(y = sin(seq(0, 2, length.out = 32L)) + rnorm(32L, sd = 0.1), x = seq(0, 1, length.out = 32L), z = ordered(rep(1:4, length.out = 32L)))",
      "fit.plot <- nplsqreg(y ~ x + z, data = dat.plot, tau = c(0.25, 0.5), tau.search = 'full', scale = rep(1, nrow(dat.plot)), nmulti = 1L, optim.control = list(maxit = 2L))",
      "pout <- plot(fit.plot, output = 'data', perspective = FALSE, neval = 6L)",
      "stopifnot(identical(names(pout), c('cd1', 'cd2')))",
      "stopifnot(inherits(pout$cd1, 'lsqregression'))",
      "stopifnot(is.matrix(pout$cd1$quantile))",
      "stopifnot(identical(colnames(pout$cd1$quantile), c('tau=0.25', 'tau=0.50')))",
      "gout <- plot(fit.plot, output = 'data', perspective = FALSE, neval = 6L, gradients = TRUE)",
      "stopifnot(identical(names(gout), c('cd1', 'cd2')))",
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
