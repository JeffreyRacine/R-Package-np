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
      "fit.ker <- nplsqreg(y ~ x + u + o, data = dat, ckertype = 'epanechnikov', ckerorder = 4L, ukertype = 'liracine', okertype = 'wangvanryzin', scale = scale0, nmulti = 1L, optim.control = list(maxit = 2L))",
      "stopifnot(identical(fit.ker$bws$reg.bws$ckertype, 'epanechnikov'))",
      "stopifnot(identical(as.integer(fit.ker$bws$reg.bws$ckerorder), 4L))",
      "stopifnot(identical(fit.ker$bws$reg.bws$ukertype, 'liracine'))",
      "stopifnot(identical(fit.ker$bws$reg.bws$okertype, 'wangvanryzin'))",
      "fit.vec <- nplsqreg(y ~ x, data = dat, tau = c(0.25, 0.5), tau.search = 'full', regtype = 'll', scale = scale0, nmulti = 1L, optim.control = list(maxit = 2L))",
      "stopifnot(identical(fit.vec$bws$tau.search, 'full'))",
      "stopifnot(identical(fit.vec$bws$fit.order, seq_along(fit.vec$tau)))",
      "stopifnot(all(vapply(fit.vec$bws$tau.bws, function(bw) identical(bw$reg.bws$regtype, 'll'), logical(1))))",
      "cat('NPLSQREG_OPTION_CONTRACT_OK\\n')"
    ),
    timeout = 90L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NPLSQREG_OPTION_CONTRACT_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
