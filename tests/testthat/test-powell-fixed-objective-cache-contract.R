test_that("np.objective.cache preserves fixed Powell regression parity under MPI", {
  env <- npRmpi_subprocess_env("NP_RMPI_NO_REUSE_SLAVES=1")
  skip_if(is.null(env))

  out <- npRmpi_run_rscript_subprocess(
    c(
      "library(npRmpi)",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "fixed_cache_regression_data <- function(kind, n = 80L) {",
      "  set.seed(101)",
      "  if (identical(kind, 'continuous')) {",
      "    dat <- data.frame(x1 = runif(n), x2 = runif(n))",
      "    dat$y <- sin(dat$x1) + dat$x2 + rnorm(n, sd = 0.2)",
      "    return(dat)",
      "  }",
      "  if (identical(kind, 'mixed')) {",
      "    dat <- data.frame(x1 = runif(n), xo = ordered(sample(1:3, n, TRUE)))",
      "    dat$y <- dat$x1 + as.integer(dat$xo) / 3 + rnorm(n, sd = 0.2)",
      "    return(dat)",
      "  }",
      "  dat <- data.frame(xo = ordered(sample(1:4, n, TRUE)), xu = factor(sample(letters[1:3], n, TRUE)))",
      "  dat$y <- as.integer(dat$xo) + as.integer(dat$xu) + rnorm(n, sd = 0.2)",
      "  dat",
      "}",
      "fixed_cache_formula <- function(kind) {",
      "  if (identical(kind, 'continuous')) return(y ~ x1 + x2)",
      "  if (identical(kind, 'mixed')) return(y ~ x1 + xo)",
      "  y ~ xo + xu",
      "}",
      "run_fixed_cache_bw <- function(kind, cache) {",
      "  old <- options(np.messages = FALSE, np.tree = FALSE, np.objective.cache = cache)",
      "  on.exit(options(old), add = TRUE)",
      "  dat <- fixed_cache_regression_data(kind)",
      "  npregbw(fixed_cache_formula(kind), data = dat, bwmethod = 'cv.ls', bwtype = 'fixed', regtype = 'lc', nmulti = 1)",
      "}",
      "for (kind in c('continuous', 'mixed', 'categorical')) {",
      "  cached <- run_fixed_cache_bw(kind, TRUE)",
      "  uncached <- run_fixed_cache_bw(kind, FALSE)",
      "  stopifnot(isTRUE(all.equal(cached$bw, uncached$bw, tolerance = 0)))",
      "  stopifnot(isTRUE(all.equal(cached$fval, uncached$fval, tolerance = 1e-12)))",
      "  stopifnot(identical(as.numeric(cached$num.feval[1L]), as.numeric(uncached$num.feval[1L])))",
      "  stopifnot(as.numeric(cached$num.feval.fast[1L]) > 0)",
      "  stopifnot(identical(as.numeric(uncached$num.feval.fast[1L]), 0))",
      "  stopifnot(identical(unname(cached$nn.cache[['objective.enabled']]), 1))",
      "  stopifnot(unname(cached$nn.cache[['objective.hits']]) > 0)",
      "  stopifnot(identical(as.numeric(cached$num.feval.fast[1L]), unname(cached$nn.cache[['objective.hits']])))",
      "  stopifnot(identical(unname(uncached$nn.cache[['objective.enabled']]), 0))",
      "  stopifnot(identical(unname(uncached$nn.cache[['objective.hits']]), 0))",
      "}",
      "npRmpi.quit()"
    ),
    timeout = 90L,
    env = env
  )

  expect_equal(out$status, 0L, info = paste(out$output, collapse = "\n"))
})
