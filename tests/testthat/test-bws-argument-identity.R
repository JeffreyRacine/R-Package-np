test_that("formula bandwidth selection does not partially match bwsolver", {
  if (!spawn_mpi_slaves()) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  set.seed(42)
  dat <- data.frame(x = runif(30, -1, 1))
  dat$y <- sin(dat$x) + rnorm(30, sd = .2)
  fit <- npRmpi::npreg(y ~ x, data = dat, regtype = "lp", degree = 2L,
    bwsolver = "mads+powell", nmulti = 1L,
    nomad.opts = list(MAX_BB_EVAL = 12L))
  bw <- npRmpi::npregbw(y ~ x, data = dat, regtype = "lp", degree = 2L,
    bwsolver = "mads+powell", nmulti = 1L,
    nomad.opts = list(MAX_BB_EVAL = 12L))
  expect_identical(fit$bws$bw, bw$bw)
  expect_identical(fitted(fit), fitted(npRmpi::npreg(bws = bw)))
})

test_that("IV derivative forwards its seed when solver is present without bws", {
  if (!spawn_mpi_slaves()) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old), add = TRUE)
  fun <- get("npregivderiv.default", asNamespace("npRmpi"))
  isolated <- new.env(parent = environment(fun))
  seen <- NULL
  isolated$.npreg_complete <- function(...) {
    seen <<- list(...)
    stop("captured initial IV stage", call. = FALSE)
  }
  environment(fun) <- isolated
  x <- seq(-1, 1, length.out = 20L)
  expect_error(fun(y = sin(x), z = x, w = x + .1,
    random.seed = 73L, bwsolver = "mads+powell", iterate.max = 2L,
    nmulti = 1L), "captured initial IV stage")
  expect_identical(seen[["random.seed", exact = TRUE]], 73L)
  expect_identical(seen[["bwsolver", exact = TRUE]], "mads+powell")
  expect_null(seen[["bws", exact = TRUE]])
})
