library(npRmpi)

test_that("provided fixed bounded cv.ls eval_only survives installed subprocess teardown", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "set.seed(42)",
      "n <- 80L",
      "x <- runif(n)",
      "y <- rbeta(n, 1, 1)",
      "tx <- data.frame(x = x)",
      "ty <- data.frame(y = y)",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "bw.ls <- npcdensbw(",
      "  xdat = tx,",
      "  ydat = ty,",
      "  bws = c(0.16, 0.14),",
      "  regtype = 'lp',",
      "  degree = 3L,",
      "  bwtype = 'fixed',",
      "  bwmethod = 'cv.ls',",
      "  cxkerbound = 'range',",
      "  cykerbound = 'range',",
      "  bandwidth.compute = FALSE",
      ")",
      "npRmpi.quit(force = TRUE)",
      "out <- npRmpi:::.npcdensbw_eval_only(tx, ty, bw.ls)",
      "stopifnot(is.list(out), is.finite(out$objective))",
      "cat('NPCDENS_CVLS_FIXED_BW_EVAL_ONLY_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NPCDENS_CVLS_FIXED_BW_EVAL_ONLY_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("public npcdensbw cv.ls lc matches the production fixed-point objective", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(222)
  n <- 32L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = rnorm(n))

  bw.lc <- npcdensbw(xdat = x, ydat = y, regtype = "lc", bwmethod = "cv.ls", nmulti = 1)

  expect_equal(
    bw.lc$fval,
    npRmpi:::.npcdensbw_eval_only(x, y$y1, bw.lc)$objective,
    tolerance = 1e-10
  )
})

test_that("public npcdensbw cv.ls lc frozen benchmark stays on the serial-aligned route", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  n <- 250L
  rho <- 0.25
  mu <- c(0, 0)
  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  data <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  mydat <- data.frame(x = data[, 2], y = data[, 1])

  bw.lc <- npcdensbw(y ~ x, data = mydat, bwmethod = "cv.ls")

  expect_equal(bw.lc$fval, 0.294769693962935, tolerance = 1e-12)
})

test_that("provided fixed lc cv.ls eval_only remains finite", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(42)
  x <- data.frame(x = runif(80L))
  y <- data.frame(y = rbeta(80L, 1, 1))

  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lc",
    bwtype = "fixed",
    bwmethod = "cv.ls",
    bws = c(0.4, 0.4),
    bandwidth.compute = FALSE
  )

  out <- npRmpi:::.npcdensbw_eval_only(x, y, bw)

  expect_true(is.list(out))
  expect_true(is.finite(out$objective))
  expect_equal(out$num.feval, 1)
})

test_that("public npcdensbw cv.ls fixed LP/LL route activates with ll == lp parity", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(141)
  n <- 36L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1 + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwmethod = "cv.ls",
    nmulti = 1
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ls",
    nmulti = 1
  )

  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})

test_that("public npcdensbw cv.ls fixed LP tree and serial evaluators agree at fixed points", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(142)
  n <- 34L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = sin(2 * pi * x$x1) + x$x2 + rnorm(n, sd = 0.12))
  degree <- rep.int(1L, ncol(x))

  bw.serial <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ls",
    nmulti = 1
  )

  old_opt <- getOption("np.tree")
  on.exit(options(np.tree = old_opt), add = TRUE)
  options(np.tree = TRUE)

  bw.tree <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwmethod = "cv.ls",
    nmulti = 1
  )

  options(np.tree = FALSE)
  serial.at.serial <- npRmpi:::.npcdensbw_eval_only(x, y$y1, bw.serial)$objective
  serial.at.tree <- npRmpi:::.npcdensbw_eval_only(x, y$y1, bw.tree)$objective

  options(np.tree = TRUE)
  tree.at.serial <- npRmpi:::.npcdensbw_eval_only(x, y$y1, bw.serial)$objective
  tree.at.tree <- npRmpi:::.npcdensbw_eval_only(x, y$y1, bw.tree)$objective

  expect_equal(tree.at.serial, serial.at.serial, tolerance = 2e-2)
  expect_equal(tree.at.tree, serial.at.tree, tolerance = 2e-2)
})

test_that("public npcdensbw cv.ls generalized-nn LP route activates with ll == lp parity", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(143)
  n <- 24L
  x <- data.frame(x1 = runif(n))
  y <- data.frame(y1 = x$x1 + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwtype = "generalized_nn",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwtype = "generalized_nn",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1
  )

  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})

test_that("public npcdensbw cv.ls adaptive-nn LP route activates with ll == lp parity", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old_opts <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(144)
  n <- 24L
  x <- data.frame(x1 = runif(n))
  y <- data.frame(y1 = x$x1 + rnorm(n, sd = 0.1))
  degree <- rep.int(1L, ncol(x))

  bw.ll <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "ll",
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1
  )
  bw.lp <- npcdensbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = degree,
    bwtype = "adaptive_nn",
    bwmethod = "cv.ls",
    nmulti = 1,
    itmax = 1
  )

  expect_true(is.finite(bw.ll$fval))
  expect_true(is.finite(bw.lp$fval))
  expect_equal(bw.ll$fval, bw.lp$fval, tolerance = 1e-8)
})
