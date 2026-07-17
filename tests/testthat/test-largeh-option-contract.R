test_that("np.largeh toggles continuous large-bandwidth shortcuts", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  } else {
    skip_if_not(isTRUE(getOption("npRmpi.mpi.initialized", FALSE)),
                "MPI pool unavailable")
  }
  on.exit(close_mpi_slaves(force = FALSE), add = TRUE)

  old.options <- options(np.tree = FALSE, np.largeh = TRUE, np.largelambda = TRUE)
  on.exit(options(old.options), add = TRUE)

  npregbw <- getFromNamespace("npregbw", "npRmpi")
  eval_only <- getFromNamespace(".npregbw_eval_only", "npRmpi")

  set.seed(42)
  n <- 120L
  xdat <- data.frame(x1 = runif(n), x2 = runif(n))
  ydat <- xdat$x1 + xdat$x2 + rnorm(n)

  bws <- npregbw(
    xdat = xdat,
    ydat = ydat,
    regtype = "ll",
    bwmethod = "cv.ls",
    ckertype = "epanechnikov",
    bandwidth.compute = FALSE,
    bws = c(1e8, 1e8)
  )

  options(np.largeh = TRUE)
  enabled <- eval_only(xdat = xdat, ydat = ydat, bws = bws)
  options(np.largeh = FALSE)
  disabled <- eval_only(xdat = xdat, ydat = ydat, bws = bws)
  options(np.largeh = TRUE)
  reenabled <- eval_only(xdat = xdat, ydat = ydat, bws = bws)

  expect_equal(enabled$objective, disabled$objective, tolerance = 1e-10)
  expect_equal(enabled$objective, reenabled$objective, tolerance = 1e-10)
  expect_equal(enabled$num.feval.fast, 1)
  expect_equal(disabled$num.feval.fast, 0)
  expect_equal(reenabled$num.feval.fast, 1)
})

test_that("np.largeh continuous shortcut eligibility rejects unsafe kernels", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  } else {
    skip_if_not(isTRUE(getOption("npRmpi.mpi.initialized", FALSE)),
                "MPI pool unavailable")
  }
  on.exit(close_mpi_slaves(force = FALSE), add = TRUE)

  old.options <- options(np.tree = FALSE, np.largeh = TRUE, np.largelambda = TRUE)
  on.exit(options(old.options), add = TRUE)

  npregbw <- getFromNamespace("npregbw", "npRmpi")
  eval_only <- getFromNamespace(".npregbw_eval_only", "npRmpi")

  set.seed(4401)
  n <- 100L
  xdat <- data.frame(x = runif(n))
  ydat <- sin(xdat$x) + rnorm(n, sd = 0.05)
  cases <- data.frame(
    ckertype = c("gaussian", "gaussian", "epanechnikov", "epanechnikov",
                 "uniform"),
    ckerorder = c(2L, 4L, 2L, 4L, 2L),
    expected.fast = c(1L, 0L, 1L, 0L, 1L),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(cases))) {
    row <- cases[i, ]
    bws <- suppressWarnings(npregbw(
      xdat = xdat,
      ydat = ydat,
      regtype = "ll",
      bwmethod = "cv.ls",
      ckertype = row$ckertype,
      ckerorder = row$ckerorder,
      bandwidth.compute = FALSE,
      bws = 1e8
    ))

    options(np.largeh = TRUE)
    enabled <- eval_only(xdat = xdat, ydat = ydat, bws = bws)
    options(np.largeh = FALSE)
    disabled <- eval_only(xdat = xdat, ydat = ydat, bws = bws)
    options(np.largeh = TRUE)

    expect_equal(enabled$objective, disabled$objective, tolerance = 1e-10)
    expect_equal(as.integer(enabled$num.feval.fast), row$expected.fast)
    expect_equal(as.integer(disabled$num.feval.fast), 0L)
  }
})

test_that("np.largelambda toggles discrete upper-lambda shortcuts", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  } else {
    skip_if_not(isTRUE(getOption("npRmpi.mpi.initialized", FALSE)),
                "MPI pool unavailable")
  }
  on.exit(close_mpi_slaves(force = FALSE), add = TRUE)

  old.options <- options(np.tree = FALSE, np.largeh = TRUE, np.largelambda = TRUE)
  on.exit(options(old.options), add = TRUE)

  npregbw <- getFromNamespace("npregbw", "npRmpi")
  eval_only <- getFromNamespace(".npregbw_eval_only", "npRmpi")

  set.seed(43)
  n <- 120L
  xdat <- data.frame(z = factor(sample(letters[1:3], n, replace = TRUE)))
  ydat <- rnorm(n)

  bws <- npregbw(
    xdat = xdat,
    ydat = ydat,
    regtype = "lc",
    bwmethod = "cv.ls",
    ukertype = "aitchisonaitken",
    bandwidth.compute = FALSE,
    bws = 2/3
  )

  options(np.largelambda = TRUE)
  enabled <- eval_only(xdat = xdat, ydat = ydat, bws = bws)
  options(np.largelambda = FALSE)
  disabled <- eval_only(xdat = xdat, ydat = ydat, bws = bws)
  options(np.largelambda = TRUE)
  reenabled <- eval_only(xdat = xdat, ydat = ydat, bws = bws)

  expect_equal(enabled$objective, disabled$objective, tolerance = 1e-10)
  expect_equal(enabled$objective, reenabled$objective, tolerance = 1e-10)
  expect_equal(enabled$num.feval.fast, 1)
  expect_equal(disabled$num.feval.fast, 0)
  expect_equal(reenabled$num.feval.fast, 1)
})

test_that("np.largeh toggles one-step continuous bandwidth-search fast counts", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  } else {
    skip_if_not(isTRUE(getOption("npRmpi.mpi.initialized", FALSE)),
                "MPI pool unavailable")
  }
  on.exit(close_mpi_slaves(force = FALSE), add = TRUE)

  old.options <- options(np.messages = FALSE, np.tree = FALSE,
                         np.largeh = TRUE, np.largelambda = TRUE)
  on.exit(options(old.options), add = TRUE)

  npregbw <- getFromNamespace("npregbw", "npRmpi")

  set.seed(101)
  n <- 80L
  xdat <- data.frame(x1 = runif(n), x2 = runif(n))
  ydat <- xdat$x1 + xdat$x2 + rnorm(n)

  bws <- npregbw(
    xdat = xdat,
    ydat = ydat,
    regtype = "ll",
    bwmethod = "cv.ls",
    ckertype = "epanechnikov",
    bandwidth.compute = FALSE,
    bws = c(1e8, 1e8)
  )

  options(np.largeh = TRUE)
  enabled <- npregbw(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    regtype = "ll",
    bwmethod = "cv.ls",
    ckertype = "epanechnikov",
    bandwidth.compute = TRUE,
    nmulti = 1L,
    itmax = 1L
  )

  options(np.largeh = FALSE)
  disabled <- npregbw(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    regtype = "ll",
    bwmethod = "cv.ls",
    ckertype = "epanechnikov",
    bandwidth.compute = TRUE,
    nmulti = 1L,
    itmax = 1L
  )

  expect_equal(enabled$fval, disabled$fval, tolerance = 1e-10)
  expect_gt(as.numeric(enabled$num.feval.fast[1L]), 0)
  expect_gt(as.numeric(enabled$num.feval.fast[1L]),
            as.numeric(disabled$num.feval.fast[1L]))
})

test_that("np.largeh and np.largelambda both gate mixed fast objective rows", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  } else {
    skip_if_not(isTRUE(getOption("npRmpi.mpi.initialized", FALSE)),
                "MPI pool unavailable")
  }
  on.exit(close_mpi_slaves(force = FALSE), add = TRUE)

  old.options <- options(np.messages = FALSE, np.tree = FALSE,
                         np.largeh = TRUE, np.largelambda = TRUE)
  on.exit(options(old.options), add = TRUE)

  npregbw <- getFromNamespace("npregbw", "npRmpi")
  npreg_eval_only <- getFromNamespace(".npregbw_eval_only", "npRmpi")
  npscoefbw <- getFromNamespace("npscoefbw", "npRmpi")
  npscoef_eval_only <- getFromNamespace(".npscoefbw_eval_only", "npRmpi")

  set.seed(107)
  n <- 80L
  xdat <- data.frame(
    x = runif(n),
    z = factor(sample(letters[1:3], n, replace = TRUE))
  )
  ydat <- rnorm(n)

  rbw <- npregbw(
    xdat = xdat,
    ydat = ydat,
    regtype = "lc",
    bwmethod = "cv.ls",
    ckertype = "epanechnikov",
    ukertype = "aitchisonaitken",
    bandwidth.compute = FALSE,
    bws = c(1e8, 2 / 3)
  )

  options(np.largeh = TRUE, np.largelambda = TRUE)
  reg.both <- npreg_eval_only(xdat = xdat, ydat = ydat, bws = rbw)
  options(np.largeh = FALSE, np.largelambda = TRUE)
  reg.noh <- npreg_eval_only(xdat = xdat, ydat = ydat, bws = rbw)
  options(np.largeh = TRUE, np.largelambda = FALSE)
  reg.nolam <- npreg_eval_only(xdat = xdat, ydat = ydat, bws = rbw)
  options(np.largeh = TRUE, np.largelambda = TRUE)
  reg.both2 <- npreg_eval_only(xdat = xdat, ydat = ydat, bws = rbw)

  expect_equal(reg.both$objective, reg.noh$objective, tolerance = 1e-10)
  expect_equal(reg.both$objective, reg.nolam$objective, tolerance = 1e-10)
  expect_equal(reg.both$objective, reg.both2$objective, tolerance = 1e-10)
  expect_equal(as.numeric(reg.both$num.feval.fast[1L]), 1)
  expect_equal(as.numeric(reg.noh$num.feval.fast[1L]), 0)
  expect_equal(as.numeric(reg.nolam$num.feval.fast[1L]), 0)
  expect_equal(as.numeric(reg.both2$num.feval.fast[1L]), 1)

  set.seed(108)
  sxdat <- data.frame(x = runif(n))
  szdat <- data.frame(
    z = runif(n),
    g = factor(sample(letters[1:3], n, replace = TRUE))
  )
  sydat <- rnorm(n)

  sbw <- npscoefbw(
    xdat = sxdat,
    zdat = szdat,
    ydat = sydat,
    bwmethod = "cv.ls",
    regtype = "lc",
    bandwidth.compute = FALSE,
    bws = c(1e8, 2 / 3)
  )

  options(np.largeh = TRUE, np.largelambda = TRUE)
  sc.both <- npscoef_eval_only(xdat = sxdat, zdat = szdat, ydat = sydat, bws = sbw)
  options(np.largeh = FALSE, np.largelambda = TRUE)
  sc.noh <- npscoef_eval_only(xdat = sxdat, zdat = szdat, ydat = sydat, bws = sbw)
  options(np.largeh = TRUE, np.largelambda = FALSE)
  sc.nolam <- npscoef_eval_only(xdat = sxdat, zdat = szdat, ydat = sydat, bws = sbw)
  options(np.largeh = TRUE, np.largelambda = TRUE)
  sc.both2 <- npscoef_eval_only(xdat = sxdat, zdat = szdat, ydat = sydat, bws = sbw)

  expect_equal(sc.both$objective, sc.noh$objective, tolerance = 1e-10)
  expect_equal(sc.both$objective, sc.nolam$objective, tolerance = 1e-10)
  expect_equal(sc.both$objective, sc.both2$objective, tolerance = 1e-10)
  expect_equal(as.numeric(sc.both$num.feval.fast[1L]), 1)
  expect_equal(as.numeric(sc.noh$num.feval.fast[1L]), 0)
  expect_equal(as.numeric(sc.nolam$num.feval.fast[1L]), 0)
  expect_equal(as.numeric(sc.both2$num.feval.fast[1L]), 1)
})

test_that("np.largeh toggles conditional density eval-only fast counts", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "old_opts <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = TRUE, np.largelambda = TRUE)",
      "on.exit(options(old_opts), add = TRUE)",
      "npRmpi.init(nslaves = 1, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "set.seed(105)",
      "n <- 80L",
      "xdat <- data.frame(x = runif(n))",
      "ydat <- data.frame(y = runif(n))",
      "bws <- npcdensbw(xdat = xdat, ydat = ydat, regtype = 'lc', bwmethod = 'cv.ml',",
      "                  cxkertype = 'epanechnikov', cykertype = 'epanechnikov',",
      "                  bandwidth.compute = FALSE, bws = c(1e8, 1e8))",
      "options(np.largeh = TRUE)",
      "enabled <- npRmpi:::.npcdensbw_eval_only(xdat = xdat, ydat = ydat, bws = bws)",
      "options(np.largeh = FALSE)",
      "disabled <- npRmpi:::.npcdensbw_eval_only(xdat = xdat, ydat = ydat, bws = bws)",
      "options(np.largeh = TRUE)",
      "reenabled <- npRmpi:::.npcdensbw_eval_only(xdat = xdat, ydat = ydat, bws = bws)",
      "stopifnot(isTRUE(all.equal(enabled$objective, disabled$objective, tolerance = 1e-10)))",
      "stopifnot(isTRUE(all.equal(enabled$objective, reenabled$objective, tolerance = 1e-10)))",
      "stopifnot(identical(as.numeric(enabled$num.feval.fast[1L]), 1))",
      "stopifnot(identical(as.numeric(disabled$num.feval.fast[1L]), 0))",
      "stopifnot(identical(as.numeric(reenabled$num.feval.fast[1L]), 1))",
      "cat('NPCDENS_LARGEH_OPTION_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NPCDENS_LARGEH_OPTION_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})

test_that("np.largeh toggles unconditional distribution bandwidth fast counts", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "old_opts <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = TRUE, np.largelambda = TRUE)",
      "on.exit(options(old_opts), add = TRUE)",
      "npRmpi.init(nslaves = 1, quiet = TRUE)",
      "on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)",
      "set.seed(104)",
      "n <- 80L",
      "dat <- data.frame(y = runif(n))",
      "bws <- npudistbw(dat = dat, bwmethod = 'cv.cdf', ckertype = 'epanechnikov',",
      "                  bandwidth.compute = FALSE, ngrid = 20L)",
      "bws$bw[] <- 1e8",
      "bws$bandwidth[] <- 1e8",
      "bws$sfactor[] <- 1e8",
      "options(np.largeh = TRUE)",
      "enabled <- npudistbw(dat = dat, bws = bws, bandwidth.compute = TRUE,",
      "                     nmulti = 1L, itmax = 1L, ngrid = 20L)",
      "options(np.largeh = FALSE)",
      "disabled <- npudistbw(dat = dat, bws = bws, bandwidth.compute = TRUE,",
      "                      nmulti = 1L, itmax = 1L, ngrid = 20L)",
      "stopifnot(isTRUE(all.equal(enabled$fval, disabled$fval, tolerance = 1e-10)))",
      "stopifnot(as.numeric(enabled$num.feval.fast[1L]) > 0)",
      "stopifnot(as.numeric(enabled$num.feval.fast[1L]) > as.numeric(disabled$num.feval.fast[1L]))",
      "cat('NPUDIST_LARGEH_OPTION_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl("NPUDIST_LARGEH_OPTION_OK", res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
