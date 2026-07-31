test_that("regression and index consumers require coherent canonical state", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  validate <- getFromNamespace(
    "npValidatedConditionalRegSpec",
    "npRmpi"
  )

  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(20260731)
  n <- 32L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- x$x1 - 0.4 * x$x2 + rnorm(n, sd = 0.1)

  for (regtype in c("lc", "ll", "lp")) {
    args <- list(
      xdat = x,
      ydat = y,
      bws = c(0.3, 0.3),
      bandwidth.compute = FALSE,
      regtype = regtype
    )
    if (identical(regtype, "lp"))
      args$degree <- c(2L, 1L)
    bw <- do.call(npregbw, args)
    spec <- validate(bw, where = "test", ncon.field = "ncon")

    expect_identical(spec$regtype, regtype)
    expect_identical(spec$regtype.engine, bw$regtype.engine)
    expect_identical(spec$basis.engine, bw$basis.engine)
    expect_identical(spec$degree.engine, bw$degree.engine)
    expect_identical(
      spec$bernstein.basis.engine,
      bw$bernstein.basis.engine
    )
  }

  bw <- npregbw(
    xdat = x,
    ydat = y,
    bws = c(0.3, 0.3),
    bandwidth.compute = FALSE
  )
  missing.state <- bw
  names(missing.state)[names(missing.state) == "regtype.engine"] <-
    "regtype.engine.shadow"

  expect_error(
    npreg(bws = missing.state, txdat = x, tydat = y),
    "npreg is missing required field 'regtype.engine'",
    fixed = TRUE
  )
  expect_error(
    npreghat(
      bws = missing.state,
      txdat = x,
      y = y,
      output = "apply"
    ),
    "npreghat is missing required field 'regtype.engine'",
    fixed = TRUE
  )
  expect_error(
    plot(
      missing.state,
      xdat = x,
      ydat = y,
      output = "data",
      perspective = FALSE,
      neval = 5L
    ),
    "plot.rbandwidth() is missing required field 'regtype.engine'",
    fixed = TRUE
  )

  mismatched <- bw
  mismatched$regtype <- "lp"
  mismatched$degree <- c(0L, 0L)
  expect_error(
    npreg(bws = mismatched, txdat = x, tydat = y),
    "incoherent public and canonical engine metadata for field 'regtype.engine'",
    fixed = TRUE
  )

  index.bw <- npindexbw(
    xdat = x,
    ydat = y,
    bws = c(1, 1, 0.3),
    bandwidth.compute = FALSE,
    method = "ichimura",
    regtype = "ll"
  )
  index.spec <- validate(index.bw, where = "test", ncon = 1L)
  expect_identical(index.spec$regtype, "ll")
  expect_identical(index.spec$regtype.engine, "lp")
  expect_identical(index.spec$degree.engine, 1L)

  index.missing <- index.bw
  names(index.missing)[names(index.missing) == "basis.engine"] <-
    "basis.engine.shadow"
  expect_error(
    npindex(bws = index.missing, txdat = x, tydat = y),
    "npindex is missing required field 'basis.engine'",
    fixed = TRUE
  )

  index.incoherent <- index.bw
  index.incoherent$regtype.engine <- "lc"
  expect_error(
    npindex(bws = index.incoherent, txdat = x, tydat = y),
    "local-constant engine state: 'degree.engine' must be all zero",
    fixed = TRUE
  )
})
