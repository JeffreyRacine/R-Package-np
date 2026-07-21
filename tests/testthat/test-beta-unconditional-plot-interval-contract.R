test_that("beta unconditional plot evaluator matches public PDF and CDF fits", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  xdat <- data.frame(x = seq(0.03, 0.97, length.out = 31L))
  exdat <- data.frame(x = seq(0, 1, length.out = 17L))

  for (order in c(2L, 4L, 6L, 8L)) {
    dens.bw <- npudensbw(
      dat = xdat,
      bws = 0.11,
      bandwidth.compute = FALSE,
      ckertype = "beta",
      ckerorder = order,
      ckerbound = "fixed",
      ckerlb = 0,
      ckerub = 1
    )
    dist.bw <- npudistbw(
      dat = xdat,
      bws = 0.11,
      bandwidth.compute = FALSE,
      ckertype = "beta",
      ckerorder = order,
      ckerbound = "fixed",
      ckerlb = 0,
      ckerub = 1
    )

    expect_equal(
      npRmpi:::.np_ksum_unconditional_eval_exact(
        xdat = xdat,
        exdat = exdat,
        bws = dens.bw,
        operator = "normal"
      ),
      as.numeric(npudens(tdat = xdat, edat = exdat, bws = dens.bw)$dens),
      tolerance = 0,
      info = paste("PDF order", order)
    )
    expect_equal(
      npRmpi:::.np_ksum_unconditional_eval_exact(
        xdat = xdat,
        exdat = exdat,
        bws = dist.bw,
        operator = "integral"
      ),
      as.numeric(npudist(tdat = xdat, edat = exdat, bws = dist.bw)$dist),
      tolerance = 0,
      info = paste("CDF order", order)
    )
  }
})

test_that("beta unconditional plot evaluator preserves all bandwidth modes", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  xdat <- data.frame(x = seq(0.03, 0.97, length.out = 31L))

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bws <- if (identical(bwtype, "fixed")) 0.11 else 5
    dens.bw <- npudensbw(
      dat = xdat,
      bws = bws,
      bandwidth.compute = FALSE,
      bwtype = bwtype,
      ckertype = "beta",
      ckerbound = "range"
    )
    dist.bw <- npudistbw(
      dat = xdat,
      bws = bws,
      bandwidth.compute = FALSE,
      bwtype = bwtype,
      ckertype = "beta",
      ckerbound = "range"
    )
    exdat <- data.frame(
      x = seq(dens.bw$ckerlb, dens.bw$ckerub, length.out = 17L)
    )

    expect_equal(
      npRmpi:::.np_ksum_unconditional_eval_exact(
        xdat = xdat,
        exdat = exdat,
        bws = dens.bw,
        operator = "normal"
      ),
      as.numeric(npudens(tdat = xdat, edat = exdat, bws = dens.bw)$dens),
      tolerance = 0,
      info = paste("PDF", bwtype)
    )
    expect_equal(
      npRmpi:::.np_ksum_unconditional_eval_exact(
        xdat = xdat,
        exdat = exdat,
        bws = dist.bw,
        operator = "integral"
      ),
      as.numeric(npudist(tdat = xdat, edat = exdat, bws = dist.bw)$dist),
      tolerance = 0,
      info = paste("CDF", bwtype)
    )
  }
})

test_that("explicit beta density and distribution interval plot routes run", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(1702)
  xdat <- data.frame(x = runif(36L))
  dens.bw <- npudensbw(
    dat = xdat,
    bws = 0.12,
    bandwidth.compute = FALSE,
    ckertype = "beta",
    ckerbound = "range"
  )
  dist.bw <- npudistbw(
    dat = xdat,
    bws = 0.12,
    bandwidth.compute = FALSE,
    ckertype = "beta",
    ckerbound = "range"
  )
  dens.fit <- npudens(tdat = xdat, bws = dens.bw)
  dist.fit <- npudist(tdat = xdat, bws = dist.bw)

  for (errors in c("none", "bootstrap")) {
    extra <- if (identical(errors, "bootstrap")) list(B = 3L) else list()
    dens.plot <- do.call(
      plot,
      c(list(
        x = dens.fit,
        output = "data",
        neval = 13L,
        errors = errors,
        random.seed = 1703L
      ), extra)
    )
    dist.plot <- do.call(
      plot,
      c(list(
        x = dist.fit,
        output = "data",
        neval = 13L,
        errors = errors,
        random.seed = 1703L
      ), extra)
    )

    expect_type(dens.plot, "list")
    expect_type(dist.plot, "list")
    expect_true(all(is.finite(dens.plot$d1$dens)), info = paste("density", errors))
    expect_true(all(is.finite(dist.plot$d1$dist)), info = paste("distribution", errors))
  }
})

test_that("beta NN frozen bootstrap plots and bivariate plot evaluation run", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(1704)
  xdat <- data.frame(x = runif(32L))
  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    dens.bw <- npudensbw(
      dat = xdat,
      bws = 5,
      bandwidth.compute = FALSE,
      bwtype = bwtype,
      ckertype = "beta",
      ckerbound = "range"
    )
    dist.bw <- npudistbw(
      dat = xdat,
      bws = 5,
      bandwidth.compute = FALSE,
      bwtype = bwtype,
      ckertype = "beta",
      ckerbound = "range"
    )

    expect_type(
      plot(
        npudens(tdat = xdat, bws = dens.bw),
        output = "data",
        neval = 11L,
        errors = "bootstrap",
        boot_control = np_boot_control(nonfixed = "frozen"),
        B = 2L,
        random.seed = 1705L
      ),
      "list"
    )
    expect_type(
      plot(
        npudist(tdat = xdat, bws = dist.bw),
        output = "data",
        neval = 11L,
        errors = "bootstrap",
        boot_control = np_boot_control(nonfixed = "frozen"),
        B = 2L,
        random.seed = 1705L
      ),
      "list"
    )
  }

  xy <- data.frame(x = runif(30L), y = runif(30L))
  xy.bw <- npudensbw(
    dat = xy,
    bws = c(0.13, 0.15),
    bandwidth.compute = FALSE,
    ckertype = "beta",
    ckerbound = "range"
  )
  xy.plot <- plot(
    npudens(tdat = xy, bws = xy.bw),
    output = "data",
    neval = 5L,
    errors = "bootstrap",
    B = 2L,
    random.seed = 1706L
  )

  expect_type(xy.plot, "list")
  expect_true(all(is.finite(xy.plot$d1$dens)))
})

test_that("ordinary-kernel plot evaluator normalization is unchanged", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  xdat <- data.frame(x = seq(-1, 1, length.out = 29L))
  exdat <- data.frame(x = seq(-0.9, 0.9, length.out = 13L))
  bws <- npudensbw(
    dat = xdat,
    bws = 0.19,
    bandwidth.compute = FALSE,
    ckertype = "gaussian"
  )

  expect_equal(
    npRmpi:::.np_ksum_unconditional_eval_exact(
      xdat = xdat,
      exdat = exdat,
      bws = bws,
      operator = "normal"
    ),
    as.numeric(npudens(tdat = xdat, edat = exdat, bws = bws)$dens),
    tolerance = 0
  )
})
