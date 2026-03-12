library(npRmpi)

test_that("fixed-bwtype plot bootstrap covers supervised wild and unsupervised inid fixed geom families", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(603106)
  n <- 48
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(z = runif(n))
  xidat <- data.frame(x1 = xdat$x, x2 = runif(n))
  yreg <- sin(xdat$x) + rnorm(n, sd = 0.15)
  ydat <- data.frame(y = rnorm(n))
  yscoef <- (0.4 + xdat$x) * cos(2 * pi * zdat$z) + rnorm(n, sd = 0.05)
  yplreg <- sin(2 * pi * zdat$z) + 1.2 * xdat$x + rnorm(n, sd = 0.05)
  yindex <- sin(xidat$x1 + 0.5 * xidat$x2) + rnorm(n, sd = 0.05)

  reg.bws <- list(
    lc = npregbw(
      xdat = xdat,
      ydat = yreg,
      bws = 0.30,
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      regtype = "lc"
    ),
    ll = npregbw(
      xdat = xdat,
      ydat = yreg,
      bws = 0.30,
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      regtype = "ll"
    ),
    lp = npregbw(
      xdat = xdat,
      ydat = yreg,
      bws = 0.30,
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      regtype = "lp",
      degree = 2L,
      basis = "glp",
      bernstein.basis = FALSE
    )
  )
  u.dens.bw <- npudensbw(dat = xdat, bws = 0.30, bandwidth.compute = FALSE, bwtype = "fixed")
  u.dist.bw <- npudistbw(dat = xdat, bws = 0.30, bandwidth.compute = FALSE, bwtype = "fixed")
  c.dens.bw <- npcdensbw(xdat = xdat, ydat = ydat, bws = c(0.35, 0.35), bandwidth.compute = FALSE, bwtype = "fixed")
  c.dist.bw <- npcdistbw(xdat = xdat, ydat = ydat, bws = c(0.35, 0.35), bandwidth.compute = FALSE, bwtype = "fixed")
  sc.bws <- list(
    lc = npscoefbw(
      xdat = xdat,
      zdat = zdat,
      ydat = yscoef,
      bws = 0.30,
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      regtype = "lc"
    ),
    ll = npscoefbw(
      xdat = xdat,
      zdat = zdat,
      ydat = yscoef,
      bws = 0.30,
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      regtype = "ll"
    ),
    lp = npscoefbw(
      xdat = xdat,
      zdat = zdat,
      ydat = yscoef,
      bws = 0.30,
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      regtype = "lp",
      degree = 2L,
      basis = "glp",
      bernstein.basis = FALSE
    )
  )
  pl.bws <- list(
    lc = npplregbw(
      xdat = xdat,
      zdat = zdat,
      ydat = yplreg,
      bws = matrix(c(0.30, 0.30), nrow = 2L, ncol = 1L),
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      regtype = "lc"
    ),
    ll = npplregbw(
      xdat = xdat,
      zdat = zdat,
      ydat = yplreg,
      bws = matrix(c(0.30, 0.30), nrow = 2L, ncol = 1L),
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      regtype = "ll"
    ),
    lp = npplregbw(
      xdat = xdat,
      zdat = zdat,
      ydat = yplreg,
      bws = matrix(c(0.30, 0.30), nrow = 2L, ncol = 1L),
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      regtype = "lp",
      degree = 2L,
      basis = "glp",
      bernstein.basis = FALSE
    )
  )
  si.bws <- list(
    lc = npindexbw(
      xdat = xidat,
      ydat = yindex,
      bws = c(1, 1, 0.30),
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      regtype = "lc"
    ),
    ll = npindexbw(
      xdat = xidat,
      ydat = yindex,
      bws = c(1, 1, 0.30),
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      regtype = "ll"
    ),
    lp = npindexbw(
      xdat = xidat,
      ydat = yindex,
      bws = c(1, 1, 0.30),
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      regtype = "lp",
      degree = 2L,
      basis = "glp",
      bernstein.basis = FALSE
    )
  )

  run_plot <- function(bw, ..., boot.method) {
    suppressWarnings(plot(
      bw,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = boot.method,
      plot.errors.boot.blocklen = 3L,
      plot.errors.boot.num = 5L,
      plot.errors.type = "pointwise",
      neval = 11L,
      ...
    ))
  }

  for (boot.method in c("inid", "fixed", "geom")) {
    for (reg.bw in reg.bws) {
      expect_type(run_plot(reg.bw, xdat = xdat, ydat = yreg, boot.method = boot.method), "list")
    }
    expect_type(run_plot(u.dens.bw, xdat = xdat, boot.method = boot.method), "list")
    expect_type(run_plot(u.dist.bw, xdat = xdat, boot.method = boot.method), "list")
    expect_type(run_plot(c.dens.bw, xdat = xdat, ydat = ydat, view = "fixed", boot.method = boot.method), "list")
    expect_type(run_plot(c.dist.bw, xdat = xdat, ydat = ydat, view = "fixed", boot.method = boot.method), "list")
    for (sc.bw in sc.bws) {
      expect_type(run_plot(sc.bw, xdat = xdat, ydat = yscoef, zdat = zdat, boot.method = boot.method), "list")
    }
    for (pl.bw in pl.bws) {
      expect_type(run_plot(pl.bw, xdat = xdat, ydat = yplreg, zdat = zdat, boot.method = boot.method), "list")
    }
    for (si.bw in si.bws) {
      expect_type(run_plot(si.bw, xdat = xidat, ydat = yindex, boot.method = boot.method), "list")
    }
  }

  for (reg.bw in reg.bws) {
    expect_type(run_plot(reg.bw, xdat = xdat, ydat = yreg, boot.method = "wild"), "list")
  }
  for (sc.bw in sc.bws) {
    expect_type(run_plot(sc.bw, xdat = xdat, ydat = yscoef, zdat = zdat, boot.method = "wild"), "list")
  }
  for (pl.bw in pl.bws) {
    expect_type(run_plot(pl.bw, xdat = xdat, ydat = yplreg, zdat = zdat, boot.method = "wild"), "list")
  }
  for (si.bw in si.bws) {
    expect_type(run_plot(si.bw, xdat = xidat, ydat = yindex, boot.method = "wild"), "list")
  }

  expect_error(
    run_plot(u.dens.bw, xdat = xdat, boot.method = "wild"),
    "not supported"
  )
  expect_error(
    run_plot(u.dist.bw, xdat = xdat, boot.method = "wild"),
    "not supported"
  )
  expect_error(
    run_plot(c.dens.bw, xdat = xdat, ydat = ydat, view = "fixed", boot.method = "wild"),
    "not supported"
  )
  expect_error(
    run_plot(c.dist.bw, xdat = xdat, ydat = ydat, view = "fixed", boot.method = "wild"),
    "not supported"
  )
})
