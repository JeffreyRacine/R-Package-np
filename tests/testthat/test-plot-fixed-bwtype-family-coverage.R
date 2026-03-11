test_that("fixed-bwtype plot bootstrap covers regression and unsupervised families for inid fixed and geom", {
  skip_if_not_installed("np")

  set.seed(603106)
  n <- 48
  xdat <- data.frame(x = rnorm(n))
  yreg <- sin(xdat$x) + rnorm(n, sd = 0.15)
  ydat <- data.frame(y = rnorm(n))

  reg.bw <- npregbw(
    xdat = xdat,
    ydat = yreg,
    bws = 0.30,
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    regtype = "ll"
  )
  u.dens.bw <- npudensbw(dat = xdat, bws = 0.30, bandwidth.compute = FALSE, bwtype = "fixed")
  u.dist.bw <- npudistbw(dat = xdat, bws = 0.30, bandwidth.compute = FALSE, bwtype = "fixed")
  c.dens.bw <- npcdensbw(xdat = xdat, ydat = ydat, bws = c(0.35, 0.35), bandwidth.compute = FALSE, bwtype = "fixed")
  c.dist.bw <- npcdistbw(xdat = xdat, ydat = ydat, bws = c(0.35, 0.35), bandwidth.compute = FALSE, bwtype = "fixed")

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
    expect_type(run_plot(reg.bw, xdat = xdat, ydat = yreg, boot.method = boot.method), "list")
    expect_type(run_plot(u.dens.bw, xdat = xdat, boot.method = boot.method), "list")
    expect_type(run_plot(u.dist.bw, xdat = xdat, boot.method = boot.method), "list")
    expect_type(run_plot(c.dens.bw, xdat = xdat, ydat = ydat, view = "fixed", boot.method = boot.method), "list")
    expect_type(run_plot(c.dist.bw, xdat = xdat, ydat = ydat, view = "fixed", boot.method = boot.method), "list")
  }
})
