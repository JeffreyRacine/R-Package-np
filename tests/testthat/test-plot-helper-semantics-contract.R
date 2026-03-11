quiet_capture <- function(expr) {
  out <- NULL
  sinkfile <- tempfile()
  on.exit(unlink(sinkfile), add = TRUE)
  invisible(capture.output(out <- eval.parent(substitute(expr)), file = sinkfile))
  out
}

flatten_plot_mean <- function(plot_obj) {
  vals <- unlist(lapply(plot_obj, function(comp) {
    keys <- c("mean", "dens", "dist", "condens", "condist")
    hit <- keys[keys %in% names(comp)]
    if (!length(hit))
      return(numeric(0))
    as.double(comp[[hit[1L]]])
  }), use.names = FALSE)
  vals[is.finite(vals)]
}

run_bootstrap_plot <- function(bw,
                               xdat = NULL,
                               ydat = NULL,
                               zdat = NULL,
                               boot_num = 7L,
                               boot_method = "inid") {
  args <- list(
    bw,
    plot.behavior = "data",
    perspective = FALSE,
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = boot_method,
    plot.errors.boot.num = boot_num
  )
  if (!is.null(xdat))
    args$xdat <- xdat
  if (!is.null(ydat))
    args$ydat <- ydat
  if (!is.null(zdat))
    args$zdat <- zdat

  quiet_capture(suppressWarnings(do.call(plot, args)))
}

with_spawned_plot_context <- function(code) {
  if (!spawn_mpi_slaves())
    skip("Could not spawn MPI slaves")

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  force(code)
}

test_that("unsupervised plot helper families use method-defining kernel options", {
  with_spawned_plot_context({
    set.seed(9401)
    n <- 48L
    xdat <- data.frame(x = runif(n))
    ydat <- data.frame(y = runif(n))

    ud.gaussian <- npudensbw(
      dat = xdat,
      bws = 0.24,
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      ckertype = "gaussian"
    )
    ud.epan <- npudensbw(
      dat = xdat,
      bws = 0.24,
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      ckertype = "epanechnikov",
      ckerorder = 2L
    )

    set.seed(9402)
    ud.out.g <- run_bootstrap_plot(ud.gaussian, xdat = xdat, boot_method = "inid", boot_num = 7L)
    set.seed(9402)
    ud.out.e <- run_bootstrap_plot(ud.epan, xdat = xdat, boot_method = "inid", boot_num = 7L)
    ud.mean.g <- flatten_plot_mean(ud.out.g)
    ud.mean.e <- flatten_plot_mean(ud.out.e)
    expect_gt(length(ud.mean.g), 0L)
    expect_equal(length(ud.mean.g), length(ud.mean.e))
    expect_gt(max(abs(ud.mean.g - ud.mean.e)), 1e-6)

    uf.gaussian <- npudistbw(
      dat = xdat,
      bws = 0.24,
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      ckertype = "gaussian"
    )
    uf.epan <- npudistbw(
      dat = xdat,
      bws = 0.24,
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      ckertype = "epanechnikov",
      ckerorder = 2L
    )

    set.seed(9403)
    uf.out.g <- run_bootstrap_plot(uf.gaussian, xdat = xdat, boot_method = "inid", boot_num = 7L)
    set.seed(9403)
    uf.out.e <- run_bootstrap_plot(uf.epan, xdat = xdat, boot_method = "inid", boot_num = 7L)
    uf.mean.g <- flatten_plot_mean(uf.out.g)
    uf.mean.e <- flatten_plot_mean(uf.out.e)
    expect_gt(length(uf.mean.g), 0L)
    expect_equal(length(uf.mean.g), length(uf.mean.e))
    expect_gt(max(abs(uf.mean.g - uf.mean.e)), 1e-6)

    cd.gaussian <- npcdensbw(
      xdat = xdat,
      ydat = ydat,
      bws = c(0.24, 0.24),
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      cxkertype = "gaussian",
      cykertype = "gaussian"
    )
    cd.epan <- npcdensbw(
      xdat = xdat,
      ydat = ydat,
      bws = c(0.24, 0.24),
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      cxkertype = "epanechnikov",
      cykertype = "epanechnikov",
      cxkerorder = 2L,
      cykerorder = 2L
    )

    set.seed(9404)
    cd.out.g <- run_bootstrap_plot(cd.gaussian, xdat = xdat, ydat = ydat, boot_method = "inid", boot_num = 7L)
    set.seed(9404)
    cd.out.e <- run_bootstrap_plot(cd.epan, xdat = xdat, ydat = ydat, boot_method = "inid", boot_num = 7L)
    cd.mean.g <- flatten_plot_mean(cd.out.g)
    cd.mean.e <- flatten_plot_mean(cd.out.e)
    expect_gt(length(cd.mean.g), 0L)
    expect_equal(length(cd.mean.g), length(cd.mean.e))
    expect_gt(max(abs(cd.mean.g - cd.mean.e)), 1e-6)

    cf.gaussian <- npcdistbw(
      xdat = xdat,
      ydat = ydat,
      bws = c(0.24, 0.24),
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      cxkertype = "gaussian",
      cykertype = "gaussian"
    )
    cf.epan <- npcdistbw(
      xdat = xdat,
      ydat = ydat,
      bws = c(0.24, 0.24),
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      cxkertype = "epanechnikov",
      cykertype = "epanechnikov",
      cxkerorder = 2L,
      cykerorder = 2L
    )

    set.seed(9405)
    cf.out.g <- run_bootstrap_plot(cf.gaussian, xdat = xdat, ydat = ydat, boot_method = "inid", boot_num = 7L)
    set.seed(9405)
    cf.out.e <- run_bootstrap_plot(cf.epan, xdat = xdat, ydat = ydat, boot_method = "inid", boot_num = 7L)
    cf.mean.g <- flatten_plot_mean(cf.out.g)
    cf.mean.e <- flatten_plot_mean(cf.out.e)
    expect_gt(length(cf.mean.g), 0L)
    expect_equal(length(cf.mean.g), length(cf.mean.e))
    expect_gt(max(abs(cf.mean.g - cf.mean.e)), 1e-6)
  })
})
