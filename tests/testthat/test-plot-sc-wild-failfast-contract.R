quiet_capture <- function(expr) {
  out <- NULL
  sinkfile <- tempfile()
  on.exit(unlink(sinkfile), add = TRUE)
  invisible(capture.output(out <- eval.parent(substitute(expr)), file = sinkfile))
  out
}

test_that("smooth-coefficient wild bootstrap returns finite plot-data errors", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(9711)
  n <- 40L
  xdat <- data.frame(x = runif(n))
  zdat <- data.frame(z = runif(n))
  ydat <- sin(2 * pi * xdat$x) * zdat$z + rnorm(n, sd = 0.1)

  bws <- npscoefbw(
    xdat = xdat,
    zdat = zdat,
    ydat = ydat,
    bws = 0.2,
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = 2L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    ckertype = "gaussian"
  )

  out <- expect_no_error(
    quiet_capture(
      suppressWarnings(
        plot(
          bws,
          xdat = xdat,
          ydat = ydat,
          zdat = zdat,
          output = "data",
          perspective = FALSE,
          errors = "bootstrap",
          bootstrap = "wild",
          B = 9L
        )
      )
    )
  )
  expect_type(out, "list")
  expect_true(all(is.finite(out$mean)))
  expect_true(all(is.finite(out$merr[, 1:2, drop = FALSE])))
})
