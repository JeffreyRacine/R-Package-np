quiet_capture <- function(expr) {
  out <- NULL
  sinkfile <- tempfile()
  on.exit(unlink(sinkfile), add = TRUE)
  invisible(capture.output(out <- eval.parent(substitute(expr)), file = sinkfile))
  out
}

test_that("smooth-coefficient wild bootstrap fails fast with explicit diagnostics", {
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

  expect_error(
    quiet_capture(
      suppressWarnings(
        plot(
          bws,
          xdat = xdat,
          ydat = ydat,
          zdat = zdat,
          plot.behavior = "data",
          perspective = FALSE,
          plot.errors.method = "bootstrap",
          plot.errors.boot.method = "wild",
          plot.errors.boot.num = 9L
        )
      )
    ),
    "unsupported for smooth coefficient bootstrap in npRmpi canonical SPMD mode",
    fixed = TRUE
  )
})
