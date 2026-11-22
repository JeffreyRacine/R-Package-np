## Profile/manual-broadcast plot-route smoke script.
## Launch with mpiexec and npRmpi profile loaded via R_PROFILE_USER.

suppressPackageStartupMessages(library(npRmpi))

if (mpi.comm.rank(0L) == 0L) {
  mpi.bcast.cmd(np.mpi.initialize(), caller.execute = TRUE)
  options(np.messages = FALSE, npRmpi.autodispatch = TRUE)

  stopifnot(isTRUE(npRmpi.init(mode = "spawn", nslaves = mpi.comm.size(1L) - 1L, quiet = TRUE)))

  set.seed(4242)
  n <- as.integer(Sys.getenv("NP_RMPI_PROFILE_PLOT_N", "40"))
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.08)
  dat <- data.frame(y = y, x = x)

  bw <- npregbw(
    y ~ x,
    data = dat,
    bws = 0.2,
    bandwidth.compute = FALSE
  )

  out <- suppressWarnings(
    plot(
      bw,
      output = "data",
      perspective = FALSE,
      errors = "bootstrap",
      bootstrap = "wild",
      B = 5,
      band = "pointwise"
    )
  )
  stopifnot(inherits(out, "np_plot_data"), nrow(out$panels[[1L]]$values) > 0L)

  err <- try(
    mpi.bcast.cmd(
      plot(
        bw,
        output = "data",
        perspective = FALSE,
        errors = "none"
      ),
      caller.execute = TRUE
    ),
    silent = TRUE
  )
  stopifnot(inherits(err, "try-error"))
  msg <- as.character(err)
  stopifnot(any(grepl(
    "plot(...) inside mpi.bcast.cmd(...) is unsupported in canonical SPMD mode",
    msg,
    fixed = TRUE
  )))

  cat("PROFILE_PLOT_ROUTE_OK\n")
  mpi.bcast.cmd(mpi.quit(), caller.execute = TRUE)
}
