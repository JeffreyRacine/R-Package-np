test_that("nearest-neighbor frozen bootstrap plots run across regression and unsupervised families", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(32908)
  n <- 24
  x <- data.frame(x = runif(n))
  y <- sin(2 * pi * x$x) + rnorm(n, sd = 0.08)
  yframe <- data.frame(y = y)

  run_plot <- function(bw, ...) {
    suppressWarnings(plot(
      bw,
      output = "data",
      perspective = FALSE,
      errors = "bootstrap",
      boot_control = np_boot_control(nonfixed = "frozen"),
      B = 3,
      ...
    ))
  }

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw.val <- if (identical(bt, "adaptive_nn")) 5 else 2
    rbw <- do.call(npregbw, list(xdat = x, ydat = y, regtype = "ll",
                                 bws = bw.val, bwtype = bt,
                                 bandwidth.compute = FALSE))

    for (boot.method in c("inid", "fixed", "geom")) {
      expect_type(
        run_plot(rbw, xdat = x, ydat = y, bootstrap = boot.method),
        "list"
      )

      ubw <- do.call(npudensbw, list(dat = x, bws = bw.val, bwtype = bt,
                                     bandwidth.compute = FALSE))
      expect_type(run_plot(ubw, bootstrap = boot.method), "list")

      dbw <- do.call(npudistbw, list(dat = x, bws = bw.val, bwtype = bt,
                                     bandwidth.compute = FALSE))
      expect_type(run_plot(dbw, bootstrap = boot.method), "list")

      cbw <- do.call(npcdensbw, list(xdat = x, ydat = yframe,
                                     bws = rep.int(bw.val, 2L), bwtype = bt,
                                     bandwidth.compute = FALSE))
      expect_type(
        run_plot(cbw, xdat = x, ydat = yframe, view = "fixed", bootstrap = boot.method),
        "list"
      )

      cdbw <- do.call(npcdistbw, list(xdat = x, ydat = yframe,
                                      bws = rep.int(bw.val, 2L), bwtype = bt,
                                      bandwidth.compute = FALSE))
      expect_type(
        run_plot(cdbw, xdat = x, ydat = yframe, view = "fixed", bootstrap = boot.method),
        "list"
      )
    }
  }
})
