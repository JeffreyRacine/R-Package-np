test_that("CDF corrected slice drawings use each panel's own coordinates", {
  if (TRUE) {
    skip_on_cran()
    if (!spawn_mpi_slaves(as.integer(Sys.getenv("N5_SLAVES", "1"))))
      stop("required MPI pool unavailable")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  dev.control("enable")
  dotted <- function(rp) {
    out <- list()
    for (op in rp[[1L]]) {
      if (!identical(op[[2L]][[1L]]$name, "C_plotXY")) next
      a <- op[[2L]]
      if (identical(a[[5L]], "dotted") || isTRUE(a[[5L]] == 3))
        out[[length(out) + 1L]] <- a[[2L]]
    }
    out
  }
  set.seed(13)
  x <- data.frame(a = runif(90), b = rnorm(90, 5, 2), c = runif(90, -8, -2))
  full <- identical(Sys.getenv("N5_FULL"), "true")
  snapshots <- list()
  for (q in if (full) 1:2 else 2L) {
    d <- x[seq_len(q)]
    bw <- npudistbw(dat = d, bws = c(.15, .8, .7)[seq_len(q)],
                   bandwidth.compute = FALSE)
    f <- npudist(bws = bw, tdat = d)
    frozen <- serialize(f, NULL)
    for (cs in c(TRUE, FALSE)) for (band in if (full) c("pointwise", "all") else "pointwise") {
      args <- list(f, neval = 6L, perspective = FALSE, errors = "asymptotic",
                   band = band, center = "bias-corrected", common_scale = cs)
      if (band == "all") args <- c(args, list(gaussian_draws = 199L, random_seed = 19L))
      set.seed(73)
      out <- do.call(plot, c(args, list(output = "plot-data")))
      rng <- .Random.seed
      dl <- dotted(recordPlot())
      expect_length(dl, q)
      for (k in seq_len(q)) {
        dat <- as.data.frame(out, panel = k, band = if (band == "all") "simultaneous" else band)
        expect_equal(dl[[k]]$x, dat[[paste0("eval_", k)]], tolerance = 0)
        expect_equal(dl[[k]]$y, dat$bias_corrected, tolerance = 0)
      }
      set.seed(73)
      data <- do.call(plot, c(args, list(output = "data")))
      expect_identical(out, data)
      expect_identical(.Random.seed, rng)
      expect_identical(serialize(f, NULL), frozen)
      snapshots[[paste(q, cs, band)]] <- list(data = data, rng = rng)
      if (full) {
        plot(out)
        replay <- dotted(recordPlot())
        expect_equal(replay, dl, tolerance = 0)
        expect_identical(.Random.seed, rng)
        args$center <- "estimate"
        if (band == "all") expect_warning(do.call(plot, args), "below the recommended")
        else do.call(plot, args)
        expect_length(dotted(recordPlot()), 0L)
      }
    }
  }
  if (full) {
    bw <- npudistbw(dat = x, bws = c(.15, .8, .7), bandwidth.compute = FALSE)
    f <- npudist(bws = bw, tdat = x)
    expect_error(plot(f, neval = 6L, perspective = FALSE, errors = "asymptotic",
                      band = "pointwise", center = "bias-corrected", output = "data"),
                 "at most two numeric")
  }
  if (nzchar(Sys.getenv("N5_SNAPSHOT"))) saveRDS(snapshots, Sys.getenv("N5_SNAPSHOT"))
})
