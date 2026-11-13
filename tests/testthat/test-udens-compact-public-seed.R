test_that("public compact bands use the recorded child seed", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  d <- data.frame(x = qbeta((seq_len(20L) - .35) / 21, 2, 3))
  for (kernel in c("uniform", "epanechnikov")) {
    b <- npudensbw(dat = d, bws = .18, bandwidth.compute = FALSE,
      ckertype = kernel, ckerbound = "fixed", ckerlb = 0, ckerub = 1)
    for (variant in c("asymptotic", "symmetric-sd", "asymmetric-sd")) {
      args <- list(x = b, xdat = d, errors = if (variant == "asymptotic")
        "asymptotic" else "bootstrap", center = "bias-corrected", band = "all",
        neval = 3L, output = "data", random.seed = 19L, alpha = .1)
      if (variant == "asymptotic") args$gaussian.draws <- 9L else {
        args$bootstrap <- "inid"; args$B <- 9L
        args$boot.control <- np_boot_control(calibration = variant)
      }
      actual <- suppressWarnings(do.call(plot, args))
      panel <- actual$panels[[1L]]
      # Direct seed-explicit replay isolates public wiring from stream generation.
      expected <- suppressWarnings(.np_udens_compact_hidden(d$x, .18,
        panel$grid[[1L]], kernel, .1, "all", 9L, 19L, variant = variant))
      expect_identical(actual$request$random_seed, 19L)
      expect_identical(panel$values$bias_corrected, expected$center)
      for (band in names(expected$bands)) {
        ci <- panel$intervals[panel$intervals$band == band, ]
        expect_identical(ci$lower, expected$bands[[band]]$lower)
        expect_identical(ci$upper, expected$bands[[band]]$upper)
      }
    }
  }
})
