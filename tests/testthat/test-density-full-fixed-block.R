test_that("corrected density refuses a structurally degenerate full fixed block", {
  if (TRUE) {
    skip_on_cran()
    if (!spawn_mpi_slaves(as.integer(Sys.getenv("N5_SLAVES", "1")))) stop("MPI unavailable")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(61); n <- 80L; d <- data.frame(x = rnorm(n))
  f <- npudens(bws = npudensbw(dat = d, bws = .5, bandwidth.compute = FALSE), tdat = d)
  full <- identical(Sys.getenv("N5_FULL"), "true")
  diagnostics <- identical(Sys.getenv("N5_CONTROLS_ONLY"), "true")
  calibrations <- if (full) c("symmetric-sd", "asymmetric-sd") else "symmetric-sd"
  bands <- if (full) c("pointwise", "bonferroni", "simultaneous", "all") else "pointwise"
  if (!diagnostics) for (cal in calibrations) for (band in bands)
    expect_error(plot(f, errors = "bootstrap", bootstrap = "fixed", center = "bias-corrected",
      band = band, B = 999L, neval = 5L,
      boot.control = np_boot_control(calibration = cal, blocklen = n), output = "data"),
      "fixed block length equal to sample size")
  ns <- asNamespace("npRmpi")
  request <- get(".np_udens_b2_boot_request", ns)
  counts <- get(".np_block_counts_drawer", ns)(n, 9L, n, "fixed")(1L, 9L)
  expect_true(all(counts == 1)) # independent structural count witness
  seed <- .Random.seed
  if (!diagnostics) expect_error(request("fixed", n, n), "no resampling variation")
  expect_identical(.Random.seed, seed)
  expect_equal(request("fixed", n - 1L, n)$blocklen, n - 1L)
  expect_equal(request("geom", n, n)$blocklen, n)
  payloads <- list()
  if (full) {
    second <- data.frame(x = seq(-2, 2, length.out = n)^3)
    g <- npudens(bws = npudensbw(dat = second, bws = .7, bandwidth.compute = FALSE), tdat = second)
    if (!diagnostics) expect_error(plot(g, errors="bootstrap",bootstrap="fixed",center="bias-corrected",
      band="pointwise",B=49L,neval=5L,boot.control=np_boot_control(blocklen=n),output="data"),
      "no resampling variation")
    for (cal in calibrations) for (cfg in list(c("inid", 1), c("fixed", 1), c("fixed", n-1L), c("geom", n))) {
      seed <- .Random.seed
      a <- plot(f, errors="bootstrap",bootstrap=cfg[1],center="bias-corrected",band="all",B=999L,neval=5L,
        boot.control=np_boot_control(calibration=cal,blocklen=as.numeric(cfg[2])),output="data")
      expect_true(all(is.finite(as.data.frame(a, band="pointwise")$std_error)))
      expect_identical(.Random.seed, seed)
      payloads[[paste(cal,cfg[1],cfg[2])]] <- a
    }
  }
  if (nzchar(Sys.getenv("N5_SNAPSHOT"))) saveRDS(payloads, Sys.getenv("N5_SNAPSHOT"))
})
