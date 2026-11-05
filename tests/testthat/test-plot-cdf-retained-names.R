test_that("corrected CDF plots preserve retained non-syntactic column names", {
  if (TRUE) {
    skip_on_cran()
    if (!spawn_mpi_slaves(as.integer(Sys.getenv("N5_SLAVES", "1"))))
      stop("required MPI pool unavailable")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(4101)
  d <- data.frame(z = rnorm(60))
  special <- d; names(special) <- "a b"
  bw <- npudistbw(dat = special, bws = .4, bandwidth.compute = FALSE)
  fit <- npudist(bws = bw, tdat = special)
  reference <- npudist(bws = npudistbw(dat = d, bws = .4, bandwidth.compute = FALSE), tdat = d)
  formula.bw <- npudistbw(~ `a b`, data = special, bws = .4, bandwidth.compute = FALSE)
  frozen <- serialize(fit, NULL)
  full <- identical(Sys.getenv("N5_FULL"), "true")
  for (method in if (full) c("asymptotic", "bootstrap") else "asymptotic") {
    args <- list(errors = method, center = "bias-corrected", band = "pointwise",
                 neval = 6L, output = "data")
    if (method == "bootstrap") args <- c(args, list(B = 999L, bootstrap = "inid", random_seed = 13L))
    expected <- do.call(plot, c(list(reference), args))
    for (object in list(bw, fit, formula.bw)) {
      got <- do.call(plot, c(list(object), args))
      expect_identical(got$panels[[1L]]$variables$name, "a b")
      expect_identical(got$panels[[1L]]$values, expected$panels[[1L]]$values)
      expect_identical(got$panels[[1L]]$intervals, expected$panels[[1L]]$intervals)
    }
  }
  explicit <- plot(bw, xdat = special, errors = "asymptotic", center = "bias-corrected",
                   band = "pointwise", neval = 6L, output = "data")
  retained <- plot(bw, errors = "asymptotic", center = "bias-corrected",
                   band = "pointwise", neval = 6L, output = "data")
  expect_identical(explicit$panels, retained$panels)
  expect_identical(serialize(fit, NULL), frozen)
  expect_error(plot(bw, xdat = d, errors = "asymptotic", center = "bias-corrected",
                    band = "pointwise", output = "data"), "invalid CDF correction training frame")
  bad <- data.frame(other = special[[1L]], second = rnorm(60))
  expect_error(plot(bw, xdat = bad, errors = "asymptotic", center = "bias-corrected",
                    band = "pointwise", output = "data"))
})
