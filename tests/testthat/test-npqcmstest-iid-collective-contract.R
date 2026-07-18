test_that("npqcmstest IID index plan matches serial draws and RNG state", {
  make_plan <- getFromNamespace(".npRmpi_qcms_iid_index_plan", "npRmpi")
  n <- 17L
  B <- 13L

  set.seed(8711)
  reference <- matrix(NA_integer_, nrow = B, ncol = n)
  for (ii in seq_len(B))
    reference[ii, ] <- sample.int(n, replace = TRUE)
  reference.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  set.seed(8711)
  candidate <- make_plan(n, B)
  candidate.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  expect_identical(candidate, reference)
  expect_identical(candidate.seed, reference.seed)
})

test_that("npqcmstest IID route declares its collective helper", {
  route <- paste(deparse(getFromNamespace("npqcmstest", "npRmpi")),
                 collapse = "\n")

  expect_match(route, ".npRmpi_qcms_iid_index_plan", fixed = TRUE)
  expect_match(route, ".npRmpi_qcms_collective_iid_bootstrap", fixed = TRUE)
  expect_match(route, 'identical(boot.method, "iid")', fixed = TRUE)
})

test_that("npqcmstest IID bootstrap uses collective transport", {
  skip_on_cran()
  skip_if_not_installed("quantreg")

  trace.file <- tempfile("npqcmstest-iid-transport-", fileext = ".tsv")
  old.trace <- Sys.getenv("NP_RMPI_BOOTSTRAP_TRANSPORT_TRACE_FILE", unset = NA)
  Sys.setenv(NP_RMPI_BOOTSTRAP_TRANSPORT_TRACE_FILE = trace.file)
  on.exit({
    if (is.na(old.trace))
      Sys.unsetenv("NP_RMPI_BOOTSTRAP_TRANSPORT_TRACE_FILE")
    else
      Sys.setenv(NP_RMPI_BOOTSTRAP_TRANSPORT_TRACE_FILE = old.trace)
  }, add = TRUE)

  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(8712)
  n <- 36L
  x <- rnorm(n)
  y <- 1 + 0.8 * x + (0.4 + 0.1 * abs(x)) * rnorm(n)
  model <- quantreg::rq(y ~ x, tau = 0.5, model = TRUE)

  old.options <- options(np.messages = FALSE)
  on.exit(options(old.options), add = TRUE)
  out <- npqcmstest(
    model = model,
    xdat = x,
    ydat = y,
    distribution = "bootstrap",
    boot.method = "iid",
    boot.num = 9L,
    random.seed = 8713,
    nmulti = 1L
  )

  expect_s3_class(out, "cmstest")
  expect_length(out$Jn.bootstrap, 9L)
  trace <- readLines(trace.file, warn = FALSE)
  expect_true(any(grepl("what=npqcmstest", trace, fixed = TRUE)))
  expect_true(any(grepl("event=fanout.collective.start", trace, fixed = TRUE)))
  expect_true(any(grepl("event=fanout.collective.done", trace, fixed = TRUE)))
  expect_true(any(grepl("method=iid", trace, fixed = TRUE)))
})
