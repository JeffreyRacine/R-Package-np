test_that("ordinary bootstrap upper-tail p-values are conservative and tie-inclusive", {
  pvalue <- getFromNamespace(".np_bootstrap_upper_tail_pvalue", "npRmpi")

  expect_identical(pvalue(rep(0, 9L), 1), 0)
  expect_identical(pvalue(rep(1, 9L), 1), 1)
  expect_identical(pvalue(c(0, 1, 1, 2), 1), 3 / 4)
  expect_error(pvalue(numeric(), 1), "at least one replication")
  expect_error(pvalue(1:9, c(1, 2)), "must be scalar")
})

test_that("categorical equality reports conservative tied bootstrap p-values", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- factor(rep(1:2, each = 6L), levels = 1:2)
  set.seed(221)
  before <- .Random.seed
  result <- npunitest(x, x, bw.x = 0.5, bw.y = 0.5, B = 9,
                     ukertype = "aitchisonaitken")
  expect_identical(.Random.seed, before)
  expect_identical(result$Srho, 0)
  expect_true(all(result$Srho.bootstrap == 0))
  expect_identical(result$P, 1)
})

test_that("all public bootstrap-inference p-value owners use the shared policy", {
  files <- c(
    "np.cmstest.R", "np.qcmstest.R", "np.deptest.R", "np.deneqtest.R",
    "np.unitest.R", "np.sdeptest.R", "np.symtest.R", "np.sigtest.R"
  )
  for (file in files) {
    source <- paste(readLines(test_path("..", "..", "R", file), warn = FALSE),
                    collapse = "\n")
    expect_match(source, ".np_bootstrap_upper_tail_pvalue", fixed = TRUE,
                 info = file)
  }
})
