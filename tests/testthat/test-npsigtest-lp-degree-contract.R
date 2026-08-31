test_that("npsigtest LP degree validator is exact by tested coordinate", {
  validate_degree <- getFromNamespace(".np_npsig_validate_lp_degree", "npRmpi")
  xdat <- data.frame(
    x0 = c(0.1, 0.4, 0.9),
    x2 = c(0.2, 0.7, 0.3),
    z = factor(c("a", "b", "a"))
  )
  lc <- list(regtype = "lc", icon = c(TRUE, TRUE, FALSE), degree.engine = c(0L, 0L))
  lp <- list(regtype = "lp", icon = c(TRUE, TRUE, FALSE), degree.engine = c(0L, 2L))

  expect_true(validate_degree(lc, xdat, 1:2))
  expect_true(validate_degree(lp, xdat, 2L))
  expect_true(validate_degree(lp, xdat, 3L))
  expect_error(validate_degree(lp, xdat, 1L), "x0", fixed = TRUE)
  expect_error(validate_degree(lp, xdat, 1:2), "degree.min = 1", fixed = TRUE)

  malformed <- lp
  malformed$degree.engine <- 2L
  expect_error(
    validate_degree(malformed, xdat, 2L),
    "inconsistent local-polynomial degree metadata",
    fixed = TRUE
  )
})

test_that("npsigtest rejects tested LP degree zero before distributed work", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_options <- options(np.messages = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_options), add = TRUE)

  set.seed(280829)
  n <- 45L
  x0 <- runif(n)
  x2 <- runif(n)
  y <- x0 + x2^2 + rnorm(n, sd = 0.1)
  bw <- npregbw(
    xdat = data.frame(x0, x2),
    ydat = y,
    bws = c(0.3, 0.3),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = c(0L, 2L)
  )
  before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  expect_error(
    npsigtest(bw, index = 1L, pivot = FALSE, boot.num = 9L),
    "degree.min = 1",
    fixed = TRUE
  )
  expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), before)

  positive <- suppressWarnings(npsigtest(
    bw,
    index = 2L,
    pivot = FALSE,
    boot.num = 9L,
    random.seed = 42L,
    warn.glp.gradient = FALSE
  ))
  expect_s3_class(positive, "sigtest")
  expect_true(is.finite(positive$In))
  expect_true(all(is.finite(positive$In.bootstrap)))
})
