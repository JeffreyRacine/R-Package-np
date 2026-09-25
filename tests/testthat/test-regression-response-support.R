# This pure metadata test is intentionally small enough for the default lane.
test_that("regression support excludes response roles without losing X levels", {
  x <- list(iuno = c(FALSE, TRUE, FALSE), iord = c(TRUE, FALSE, FALSE),
            all.dlev = list(c(-2, 0, 4, 8), as.double(1:5), numeric()))
  for (ordered in c(FALSE, TRUE)) {
    y <- list(iuno = !ordered, iord = ordered,
              all.dlev = list(as.double(1:7)))
    b <- structure(list(dati = list(x = x, y = y)), class = "rbandwidth")
    expected <- list(as.double(1:5), c(-2, 0, 4, 8))
    expect_identical(.np_native_categorical_support(b), expected)
    class(b) <- c("response_support_subclass", "rbandwidth")
    expect_identical(.np_native_categorical_support(b), expected)
    b$dati$x <- list(iuno = FALSE, iord = FALSE, all.dlev = list(numeric()))
    expect_identical(.np_native_categorical_support(b), list())
    names(b$dati$x$all.dlev) <- "continuous"
    expect_identical(.np_native_categorical_support(b), list())
  }
})

test_that("conditional support retains Y then X and unconditional retains X", {
  x <- list(iuno = TRUE, iord = FALSE, all.dlev = list(as.double(1:4)))
  y <- list(iuno = FALSE, iord = TRUE, all.dlev = list(c(-3, 0, 8)))
  for (type in c("conbandwidth", "condbandwidth")) {
    b <- structure(list(dati = list(x = x, y = y)), class = type)
    expect_identical(.np_native_categorical_support(b),
                     list(c(-3, 0, 8), as.double(1:4)))
  }
  for (type in c("bandwidth", "dbandwidth")) {
    b <- structure(list(dati = list(x = x)), class = type)
    expect_identical(.np_native_categorical_support(b), list(as.double(1:4)))
  }
})

test_that("factor-response bandwidth searches and partially linear children run", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(42)
  n <- 48L
  x1 <- rnorm(n)
  x2 <- ordered(rbinom(n, 5, .3))
  z1 <- ordered(rbinom(n, 2, .3))
  z2 <- rnorm(n)
  y <- 1 + x1 + as.numeric(as.character(x2)) + as.numeric(as.character(z1)) +
    sin(z2) + rnorm(n)
  z <- data.frame(z1 = z1, z2 = z2)
  for (response in list(x2, factor(x2))) {
    encoded <- dlev(response)[as.integer(response)]
    set.seed(23)
    b <- npregbw(xdat = z, ydat = response, nmulti = 1L)
    set.seed(23)
    reference <- npregbw(xdat = z, ydat = encoded, nmulti = 1L)
    expect_equal(b$bw, reference$bw, tolerance = 2e-10)
    expect_equal(b$fval, reference$fval, tolerance = 2e-10)
    expect_equal(as.numeric(fitted(npreg(b))), as.numeric(fitted(npreg(reference))),
                 tolerance = 2e-10)
  }
  b <- npplregbw(y ~ x1 + x2 | z1 + z2, nmulti = 1L)
  expect_s3_class(b, "plbandwidth")
  expect_true(all(is.finite(fitted(npplreg(b)))))
})
