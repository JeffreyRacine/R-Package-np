test_that("named automatic regression shares the positional training transaction", {
  skip_on_cran()
  spawn_mpi_slaves()
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(870L)
  d <- data.frame(x = rnorm(24), z = rnorm(24))
  d$y <- sin(d$x) + rnorm(24, sd = .2)
  state <- new.env(parent = emptyenv())
  state$n <- 0L
  counted <- function(x) { state$n <- state$n + 1L; x }
  f <- y ~ counted(x) + z
  set.seed(871L)
  positional <- npreg(f, data = d, nmulti = 1L, itmax = 1L,
                      powell.remin = FALSE, se = TRUE, gradients = TRUE)
  rng <- .Random.seed
  expect_identical(state$n, 1L)
  state$n <- 0L
  set.seed(871L)
  named <- npreg(formula = f, data = d, nmulti = 1L, itmax = 1L,
                powell.remin = FALSE, se = TRUE, gradients = TRUE)
  expect_identical(state$n, 1L)
  expect_identical(.Random.seed, rng)
  for (field in c("mean", "merr", "grad", "gerr"))
    expect_identical(named[[field]], positional[[field]])
  for (field in c("bw", "fval", "num.feval"))
    expect_identical(named$bws[[field]], positional$bws[[field]])
  state$n <- 0L
  refit <- npreg(bws = named$bws, se = TRUE, gradients = TRUE)
  expect_identical(state$n, 1L)
  expect_identical(refit$mean, named$mean)
})
