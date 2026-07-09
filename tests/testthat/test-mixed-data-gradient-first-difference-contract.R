library(npRmpi)

test_that("categorical-only regression gradients are first differences under autodispatch", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260709)
  n <- 36L
  dat <- data.frame(
    y = rnorm(n),
    xo = ordered(rep(letters[1:3], length.out = n)),
    xu = factor(rep(c("lo", "hi"), length.out = n), levels = c("lo", "hi"))
  )
  dat$y <- 0.4 * as.numeric(dat$xo) - 0.3 * (dat$xu == "hi") +
    rnorm(n, sd = 0.05)

  fit.ord <- npRmpi::npreg(y ~ xo, data = dat, gradients = TRUE)
  expect_equal(dim(fit.ord$grad), c(n, 1L))
  expect_true(all(is.finite(fit.ord$grad)))
  expect_equal(npRmpi::gradients(fit.ord, gradient.order = 1L),
               fit.ord$grad, tolerance = 0)
  expect_error(npRmpi::gradients(fit.ord, gradient.order = 2L),
               "first differences")

  fit.uno <- npRmpi::npreg(y ~ xu, data = dat, gradients = TRUE)
  expect_equal(dim(fit.uno$grad), c(n, 1L))
  expect_true(all(is.finite(fit.uno$grad)))

  fit.mix <- npRmpi::npreg(y ~ xo + xu, data = dat, gradients = TRUE)
  expect_equal(dim(fit.mix$grad), c(n, 2L))
  expect_true(all(is.finite(fit.mix$grad)))
})

test_that("categorical-only regression gradient plots are allowed under autodispatch", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260709)
  n <- 24L
  dat <- data.frame(
    y = rnorm(n),
    xo = ordered(rep(letters[1:3], length.out = n)),
    xu = factor(rep(c("lo", "hi"), length.out = n), levels = c("lo", "hi"))
  )
  dat$y <- 0.4 * as.numeric(dat$xo) - 0.3 * (dat$xu == "hi") +
    rnorm(n, sd = 0.05)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_no_error(plot(npRmpi::npreg(y ~ xo, data = dat), gradients = TRUE))
  expect_no_error(plot(npRmpi::npreg(y ~ xu, data = dat), gradients = TRUE))
})

test_that("categorical-only conditional gradient accessors validate first differences under autodispatch", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260709)
  n <- 32L
  dat <- data.frame(
    y = rnorm(n),
    xo = ordered(rep(letters[1:3], length.out = n))
  )
  dat$y <- 0.4 * as.numeric(dat$xo) + rnorm(n, sd = 0.05)

  dens <- npRmpi::npcdens(y ~ xo, data = dat, gradients = TRUE)
  expect_equal(npRmpi::gradients(dens, gradient.order = 1L), dens$congrad,
               tolerance = 0)
  expect_error(npRmpi::gradients(dens, gradient.order = 2L),
               "first differences")

  dist <- npRmpi::npcdist(y ~ xo, data = dat, gradients = TRUE)
  expect_equal(npRmpi::gradients(dist, gradient.order = 1L), dist$congrad,
               tolerance = 0)
  expect_error(npRmpi::gradients(dist, gradient.order = 2L),
               "first differences")
})

test_that("categorical-only nplsqreg gradients route through regression effects under autodispatch", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260709)
  n <- 32L
  dat <- data.frame(
    y = rnorm(n),
    xo = ordered(rep(letters[1:3], length.out = n))
  )
  dat$y <- 0.4 * as.numeric(dat$xo) + rnorm(n, sd = 0.05)

  fit <- npRmpi::nplsqreg(y ~ xo, data = dat, tau = 0.5, gradients = TRUE)
  expect_equal(dim(npRmpi::gradients(fit)), c(n, 1L))
  expect_true(all(is.finite(npRmpi::gradients(fit))))
})
