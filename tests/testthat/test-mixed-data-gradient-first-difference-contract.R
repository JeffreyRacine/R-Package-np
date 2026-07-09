library(np)

test_that("categorical-only regression gradients are first differences", {
  set.seed(20260709)
  n <- 36L
  dat <- data.frame(
    y = rnorm(n),
    xo = ordered(rep(letters[1:3], length.out = n)),
    xu = factor(rep(c("lo", "hi"), length.out = n), levels = c("lo", "hi"))
  )
  dat$y <- 0.4 * as.numeric(dat$xo) - 0.3 * (dat$xu == "hi") +
    rnorm(n, sd = 0.05)

  fit.ord <- np::npreg(y ~ xo, data = dat, gradients = TRUE)
  expect_equal(dim(fit.ord$grad), c(n, 1L))
  expect_true(all(is.finite(fit.ord$grad)))
  expect_equal(np::gradients(fit.ord, gradient.order = 1L), fit.ord$grad,
               tolerance = 0)
  expect_error(np::gradients(fit.ord, gradient.order = 2L),
               "first differences")

  fit.uno <- np::npreg(y ~ xu, data = dat, gradients = TRUE)
  expect_equal(dim(fit.uno$grad), c(n, 1L))
  expect_true(all(is.finite(fit.uno$grad)))

  fit.mix <- np::npreg(y ~ xo + xu, data = dat, gradients = TRUE)
  expect_equal(dim(fit.mix$grad), c(n, 2L))
  expect_true(all(is.finite(fit.mix$grad)))
})

test_that("categorical-only regression gradient plots are allowed", {
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

  fit.ord <- np::npreg(y ~ xo, data = dat)
  ord.eval <- data.frame(xo = ordered(levels(dat$xo), levels = levels(dat$xo)))
  ord.direct <- np::npreg(bws = fit.ord$bws, exdat = ord.eval,
                          gradients = TRUE)$grad[, 1L]
  ord.plot <- plot(fit.ord, gradients = TRUE, plot.behavior = "data")
  expect_equal(as.numeric(ord.plot[[1L]]$grad), as.numeric(ord.direct),
               tolerance = 1e-12)

  fit.uno <- np::npreg(y ~ xu, data = dat)
  uno.eval <- data.frame(xu = factor(levels(dat$xu), levels = levels(dat$xu)))
  uno.direct <- np::npreg(bws = fit.uno$bws, exdat = uno.eval,
                          gradients = TRUE)$grad[, 1L]
  uno.plot <- plot(fit.uno, gradients = TRUE, plot.behavior = "data")
  expect_equal(as.numeric(uno.plot[[1L]]$grad), as.numeric(uno.direct),
               tolerance = 1e-12)

  expect_no_error(plot(fit.ord, gradients = TRUE))
  expect_no_error(plot(fit.uno, gradients = TRUE))
})

test_that("categorical-only conditional gradient accessors validate first differences", {
  set.seed(20260709)
  n <- 32L
  dat <- data.frame(
    y = rnorm(n),
    xo = ordered(rep(letters[1:3], length.out = n))
  )
  dat$y <- 0.4 * as.numeric(dat$xo) + rnorm(n, sd = 0.05)

  dens <- np::npcdens(y ~ xo, data = dat, gradients = TRUE)
  expect_equal(np::gradients(dens, gradient.order = 1L), dens$congrad,
               tolerance = 0)
  expect_error(np::gradients(dens, gradient.order = 2L),
               "first differences")

  dist <- np::npcdist(y ~ xo, data = dat, gradients = TRUE)
  expect_equal(np::gradients(dist, gradient.order = 1L), dist$congrad,
               tolerance = 0)
  expect_error(np::gradients(dist, gradient.order = 2L),
               "first differences")
})

test_that("categorical-only nplsqreg gradients route through regression effects", {
  set.seed(20260709)
  n <- 32L
  dat <- data.frame(
    y = rnorm(n),
    xo = ordered(rep(letters[1:3], length.out = n))
  )
  dat$y <- 0.4 * as.numeric(dat$xo) + rnorm(n, sd = 0.05)

  fit <- np::nplsqreg(y ~ xo, data = dat, tau = 0.5, gradients = TRUE)
  expect_equal(dim(np::gradients(fit)), c(n, 1L))
  expect_true(all(is.finite(np::gradients(fit))))
})
