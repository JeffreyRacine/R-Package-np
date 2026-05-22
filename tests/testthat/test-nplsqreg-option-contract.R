test_that("nplsqreg honors estimator bandwidth options separately from pilot options", {
  options(np.messages = FALSE)
  set.seed(20260521)
  dat <- data.frame(
    y = sin(seq(0, 2 * pi, length.out = 30L)) + rnorm(30L, sd = 0.1),
    x = seq(0, 1, length.out = 30L),
    u = factor(rep(letters[1:3], length.out = 30L)),
    o = ordered(rep(1:4, length.out = 30L))
  )
  scale0 <- rep(1, nrow(dat))

  fit.ll <- nplsqreg(y ~ x, data = dat, regtype = "ll",
                     scale = scale0, nmulti = 1L,
                     optim.control = list(maxit = 2L))
  expect_identical(fit.ll$bws$reg.bws$regtype, "ll")

  fit.lp <- nplsqreg(y ~ x, data = dat, regtype = "lp", degree = 1L,
                     scale = scale0, nmulti = 1L,
                     optim.control = list(maxit = 2L))
  expect_identical(fit.lp$bws$reg.bws$regtype, "lp")
  expect_identical(as.integer(fit.lp$bws$reg.bws$degree), 1L)

  fit.nn <- nplsqreg(bws = 4L, txdat = dat["x"], tydat = dat$y,
                     bwtype = "generalized_nn", scale = scale0,
                     bandwidth.compute = FALSE)
  expect_identical(fit.nn$bws$reg.bws$type, "generalized_nn")

  fit.ker <- nplsqreg(
    y ~ x + u + o,
    data = dat,
    ckertype = "epanechnikov",
    ckerorder = 4L,
    ukertype = "liracine",
    okertype = "wangvanryzin",
    scale = scale0,
    nmulti = 1L,
    optim.control = list(maxit = 2L)
  )
  expect_identical(fit.ker$bws$reg.bws$ckertype, "epanechnikov")
  expect_identical(as.integer(fit.ker$bws$reg.bws$ckerorder), 4L)
  expect_identical(fit.ker$bws$reg.bws$ukertype, "liracine")
  expect_identical(fit.ker$bws$reg.bws$okertype, "wangvanryzin")
})

test_that("nplsqreg vector tau applies estimator options to each full-search tau", {
  options(np.messages = FALSE)
  set.seed(20260521)
  dat <- data.frame(
    y = cos(seq(0, 2 * pi, length.out = 24L)) + rnorm(24L, sd = 0.1),
    x = seq(0, 1, length.out = 24L)
  )

  fit <- nplsqreg(
    y ~ x,
    data = dat,
    tau = c(0.25, 0.5),
    tau.search = "full",
    regtype = "ll",
    scale = rep(1, nrow(dat)),
    nmulti = 1L,
    optim.control = list(maxit = 2L)
  )

  expect_identical(fit$bws$tau.search, "full")
  expect_identical(fit$bws$fit.order, seq_along(fit$tau))
  expect_true(all(vapply(
    fit$bws$tau.bws,
    function(bw) identical(bw$reg.bws$regtype, "ll"),
    logical(1)
  )))
})
