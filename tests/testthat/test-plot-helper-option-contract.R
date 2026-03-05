test_that("semihat regbw args forward index LP/kernel options with bound collapse", {
  set.seed(8311)
  n <- 64
  xdat <- data.frame(x1 = runif(n), x2 = runif(n))
  ydat <- rnorm(n)

  bw <- npindexbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.2, 0.2, 0.3),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = 2L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    ckertype = "epanechnikov",
    ckerorder = 2L,
    ckerbound = "fixed",
    ckerlb = c(0, 0),
    ckerub = c(1, 1)
  )
  idx.train <- data.frame(index = as.vector(toMatrix(xdat) %*% bw$beta))

  make.args <- getFromNamespace(".np_semihat_make_regbw_args", "np")
  args <- make.args(
    source = bw,
    xdat = idx.train,
    ydat = rep.int(0.0, nrow(idx.train)),
    bw = bw$bw
  )

  expect_identical(args$regtype, as.character(bw$regtype))
  expect_identical(args$basis, bw$basis)
  expect_equal(args$degree, bw$degree)
  expect_identical(isTRUE(args$bernstein.basis), isTRUE(bw$bernstein.basis))
  expect_identical(args$bwtype, bw$type)
  expect_identical(args$ckertype, bw$ckertype)
  expect_identical(args$ckerorder, bw$ckerorder)
  expect_identical(args$ckerbound, bw$ckerbound)
  expect_equal(as.double(args$ckerlb), 0)
  expect_equal(as.double(args$ckerub), 1)
  expect_identical(args$bandwidth.compute, FALSE)
})

test_that("semihat regbw args forward smooth-coef LP/kernel/scaling options", {
  set.seed(8312)
  n <- 58
  xdat <- data.frame(x = runif(n))
  zdat <- data.frame(z = runif(n))
  ydat <- rnorm(n)

  bw <- npscoefbw(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    bws = 0.25,
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = 2L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    bwscaling = TRUE,
    ckertype = "epanechnikov",
    ckerorder = 2L,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1
  )

  make.args <- getFromNamespace(".np_semihat_make_regbw_args", "np")
  args <- make.args(
    source = bw,
    xdat = zdat,
    ydat = rep.int(0.0, nrow(zdat)),
    bw = bw$bw
  )

  expect_identical(args$regtype, as.character(bw$regtype))
  expect_identical(args$basis, bw$basis)
  expect_equal(args$degree, bw$degree)
  expect_identical(isTRUE(args$bernstein.basis), isTRUE(bw$bernstein.basis))
  expect_identical(args$bwscaling, isTRUE(bw$scaling))
  expect_identical(args$bwtype, bw$type)
  expect_identical(args$ckertype, bw$ckertype)
  expect_identical(args$ckerorder, bw$ckerorder)
  expect_identical(args$ckerbound, bw$ckerbound)
  expect_equal(as.double(args$ckerlb), 0)
  expect_equal(as.double(args$ckerub), 1)
  expect_identical(args$bandwidth.compute, FALSE)

  if (!is.null(bw$method) && bw$method %in% c("cv.ls", "cv.aic")) {
    expect_identical(args$bwmethod, bw$method)
  }
})
