test_that("npcdens LC fixed no-error plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npcdens_lc_fixed_none_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(123)

  n <- 80L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rnorm(n))
  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    nmulti = 1L,
    regtype = "lc",
    bwtype = "fixed"
  )

  old <- plot(
    bw,
    xdat = x,
    ydat = y,
    plot.behavior = "data",
    plot.errors.method = "none",
    view = "fixed",
    neval = 9L,
    perspective = TRUE
  )
  candidate <- proto(bw, xdat = x, ydat = y, neval = 9L)

  expect_named(candidate, names(old))
  expect_s3_class(candidate$cd1, "condensity")
  expect_named(candidate$cd1, names(old$cd1))
  expect_equal(candidate$cd1$xeval, old$cd1$xeval)
  expect_equal(candidate$cd1$yeval, old$cd1$yeval)
  expect_equal(candidate$cd1$condens, old$cd1$condens)
  expect_equal(candidate$cd1$conderr, old$cd1$conderr)
  expect_equal(candidate$cd1$bias, old$cd1$bias)
  expect_identical(candidate$cd1$proper.requested, old$cd1$proper.requested)
  expect_identical(candidate$cd1$proper.applied, old$cd1$proper.applied)
})

test_that("npcdens plot prototype fails early outside its vertical slice", {
  proto <- getFromNamespace(".np_plot_proto_npcdens_lc_fixed_none_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(124)

  n <- 60L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rnorm(n))
  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    nmulti = 1L,
    regtype = "ll",
    bwtype = "fixed"
  )

  expect_error(
    proto(bw, xdat = x, ydat = y, neval = 9L),
    "regtype='lc'",
    fixed = TRUE
  )

  expect_error(
    proto(bw, neval = 9L),
    "explicit xdat and ydat",
    fixed = TRUE
  )
})
