test_that("lp entry points require explicit degree", {
  set.seed(1)
  n <- 24L
  x <- runif(n)
  z <- rnorm(n)
  y <- x + z + rnorm(n, sd = 0.1)

  calls <- list(
    npregbw = quote(npregbw(xdat = data.frame(x = x), ydat = y, regtype = "lp")),
    npcdensbw = quote(npcdensbw(xdat = data.frame(x = x), ydat = y, regtype = "lp")),
    npcdistbw = quote(npcdistbw(xdat = data.frame(x = x), ydat = y, regtype = "lp")),
    npindexbw = quote(npindexbw(xdat = data.frame(x = x), ydat = y, regtype = "lp")),
    npplregbw = quote(npplregbw(xdat = data.frame(x = x), ydat = y, zdat = data.frame(z = z), regtype = "lp")),
    npscoefbw = quote(npscoefbw(xdat = data.frame(x = x), ydat = y, zdat = data.frame(z = z), regtype = "lp")),
    npreg = quote(npreg(y ~ x, regtype = "lp")),
    npcdens = quote(npcdens(y ~ x, regtype = "lp")),
    npcdist = quote(npcdist(y ~ x, regtype = "lp")),
    npindex = quote(npindex(y ~ x, regtype = "lp")),
    npplreg = quote(npplreg(y ~ x | z, regtype = "lp")),
    npscoef = quote(npscoef(y ~ x | z, regtype = "lp"))
  )

  for (nm in names(calls)) {
    expect_error(
      suppressWarnings(eval(calls[[nm]])),
      regexp = "degree must be supplied explicitly when regtype='lp'",
      info = nm
    )
  }
})

test_that("lp degree length validation remains in force", {
  set.seed(1)
  n <- 20L
  x1 <- runif(n)
  x2 <- runif(n)
  y <- x1 + x2 + rnorm(n, sd = 0.1)

  expect_error(
    npregbw(xdat = data.frame(x1 = x1, x2 = x2), ydat = y, regtype = "lp", degree = 1L),
    regexp = "degree must have one entry per continuous predictor",
    fixed = FALSE
  )

  expect_error(
    npreg(y ~ x1 + x2, regtype = "lp", degree = 1L),
    regexp = "degree must have one entry per continuous predictor",
    fixed = FALSE
  )
})

test_that("automatic lp degree search still works without explicit degree", {
  set.seed(7)
  n <- 18L
  x1 <- runif(n)
  x2 <- runif(n)
  z <- rnorm(n)
  y <- x1 + x2 + z + rnorm(n, sd = 0.1)

  common.args <- list(
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "cell",
    degree.min = 1L,
    degree.max = 1L,
    bwtype = "fixed"
  )

  fits <- list(
    npregbw = do.call(npregbw, c(list(
      xdat = data.frame(x = x1),
      ydat = y,
      nmulti = 1L
    ), common.args)),
    npcdensbw = do.call(npcdensbw, c(list(
      xdat = data.frame(x = x1),
      ydat = y
    ), common.args)),
    npcdistbw = do.call(npcdistbw, c(list(
      xdat = data.frame(x = x1),
      ydat = y
    ), common.args)),
    npindexbw = do.call(npindexbw, c(list(
      xdat = data.frame(x1 = x1, x2 = x2),
      ydat = y,
      nmulti = 1L
    ), common.args)),
    npplregbw = do.call(npplregbw, c(list(
      xdat = data.frame(x = x1),
      ydat = y,
      zdat = data.frame(z = z),
      nmulti = 1L
    ), common.args)),
    npscoefbw = do.call(npscoefbw, c(list(
      xdat = data.frame(x = x1),
      ydat = y,
      zdat = data.frame(z = z),
      nmulti = 1L
    ), common.args))
  )

  expect_equal(length(as.integer(fits$npregbw$degree)), 1L)
  expect_equal(length(as.integer(fits$npcdensbw$degree)), 1L)
  expect_equal(length(as.integer(fits$npcdistbw$degree)), 1L)
  expect_equal(length(as.integer(fits$npindexbw$degree)), 1L)
  expect_equal(length(as.integer(fits$npplregbw$degree)), 1L)
  expect_equal(length(as.integer(fits$npscoefbw$degree)), 1L)

  for (fit in fits)
    expect_true(all(as.integer(fit$degree) >= 0L))
})
