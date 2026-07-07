test_that("np.objective.cache controls npscoef continuous NN R optimizer caching", {
  old <- options(np.messages = FALSE, np.objective.cache = TRUE)
  on.exit(options(old), add = TRUE)

  run_bw <- function(bwtype, cache) {
    set.seed(20260523)
    n <- 80L
    z1 <- runif(n)
    z2 <- runif(n)
    x <- rnorm(n)
    y <- 1 + (0.5 + sin(2 * pi * z1)) * x +
      0.25 * cos(2 * pi * z2) + rnorm(n, sd = 0.2)
    options(np.objective.cache = cache)
    np::npscoefbw(
      xdat = data.frame(x = x),
      ydat = y,
      zdat = data.frame(z1 = z1, z2 = z2),
      regtype = "lc",
      bwtype = bwtype,
      ckertype = "epanechnikov",
      nmulti = 1L,
      optim.maxit = 35L,
      optim.maxattempts = 1L
    )
  }

  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    cached <- run_bw(bwtype, TRUE)
    uncached <- run_bw(bwtype, FALSE)

    expect_equal(cached$bw, uncached$bw)
    expect_equal(cached$fval, uncached$fval, tolerance = 0)
    expect_equal(cached$num.feval, uncached$num.feval)

    expect_equal(unname(cached$nn.cache[["enabled"]]), 1)
    expect_gt(unname(cached$nn.cache[["hits"]]), 0)
    expect_gte(as.numeric(cached$num.feval.fast[1L]),
               unname(cached$nn.cache[["hits"]]))
    expect_lte(as.numeric(cached$num.feval.fast[1L]),
               as.numeric(cached$num.feval[1L]))

    expect_equal(unname(uncached$nn.cache[["enabled"]]), 0)
    expect_equal(unname(uncached$nn.cache[["hits"]]), 0)
  }
})

test_that("npscoef R NN cache leaves fixed bandwidth searches unmarked", {
  old <- options(np.messages = FALSE, np.objective.cache = TRUE)
  on.exit(options(old), add = TRUE)

  set.seed(20260524)
  n <- 70L
  z1 <- runif(n)
  z2 <- runif(n)
  x <- rnorm(n)
  y <- 1 + (0.5 + z1) * x + rnorm(n, sd = 0.25)

  bw <- np::npscoefbw(
    xdat = data.frame(x = x),
    ydat = y,
    zdat = data.frame(z1 = z1, z2 = z2),
    regtype = "lc",
    bwtype = "fixed",
    ckertype = "epanechnikov",
    nmulti = 1L,
    optim.maxit = 25L,
    optim.maxattempts = 1L
  )

  expect_gte(unname(bw$nn.cache[["enabled"]]), 0)
  expect_gte(unname(bw$nn.cache[["hits"]]), 0)
})
