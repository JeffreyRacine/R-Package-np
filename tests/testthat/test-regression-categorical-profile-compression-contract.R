test_that("all-categorical regression tree route preserves cv.ls and cv.aic bandwidth search", {
  skip_on_cran()
  npRmpi.init(nslaves = 1L, quiet = TRUE)
  on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(42)
  n <- 384L
  dat <- data.frame(
    y = rnorm(n),
    u1 = factor(rbinom(n, 1L, 0.5)),
    u2 = factor(sample(letters[1:3], n, TRUE)),
    o1 = ordered(sample(1:4, n, TRUE))
  )
  dat$y <- as.numeric(dat$u1) + 0.5 * as.numeric(dat$u2) +
    sin(as.numeric(dat$o1)) + 0.1 * dat$y

  for (method in c("cv.ls", "cv.aic")) {
    options(np.tree = FALSE)
    bw_dense <- npregbw(
      y ~ u1 + u2 + o1,
      data = dat,
      bwmethod = method,
      nmulti = 1,
      ukertype = "aitchisonaitken",
      okertype = "liracine"
    )

    options(np.tree = TRUE)
    bw_profile <- npregbw(
      y ~ u1 + u2 + o1,
      data = dat,
      bwmethod = method,
      nmulti = 1,
      ukertype = "aitchisonaitken",
      okertype = "liracine"
    )

    expect_equal(bw_profile$fval, bw_dense$fval, tolerance = 1e-10)
    expect_lt(max(abs(as.numeric(bw_profile$bw) - as.numeric(bw_dense$bw))), 1e-7)
    expect_equal(
      fitted(npreg(bws = bw_profile)),
      fitted(npreg(bws = bw_dense)),
      tolerance = 1e-8
    )
  }
})

test_that("all-categorical regression tree route preserves ordered Racine-Li-Yan route", {
  skip_on_cran()
  npRmpi.init(nslaves = 1L, quiet = TRUE)
  on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(314)
  n <- 256L
  dat <- data.frame(
    y = rnorm(n),
    u1 = factor(sample(letters[1:2], n, TRUE)),
    o1 = ordered(sample(1:5, n, TRUE))
  )
  dat$y <- as.numeric(dat$u1) + cos(as.numeric(dat$o1)) + 0.1 * dat$y

  for (method in c("cv.ls", "cv.aic")) {
    options(np.tree = FALSE)
    bw_dense <- npregbw(
      y ~ u1 + o1,
      data = dat,
      bwmethod = method,
      nmulti = 1,
      ukertype = "liracine",
      okertype = "racineliyan"
    )

    options(np.tree = TRUE)
    bw_profile <- npregbw(
      y ~ u1 + o1,
      data = dat,
      bwmethod = method,
      nmulti = 1,
      ukertype = "liracine",
      okertype = "racineliyan"
    )

    expect_equal(bw_profile$fval, bw_dense$fval, tolerance = 1e-10)
    expect_lt(max(abs(as.numeric(bw_profile$bw) - as.numeric(bw_dense$bw))), 1e-8)
  }
})
