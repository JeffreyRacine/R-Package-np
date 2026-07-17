library(npRmpi)

test_that("tree and categorical compression predicates are orthogonal", {
  old_opts <- options(np.messages = FALSE,
                      np.tree = FALSE,
                      np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  cases <- expand.grid(tree = c(FALSE, TRUE),
                       compress = c(FALSE, TRUE))
  observed <- lapply(seq_len(nrow(cases)), function(i) {
    options(np.tree = cases$tree[i],
            np.categorical.compress = cases$compress[i])
    c(
      categorical = npRmpi:::npUseCategoricalCompress(ncon = 0L, ncat = 2L),
      continuous_tree = npRmpi:::npUseContinuousTree(ncon = 2L),
      mixed_flag = npRmpi:::npUseKernelAccelerationFlag(ncon = 1L, ncat = 2L),
      allcat_flag = npRmpi:::npUseKernelAccelerationFlag(ncon = 0L, ncat = 2L)
    )
  })
  observed <- do.call(rbind, observed)

  expected <- rbind(
    c(FALSE, FALSE, FALSE, FALSE),
    c(FALSE, TRUE, TRUE, FALSE),
    c(TRUE, FALSE, FALSE, TRUE),
    c(TRUE, TRUE, TRUE, TRUE)
  )
  colnames(expected) <- colnames(observed)

  expect_identical(unname(observed), unname(expected))
})

test_that("np.tree auto mode is bounded-kernel continuous tree only", {
  old_opts <- options(np.messages = FALSE,
                      np.tree = "auto",
                      np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  expect_true(npRmpi:::npUseContinuousTree(ncon = 1L, ckertype = "epanechnikov"))
  expect_true(npRmpi:::npUseContinuousTree(ncon = 1L, ckertype = "uniform"))
  options(np.tree = c(mode = "auto"))
  expect_true(npRmpi:::npUseContinuousTree(ncon = 1L, ckertype = "epanechnikov"))
  options(np.tree = "auto")
  expect_false(npRmpi:::npUseContinuousTree(ncon = 1L, ckertype = "gaussian"))
  expect_false(npRmpi:::npUseContinuousTree(ncon = 1L))
  expect_false(npRmpi:::npUseContinuousTree(ncon = 0L, ckertype = "epanechnikov"))

  expect_false(npRmpi:::npUseCategoricalCompress(ncon = 0L, ncat = 2L))
  options(np.categorical.compress = TRUE)
  expect_true(npRmpi:::npUseCategoricalCompress(ncon = 0L, ncat = 2L))

  options(np.tree = TRUE)
  expect_true(npRmpi:::npUseContinuousTree(ncon = 1L, ckertype = "gaussian"))

  options(np.tree = FALSE)
  expect_false(npRmpi:::npUseContinuousTree(ncon = 1L, ckertype = "epanechnikov"))
})

test_that("np.tree alone does not enable all-categorical profile helpers", {
  old_opts <- options(np.messages = FALSE,
                      np.tree = FALSE,
                      np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260522)
  n <- 48L
  xdat <- data.frame(
    u1 = factor(sample(letters[1:3], n, TRUE)),
    o1 = ordered(sample(1:4, n, TRUE))
  )
  ydat <- as.numeric(xdat$u1) + as.numeric(xdat$o1) + rnorm(n, sd = 0.1)
  bws <- npRmpi:::rbandwidth(
    bw = c(0.2, 0.3),
    nobs = nrow(xdat),
    xdati = npRmpi:::untangle(xdat),
    ydati = npRmpi:::untangle(data.frame(ydat)),
    regtype = "lc",
    bandwidth.compute = FALSE
  )

  options(np.tree = TRUE, np.categorical.compress = FALSE)
  expect_null(npRmpi:::.np_regression_cat_profile_mean(bws, xdat, ydat))

  options(np.tree = FALSE, np.categorical.compress = TRUE)
  profile.mean <- npRmpi:::.np_regression_cat_profile_mean(bws, xdat, ydat)
  expect_type(profile.mean, "double")
  expect_length(profile.mean, n)

  options(np.tree = TRUE, np.categorical.compress = FALSE)
  expect_equal(
    npRmpi:::.npreg_fit_tree_code(bws, ncon = bws$ncon, ncat = bws$nuno + bws$nord),
    npRmpi:::DO_TREE_NO
  )

  options(np.tree = FALSE, np.categorical.compress = TRUE)
  expect_equal(
    npRmpi:::.npregbw_tree_code(bws, ncon = bws$ncon, ncat = bws$nuno + bws$nord),
    npRmpi:::DO_TREE_YES
  )
  expect_equal(
    npRmpi:::.npreg_fit_tree_code(bws, ncon = bws$ncon, ncat = bws$nuno + bws$nord),
    npRmpi:::DO_TREE_YES
  )
})

test_that("np.tree auto inspects active bandwidth-object kernels", {
  old_opts <- options(np.messages = FALSE,
                      np.tree = "auto",
                      np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260523)
  n <- 48L
  xdat <- data.frame(x1 = runif(n), x2 = runif(n))
  ydat <- xdat$x1 - xdat$x2 + rnorm(n, sd = 0.1)
  xdati <- npRmpi:::untangle(xdat)
  ydati <- npRmpi:::untangle(data.frame(ydat))

  bw.epan <- npRmpi:::rbandwidth(
    bw = c(0.35, 0.35),
    nobs = nrow(xdat),
    xdati = xdati,
    ydati = ydati,
    regtype = "lc",
    ckertype = "epanechnikov",
    bandwidth.compute = FALSE
  )
  bw.gauss <- npRmpi:::rbandwidth(
    bw = c(0.35, 0.35),
    nobs = nrow(xdat),
    xdati = xdati,
    ydati = ydati,
    regtype = "lc",
    ckertype = "gaussian",
    bandwidth.compute = FALSE
  )

  expect_equal(
    npRmpi:::.npreg_fit_tree_code(bw.epan, ncon = bw.epan$ncon, ncat = 0L),
    npRmpi:::DO_TREE_YES
  )
  expect_equal(
    npRmpi:::.npreg_fit_tree_code(bw.gauss, ncon = bw.gauss$ncon, ncat = 0L),
    npRmpi:::DO_TREE_NO
  )
})

test_that("continuous and mixed tree route helpers still honor np.tree", {
  old_opts <- options(np.messages = FALSE,
                      np.tree = TRUE,
                      np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260522)
  n <- 48L
  xdat <- data.frame(
    x1 = runif(n),
    u1 = factor(sample(letters[1:3], n, TRUE))
  )
  ydat <- sin(2 * pi * xdat$x1) + as.numeric(xdat$u1) + rnorm(n, sd = 0.1)
  bws <- npRmpi:::rbandwidth(
    bw = c(0.35, 0.2),
    nobs = nrow(xdat),
    xdati = npRmpi:::untangle(xdat),
    ydati = npRmpi:::untangle(data.frame(ydat)),
    regtype = "lc",
    bandwidth.compute = FALSE
  )

  expect_true(npRmpi:::npUseContinuousTree(ncon = bws$ncon))
  expect_false(npRmpi:::npUseCategoricalCompress(ncon = bws$ncon,
                                                 ncat = bws$nuno + bws$nord))
  expect_equal(
    npRmpi:::.npreg_fit_tree_code(bws, ncon = bws$ncon, ncat = bws$nuno + bws$nord),
    npRmpi:::DO_TREE_YES
  )
})
