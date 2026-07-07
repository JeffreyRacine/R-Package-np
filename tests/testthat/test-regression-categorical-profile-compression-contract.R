test_that("all-categorical regression tree route preserves cv.ls and cv.aic bandwidth search", {
  skip_on_cran()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE, np.categorical.compress = FALSE)
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
    options(np.tree = FALSE, np.categorical.compress = FALSE)
    bw_dense <- npregbw(
      y ~ u1 + u2 + o1,
      data = dat,
      bwmethod = method,
      nmulti = 1,
      ukertype = "aitchisonaitken",
      okertype = "liracine"
    )

    options(np.tree = FALSE, np.categorical.compress = TRUE)
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
  old_opts <- options(np.messages = FALSE, np.tree = FALSE, np.categorical.compress = FALSE)
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
    options(np.tree = FALSE, np.categorical.compress = FALSE)
    bw_dense <- npregbw(
      y ~ u1 + o1,
      data = dat,
      bwmethod = method,
      nmulti = 1,
      ukertype = "liracine",
      okertype = "racineliyan"
    )

    options(np.tree = FALSE, np.categorical.compress = TRUE)
    bw_profile <- npregbw(
      y ~ u1 + o1,
      data = dat,
      bwmethod = method,
      nmulti = 1,
      ukertype = "liracine",
      okertype = "racineliyan"
    )

    if (identical(method, "cv.ls")) {
      expect_equal(bw_profile$fval, bw_dense$fval, tolerance = 1e-10)
      expect_lt(max(abs(as.numeric(bw_profile$bw) - as.numeric(bw_dense$bw))), 1e-8)
    } else {
      expect_equal(bw_profile$fval, bw_dense$fval, tolerance = 1e-4)
      expect_lt(max(abs(as.numeric(bw_profile$bw) - as.numeric(bw_dense$bw))), 1e-3)
      expect_equal(
        fitted(npreg(bws = bw_profile)),
        fitted(npreg(bws = bw_dense)),
        tolerance = 1e-3
      )
    }
  }
})

test_that("all-categorical regression RLY CV matches hat leave-one-out objective", {
  skip_on_cran()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE, np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260523)
  n <- 256L
  dat <- data.frame(
    y = rnorm(n),
    u1 = factor(sample(letters[1:2], n, TRUE)),
    o1 = ordered(sample(1:5, n, TRUE))
  )
  dat$y <- as.numeric(dat$u1) + cos(as.numeric(dat$o1)) + 0.1 * dat$y

  bw <- npregbw(
    y ~ u1 + o1,
    data = dat,
    bwmethod = "cv.ls",
    nmulti = 1,
    ukertype = "liracine",
    okertype = "racineliyan"
  )
  xdat <- dat[c("u1", "o1")]
  H <- npreghat(bws = bw, txdat = xdat, output = "matrix")
  fitted.full <- as.vector(H %*% dat$y)
  hii <- diag(H)
  fitted.loo <- (fitted.full - hii * dat$y) / (1 - hii)
  expect_equal(bw$fval, mean((dat$y - fitted.loo)^2), tolerance = 1e-8)
})

test_that("all-categorical regression tree route preserves fitted values and errors", {
  skip_on_cran()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE, np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260516)
  n <- 512L
  dat <- data.frame(
    y = rnorm(n),
    u1 = factor(rbinom(n, 1L, 0.5)),
    u2 = factor(sample(letters[1:3], n, TRUE)),
    u3 = factor(sample(LETTERS[1:2], n, TRUE)),
    o1 = ordered(sample(1:4, n, TRUE))
  )
  dat$y <- 1.5 * as.numeric(dat$u1) - 0.3 * as.numeric(dat$u2) +
    0.25 * as.numeric(dat$u3) + sin(as.numeric(dat$o1)) + 0.1 * dat$y

  options(np.tree = FALSE, np.categorical.compress = FALSE)
  bw <- npregbw(
    y ~ u1 + u2 + u3 + o1,
    data = dat,
    bwmethod = "cv.ls",
    nmulti = 1,
    ukertype = "aitchisonaitken",
    okertype = "liracine"
  )
  fit_dense <- npreg(bws = bw)

  options(np.tree = FALSE, np.categorical.compress = TRUE)
  fit_profile <- npreg(bws = bw)

  expect_equal(fitted(fit_profile), fitted(fit_dense), tolerance = 1e-8)
  expect_equal(fit_profile$merr, fit_dense$merr, tolerance = 1e-8)
  expect_equal(fit_profile$MSE, fit_dense$MSE, tolerance = 1e-10)
})

test_that("all-categorical regression tree route preserves native and predict evaluation", {
  skip_on_cran()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE, np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260517)
  n <- 640L
  dat <- data.frame(
    y = rnorm(n),
    u1 = factor(rbinom(n, 1L, 0.5)),
    u2 = factor(sample(letters[1:3], n, TRUE)),
    o1 = ordered(sample(1:4, n, TRUE))
  )
  dat$y <- as.numeric(dat$u1) + 0.5 * as.numeric(dat$u2) +
    sin(as.numeric(dat$o1)) + 0.1 * dat$y
  ex <- dat[c(seq_len(40L), seq_len(40L), sample(seq_len(n), 80L, TRUE)),
            c("u1", "u2", "o1"), drop = FALSE]

  options(np.tree = FALSE, np.categorical.compress = FALSE)
  bw <- npregbw(
    y ~ u1 + u2 + o1,
    data = dat,
    bwmethod = "cv.ls",
    nmulti = 1,
    ukertype = "aitchisonaitken",
    okertype = "liracine"
  )
  fit.base <- npreg(bws = bw)
  fit.dense <- npreg(bws = bw, exdat = ex)
  pred.dense <- predict(fit.base, newdata = ex, se.fit = TRUE)

  options(np.tree = FALSE, np.categorical.compress = TRUE)
  fit.profile <- npreg(bws = bw, exdat = ex)
  pred.profile <- predict(fit.base, newdata = ex, se.fit = TRUE)

  expect_equal(fitted(fit.profile), fitted(fit.dense), tolerance = 1e-8)
  expect_equal(fit.profile$merr, fit.dense$merr, tolerance = 1e-8)
  expect_equal(as.numeric(pred.profile$fit), as.numeric(pred.dense$fit),
               tolerance = 1e-8)
  expect_equal(as.numeric(pred.profile$se.fit), as.numeric(pred.dense$se.fit),
               tolerance = 1e-8)
})

test_that("all-categorical regression tree route preserves RLY evaluation", {
  skip_on_cran()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE, np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260518)
  n <- 256L
  dat <- data.frame(
    y = rnorm(n),
    u1 = factor(sample(letters[1:2], n, TRUE)),
    o1 = ordered(sample(1:5, n, TRUE))
  )
  dat$y <- as.numeric(dat$u1) + cos(as.numeric(dat$o1)) + 0.1 * dat$y
  ex <- dat[c(seq_len(30L), seq_len(30L)), c("u1", "o1"), drop = FALSE]

  options(np.tree = FALSE, np.categorical.compress = FALSE)
  bw <- npregbw(
    y ~ u1 + o1,
    data = dat,
    bwmethod = "cv.ls",
    nmulti = 1,
    ukertype = "liracine",
    okertype = "racineliyan"
  )
  fit.dense <- npreg(bws = bw, exdat = ex)

  options(np.tree = FALSE, np.categorical.compress = TRUE)
  fit.profile <- npreg(bws = bw, exdat = ex)

  expect_equal(fitted(fit.profile), fitted(fit.dense), tolerance = 1e-8)
  expect_equal(fit.profile$merr, fit.dense$merr, tolerance = 1e-8)
})

test_that("all-categorical regression tree route preserves npreghat apply output", {
  skip_on_cran()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE, np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260519)
  n <- 384L
  dat <- data.frame(
    y = rnorm(n),
    u1 = factor(rbinom(n, 1L, 0.5)),
    u2 = factor(sample(letters[1:3], n, TRUE)),
    o1 = ordered(sample(1:4, n, TRUE))
  )
  dat$y <- as.numeric(dat$u1) + 0.5 * as.numeric(dat$u2) +
    sin(as.numeric(dat$o1)) + 0.1 * dat$y
  xdat <- dat[c("u1", "u2", "o1")]
  ex <- xdat[c(seq_len(20L), seq_len(20L), sample(seq_len(n), 40L, TRUE)),
             , drop = FALSE]

  bw <- npregbw(
    y ~ u1 + u2 + o1,
    data = dat,
    bwmethod = "cv.ls",
    nmulti = 1,
    ukertype = "aitchisonaitken",
    okertype = "liracine"
  )

  H.train <- npreghat(bws = bw, txdat = xdat, output = "matrix")
  H.eval <- npreghat(bws = bw, txdat = xdat, exdat = ex, output = "matrix")

  options(np.tree = FALSE, np.categorical.compress = TRUE)
  apply.train <- npreghat(bws = bw, txdat = xdat, y = dat$y,
                          output = "apply")
  apply.eval <- npreghat(bws = bw, txdat = xdat, exdat = ex, y = dat$y,
                         output = "apply")

  expect_s3_class(H.train, "npreghat")
  expect_equal(apply.train, as.vector(H.train %*% dat$y), tolerance = 1e-8)
  expect_equal(apply.eval, as.vector(H.eval %*% dat$y), tolerance = 1e-8)
})

test_that("all-categorical regression tree route preserves deterministic plot payloads", {
  skip_on_cran()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE, np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260520)
  n <- 384L
  dat <- data.frame(
    y = rnorm(n),
    u1 = factor(rbinom(n, 1L, 0.5)),
    u2 = factor(sample(letters[1:3], n, TRUE)),
    o1 = ordered(sample(1:4, n, TRUE))
  )
  dat$y <- as.numeric(dat$u1) - 0.25 * as.numeric(dat$u2) +
    cos(as.numeric(dat$o1)) + 0.1 * dat$y

  bw <- npregbw(
    y ~ u1 + u2 + o1,
    data = dat,
    bwmethod = "cv.ls",
    nmulti = 1,
    ukertype = "aitchisonaitken",
    okertype = "liracine"
  )
  fit <- npreg(bws = bw)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  options(np.tree = FALSE, np.categorical.compress = FALSE)
  plot.dense <- plot(fit, plot.behavior = "data",
                     plot.errors.method = "none")

  options(np.tree = FALSE, np.categorical.compress = TRUE)
  plot.profile <- plot(fit, plot.behavior = "data",
                       plot.errors.method = "none")

  expect_equal(names(plot.profile), names(plot.dense))
  for (j in seq_along(plot.dense)) {
    expect_equal(plot.profile[[j]]$mean, plot.dense[[j]]$mean,
                 tolerance = 1e-8)
    expect_equal(plot.profile[[j]]$eval, plot.dense[[j]]$eval)
  }
})

test_that("all-categorical regression profile bootstrap matches hat algebra with fixed draws", {
  skip_on_cran()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE, np.categorical.compress = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260521)
  n <- 256L
  dat <- data.frame(
    y = rnorm(n),
    u1 = factor(rbinom(n, 1L, 0.5)),
    u2 = factor(sample(letters[1:3], n, TRUE)),
    o1 = ordered(sample(1:4, n, TRUE))
  )
  dat$y <- as.numeric(dat$u1) - 0.25 * as.numeric(dat$u2) +
    cos(as.numeric(dat$o1)) + 0.1 * dat$y
  xdat <- dat[c("u1", "u2", "o1")]
  ex <- xdat[c(seq_len(12L), seq_len(12L)), , drop = FALSE]

  bw <- npregbw(
    y ~ u1 + u2 + o1,
    data = dat,
    bwmethod = "cv.ls",
    nmulti = 1,
    ukertype = "aitchisonaitken",
    okertype = "liracine"
  )
  H.eval <- npreghat(bws = bw, txdat = xdat, exdat = ex, output = "matrix")
  H.train <- npreghat(bws = bw, txdat = xdat, output = "matrix")
  setup <- getFromNamespace(".np_regression_cat_profile_boot_setup", "np")(
    xdat = xdat, exdat = ex, ydat = dat$y, bws = bw
  )
  expect_false(is.null(setup))

  set.seed(551L)
  counts <- replicate(13L, tabulate(sample.int(n, n, TRUE), nbins = n))
  dense.inid <- getFromNamespace(".np_inid_lc_boot_from_hat", "np")(
    H = H.eval, ydat = dat$y, B = 13L, counts = counts
  )
  profile.inid <- getFromNamespace(".np_inid_boot_from_regression_cat_profile", "np")(
    setup = setup, B = 13L, counts = counts
  )
  expect_equal(profile.inid$t0, dense.inid$t0, tolerance = 1e-8)
  expect_equal(profile.inid$t, dense.inid$t, tolerance = 1e-8)

  fit.train <- as.vector(H.train %*% dat$y)
  set.seed(552L)
  dense.wild <- getFromNamespace(".np_plot_boot_from_hat_wild", "np")(
    H = H.eval, ydat = dat$y, fit.mean = fit.train,
    B = 13L, wild = "rademacher"
  )
  set.seed(552L)
  profile.wild <- getFromNamespace(".np_wild_boot_from_regression_cat_profile", "np")(
    setup = setup, B = 13L, wild = "rademacher"
  )
  expect_equal(profile.wild$t0, dense.wild$t0, tolerance = 1e-8)
  expect_equal(profile.wild$t, dense.wild$t, tolerance = 1e-8)
})

test_that("all-categorical regression tree route preserves bootstrap plot payloads", {
  skip_on_cran()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE, np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260521)
  n <- 512L
  dat <- data.frame(
    y = rnorm(n),
    u1 = factor(rbinom(n, 1L, 0.5)),
    u2 = factor(sample(letters[1:3], n, TRUE)),
    o1 = ordered(sample(1:4, n, TRUE))
  )
  dat$y <- as.numeric(dat$u1) - 0.25 * as.numeric(dat$u2) +
    cos(as.numeric(dat$o1)) + 0.1 * dat$y

  bw <- npregbw(
    y ~ u1 + u2 + o1,
    data = dat,
    bwmethod = "cv.ls",
    nmulti = 1,
    ukertype = "aitchisonaitken",
    okertype = "liracine"
  )
  fit <- npreg(bws = bw)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  compare_bootstrap_payloads <- function(method, seed) {
    common <- list(
      x = fit,
      plot.behavior = "data",
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = method,
      plot.errors.boot.num = 41L,
      plot.errors.boot.blocklen = 3L,
      plot.errors.type = "pointwise"
    )

    set.seed(seed)
    options(np.tree = FALSE, np.categorical.compress = FALSE)
    dense <- do.call(plot, common)

    set.seed(seed)
    options(np.tree = FALSE, np.categorical.compress = TRUE)
    profile <- do.call(plot, common)

    set.seed(seed)
    options(np.tree = FALSE, np.categorical.compress = TRUE)
    profile.repeat <- do.call(plot, common)

    expect_equal(names(profile), names(dense))
    expect_equal(profile.repeat, profile, tolerance = 1e-8)
    for (j in seq_along(profile)) {
      expect_equal(profile[[j]]$mean, dense[[j]]$mean, tolerance = 1e-8)
      expect_equal(profile[[j]]$eval, dense[[j]]$eval)
      expect_true(all(is.finite(profile[[j]]$merr)))
      expect_identical(is.finite(profile[[j]]$bias), is.finite(dense[[j]]$bias))
      expect_equal(profile[[j]]$bias, dense[[j]]$bias, tolerance = 1e-8)
    }
  }

  compare_bootstrap_payloads("inid", 771L)
  compare_bootstrap_payloads("fixed", 772L)
  compare_bootstrap_payloads("geom", 773L)
  compare_bootstrap_payloads("wild", 774L)
})

test_that("all-categorical regression profile bootstrap matches RLY hat algebra", {
  skip_on_cran()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE, np.categorical.compress = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260522)
  n <- 256L
  dat <- data.frame(
    y = rnorm(n),
    u1 = factor(sample(letters[1:2], n, TRUE)),
    o1 = ordered(sample(1:5, n, TRUE))
  )
  dat$y <- as.numeric(dat$u1) + cos(as.numeric(dat$o1)) + 0.1 * dat$y
  xdat <- dat[c("u1", "o1")]
  ex <- xdat[c(seq_len(12L), seq_len(12L)), , drop = FALSE]

  bw <- npregbw(
    y ~ u1 + o1,
    data = dat,
    bwmethod = "cv.ls",
    nmulti = 1,
    ukertype = "liracine",
    okertype = "racineliyan"
  )
  H.eval <- npreghat(bws = bw, txdat = xdat, exdat = ex, output = "matrix")
  setup <- getFromNamespace(".np_regression_cat_profile_boot_setup", "np")(
    xdat = xdat, exdat = ex, ydat = dat$y, bws = bw
  )
  expect_false(is.null(setup))

  set.seed(886L)
  counts <- replicate(13L, tabulate(sample.int(n, n, TRUE), nbins = n))
  dense.inid <- getFromNamespace(".np_inid_lc_boot_from_hat", "np")(
    H = H.eval, ydat = dat$y, B = 13L, counts = counts
  )
  profile.inid <- getFromNamespace(".np_inid_boot_from_regression_cat_profile", "np")(
    setup = setup, B = 13L, counts = counts
  )
  expect_equal(profile.inid$t0, dense.inid$t0, tolerance = 1e-8)
  expect_equal(profile.inid$t, dense.inid$t, tolerance = 1e-8)
})

test_that("all-categorical regression tree route preserves finite RLY bootstrap plots", {
  skip_on_cran()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE, np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260522)
  n <- 256L
  dat <- data.frame(
    y = rnorm(n),
    u1 = factor(sample(letters[1:2], n, TRUE)),
    o1 = ordered(sample(1:5, n, TRUE))
  )
  dat$y <- as.numeric(dat$u1) + cos(as.numeric(dat$o1)) + 0.1 * dat$y

  bw <- npregbw(
    y ~ u1 + o1,
    data = dat,
    bwmethod = "cv.ls",
    nmulti = 1,
    ukertype = "liracine",
    okertype = "racineliyan"
  )
  fit <- npreg(bws = bw)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  set.seed(885L)
  options(np.tree = FALSE, np.categorical.compress = FALSE)
  dense <- plot(fit,
                plot.behavior = "data",
                plot.errors.method = "bootstrap",
                plot.errors.boot.method = "inid",
                plot.errors.boot.num = 41L,
                plot.errors.type = "pointwise")

  set.seed(885L)
  options(np.tree = FALSE, np.categorical.compress = TRUE)
  tree <- plot(fit,
               plot.behavior = "data",
               plot.errors.method = "bootstrap",
               plot.errors.boot.method = "inid",
               plot.errors.boot.num = 41L,
               plot.errors.type = "pointwise")

  expect_equal(names(tree), names(dense))
  for (j in seq_along(dense)) {
    expect_equal(tree[[j]]$mean, dense[[j]]$mean, tolerance = 1e-8)
    expect_equal(tree[[j]]$eval, dense[[j]]$eval)
    expect_true(all(is.finite(tree[[j]]$merr)))
    expect_identical(is.finite(tree[[j]]$bias), is.finite(dense[[j]]$bias))
    expect_equal(tree[[j]]$bias, dense[[j]]$bias, tolerance = 1e-8)
  }
})

test_that("categorical compression option is isolated from continuous tree option", {
  skip_on_cran()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE,
                      np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260620L)
  n <- 256L
  dat_cat <- data.frame(
    y = rnorm(n),
    u1 = factor(rbinom(n, 1L, 0.5)),
    o1 = ordered(sample(1:4, n, TRUE))
  )
  dat_cat$y <- as.numeric(dat_cat$u1) + sin(as.numeric(dat_cat$o1)) +
    0.1 * dat_cat$y

  options(np.tree = FALSE, np.categorical.compress = FALSE)
  bw_dense <- npregbw(y ~ u1 + o1, data = dat_cat, nmulti = 1)
  fit_dense <- npreg(bws = bw_dense)

  options(np.tree = FALSE, np.categorical.compress = TRUE)
  bw_compress <- npregbw(y ~ u1 + o1, data = dat_cat, nmulti = 1)
  fit_compress <- npreg(bws = bw_compress)

  options(np.tree = TRUE, np.categorical.compress = FALSE)
  bw_tree <- npregbw(y ~ u1 + o1, data = dat_cat, nmulti = 1)
  fit_tree <- npreg(bws = bw_tree)

  expect_equal(bw_compress$fval, bw_dense$fval, tolerance = 1e-10)
  expect_equal(bw_tree$fval, bw_dense$fval, tolerance = 1e-10)
  expect_equal(fitted(fit_compress), fitted(fit_dense), tolerance = 1e-8)
  expect_equal(fitted(fit_tree), fitted(fit_dense), tolerance = 1e-8)

  dat_mixed <- data.frame(
    y = rnorm(n),
    x = runif(n),
    u1 = factor(rbinom(n, 1L, 0.5))
  )
  dat_mixed$y <- sin(2 * pi * dat_mixed$x) + as.numeric(dat_mixed$u1) +
    0.1 * dat_mixed$y

  options(np.tree = FALSE, np.categorical.compress = FALSE)
  bw_mixed_dense <- npregbw(y ~ x + u1, data = dat_mixed, nmulti = 1)
  fit_mixed_dense <- npreg(bws = bw_mixed_dense)

  options(np.tree = FALSE, np.categorical.compress = TRUE)
  bw_mixed_compress <- npregbw(y ~ x + u1, data = dat_mixed, nmulti = 1)
  fit_mixed_compress <- npreg(bws = bw_mixed_compress)

  expect_equal(bw_mixed_compress$fval, bw_mixed_dense$fval, tolerance = 1e-10)
  expect_equal(as.numeric(bw_mixed_compress$bw), as.numeric(bw_mixed_dense$bw),
               tolerance = 1e-8)
  expect_equal(fitted(fit_mixed_compress), fitted(fit_mixed_dense),
               tolerance = 1e-8)
})
