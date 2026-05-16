.npreg_cat_profile_test_env <- new.env(parent = emptyenv())

.ensure_npreg_cat_profile_pool <- function() {
  if (!isTRUE(.npreg_cat_profile_test_env$started)) {
    npRmpi.init(nslaves = 1L, quiet = TRUE)
    .npreg_cat_profile_test_env$started <- TRUE
    withr::defer({
      if (isTRUE(.npreg_cat_profile_test_env$started)) {
        try(npRmpi.quit(force = TRUE), silent = TRUE)
        .npreg_cat_profile_test_env$started <- FALSE
      }
    }, envir = testthat::teardown_env())
  }
}

test_that("all-categorical regression tree route preserves cv.ls and cv.aic bandwidth search", {
  skip_on_cran()
  .ensure_npreg_cat_profile_pool()
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
  .ensure_npreg_cat_profile_pool()
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

test_that("all-categorical regression tree route preserves fitted values and errors", {
  skip_on_cran()
  .ensure_npreg_cat_profile_pool()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
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

  options(np.tree = FALSE)
  bw <- npregbw(
    y ~ u1 + u2 + u3 + o1,
    data = dat,
    bwmethod = "cv.ls",
    nmulti = 1,
    ukertype = "aitchisonaitken",
    okertype = "liracine"
  )
  fit_dense <- npreg(bws = bw)

  options(np.tree = TRUE)
  fit_profile <- npreg(bws = bw)

  expect_equal(fitted(fit_profile), fitted(fit_dense), tolerance = 1e-8)
  expect_equal(fit_profile$merr, fit_dense$merr, tolerance = 1e-8)
  expect_equal(fit_profile$MSE, fit_dense$MSE, tolerance = 1e-10)
})

test_that("all-categorical regression tree route preserves native and predict evaluation", {
  skip_on_cran()
  .ensure_npreg_cat_profile_pool()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
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

  options(np.tree = FALSE)
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

  options(np.tree = TRUE)
  fit.profile <- npreg(bws = bw, exdat = ex)
  pred.profile <- predict(fit.base, newdata = ex, se.fit = TRUE)

  expect_equal(fitted(fit.profile), fitted(fit.dense), tolerance = 1e-8)
  expect_equal(fit.profile$merr, fit.dense$merr, tolerance = 1e-8)
  expect_equal(as.numeric(pred.profile$fit), as.numeric(pred.dense$fit),
               tolerance = 1e-8)
  expect_equal(as.numeric(pred.profile$se.fit), as.numeric(pred.dense$se.fit),
               tolerance = 1e-8)
})

test_that("all-categorical regression tree route leaves RLY evaluation on dense path", {
  skip_on_cran()
  .ensure_npreg_cat_profile_pool()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
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

  options(np.tree = FALSE)
  bw <- npregbw(
    y ~ u1 + o1,
    data = dat,
    bwmethod = "cv.ls",
    nmulti = 1,
    ukertype = "liracine",
    okertype = "racineliyan"
  )
  fit.dense <- npreg(bws = bw, exdat = ex)

  options(np.tree = TRUE)
  fit.profile <- npreg(bws = bw, exdat = ex)

  expect_equal(fitted(fit.profile), fitted(fit.dense), tolerance = 1e-8)
  expect_equal(fit.profile$merr, fit.dense$merr, tolerance = 1e-8)
})

test_that("all-categorical regression tree route preserves npreghat apply output", {
  skip_on_cran()
  .ensure_npreg_cat_profile_pool()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
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

  options(np.tree = TRUE)
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
  .ensure_npreg_cat_profile_pool()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
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

  options(np.tree = FALSE)
  plot.dense <- plot(fit, plot.behavior = "data",
                     plot.errors.method = "none")

  options(np.tree = TRUE)
  plot.profile <- plot(fit, plot.behavior = "data",
                       plot.errors.method = "none")

  expect_equal(names(plot.profile), names(plot.dense))
  for (j in seq_along(plot.dense)) {
    expect_equal(plot.profile[[j]]$mean, plot.dense[[j]]$mean,
                 tolerance = 1e-8)
    expect_equal(plot.profile[[j]]$eval, plot.dense[[j]]$eval)
  }
})
