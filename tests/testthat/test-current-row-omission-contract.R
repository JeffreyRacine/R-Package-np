co_strip <- function(x) { attr(x, "na.action") <- NULL; x }

test_that("quantile evaluation and its plot owner ignore obsolete input history", {
  
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(9233)
  x <- data.frame(x = runif(36))
  y <- x$x + rnorm(36, sd = .4)
  eraw <- data.frame(x = c(.15, NA, .4, .6, NA, .8))
  marked <- stats::na.exclude(eraw); plain <- co_strip(marked)
  cbw <- npcdistbw(xdat = x, ydat = y, bws = c(.5,.5),
    bandwidth.compute = FALSE)
  lbw <- nplsqregbw(xdat = x, ydat = y, scale = rep(1,36),
    bw = .5, delta = .4, bandwidth.compute = FALSE, tau = c(.3,.7))
  plot.owner <- getFromNamespace(".np_plot_quantile_eval", "np")
  for (tau in list(.5, c(.3,.7))) {
    for (plot in c(FALSE, TRUE)) {
      fun <- if (plot) plot.owner else npqreg
      args <- list(bws = cbw, txdat = x, tydat = y, tau = tau, gradients = TRUE)
      args[[if (plot) "need.errors" else "se"]] <- TRUE
      a <- do.call(fun, c(args, list(exdat = marked)))
      b <- do.call(fun, c(args, list(exdat = plain)))
      for (extract in list(fitted, se, gradients)) {
        av <- extract(a); bv <- extract(b)
        expect_identical(dim(av), dim(bv))
        expect_equal(as.numeric(av), as.numeric(bv), tolerance = 2e-12)
      }
    }
  }
  a <- nplsqreg(lbw, exdat = marked, se = TRUE, gradients = TRUE)
  b <- nplsqreg(lbw, exdat = plain, se = TRUE, gradients = TRUE)
  for (extract in list(fitted, se, gradients)) {
    expect_identical(dim(extract(a)), dim(extract(b)))
    expect_equal(as.numeric(extract(a)), as.numeric(extract(b)), tolerance = 2e-12)
  }
})


test_that("current-row omissions do not recycle historical restoration indices", {
  helper <- getFromNamespace(".np_current_rows_omit", "np")
  d <- data.frame(x = c(1, NA, 3, 4, NA, 6),
    u = factor(c("a", "a", "b", "b", "a", "b")))
  rownames(d) <- letters[seq_len(nrow(d))]
  for (omit in list(stats::na.omit, stats::na.exclude)) {
    marked <- omit(d); saved <- marked
    expect_null(helper(marked))
    expect_identical(marked, saved)
    marked$x[3L] <- NA_real_
    plain <- co_strip(marked)
    expected <- attr(stats::na.omit(plain), "na.action")
    expect_identical(helper(marked), expected)
    expect_identical(helper(d), attr(stats::na.omit(d), "na.action"))
  }
  # Fresh composites already have the correct policy; do not edit their owners.
  expect_null(attr(data.frame(saved), "na.action"))
  mat <- matrix(c(1, NA, 3, 4, 5, NA), ncol = 1L)
  clean <- stats::na.omit(mat)
  expect_null(helper(clean))
  expect_null(attr(as.data.frame(clean), "na.action"))
  expect_identical(helper(mat), attr(stats::na.omit(mat), "na.action"))
})

test_that("kernel sums and hat owners use the supplied current rows", {
  
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(9231)
  raw <- data.frame(x = rnorm(30))
  raw$x[c(2L, 7L)] <- NA_real_
  marked <- stats::na.exclude(raw); plain <- co_strip(marked)
  eraw <- data.frame(x = seq(-.8, .8, length.out = 7L))
  eraw$x[c(2L,7L)] <- NA_real_
  emarked <- stats::na.omit(eraw); eplain <- co_strip(emarked)
  saved <- marked; esaved <- emarked
  y <- sin(plain$x) + rnorm(nrow(plain), sd = .2)
  rbw <- npregbw(xdat = plain, ydat = y, bws = .6, bandwidth.compute = FALSE)
  ub <- npudensbw(dat = plain, bws = .6, bandwidth.compute = FALSE)
  db <- npudistbw(dat = plain, bws = .6, bandwidth.compute = FALSE)
  for (tree in c(FALSE, TRUE)) {
    options(np.tree = tree)
    for (external in c(FALSE, TRUE)) {
      ea <- if (external) list(exdat = emarked) else list()
      eb <- if (external) list(exdat = eplain) else list()
      ka <- do.call(npksum, c(list(txdat = marked, bws = .6), ea))
      kb <- do.call(npksum, c(list(txdat = plain, bws = .6), eb))
      expect_identical(ka$ksum, kb$ksum)
      for (apply in c(FALSE, TRUE)) {
        extra <- if (apply) list(y = y, output = "apply") else list()
        a <- do.call(npreghat, c(list(bws = rbw, txdat = marked), ea, extra))
        b <- do.call(npreghat, c(list(bws = rbw, txdat = plain), eb, extra))
        expect_identical(dim(a), dim(b))
        expect_equal(as.numeric(a), as.numeric(b), tolerance = 2e-12)
        for (type in c("density", "distribution")) {
          fun <- if (type == "density") npudenshat else npudisthat
          bw <- if (type == "density") ub else db
          ua <- if (external) list(edat = emarked) else list()
          ubargs <- if (external) list(edat = eplain) else list()
          a <- do.call(fun, c(list(bws = bw, tdat = marked), ua, extra))
          b <- do.call(fun, c(list(bws = bw, tdat = plain), ubargs, extra))
          expect_identical(dim(a), dim(b))
          expect_equal(as.numeric(a), as.numeric(b), tolerance = 2e-12)
        }
      }
    }
  }
  # Matrix response and weight histories are independent omission owners.
  response <- matrix(y, ncol = 1L)
  attr(response, "na.action") <- attr(marked, "na.action")
  weight <- matrix(seq_len(nrow(plain))/nrow(plain), ncol = 1L)
  attr(weight, "na.action") <- attr(marked, "na.action")
  a <- npksum(txdat = marked, tydat = response, weights = weight, bws = .6)
  b <- npksum(txdat = plain, tydat = co_strip(response),
    weights = co_strip(weight), bws = .6)
  expect_identical(a$ksum, b$ksum)
  expect_identical(marked, saved); expect_identical(emarked, esaved)
})

test_that("formula NA restoration survives current-row normalization", {
  
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(9232)
  d <- data.frame(x = rnorm(24), y = rnorm(24))
  d$x[c(2L,9L)] <- NA_real_
  a <- npreg(y ~ x, data = d, bws = .6, na.action = na.exclude, se = TRUE)
  clean <- stats::na.omit(d)
  b <- npreg(y ~ x, data = clean, bws = .6, se = TRUE)
  expect_identical(which(is.na(fitted(a))), c(2L,9L))
  expect_equal(as.numeric(fitted(a)[-c(2L,9L)]), as.numeric(fitted(b)),
    tolerance = 2e-12)
  expect_length(predict(a, newdata = d, na.action = na.exclude), 24L)
})
