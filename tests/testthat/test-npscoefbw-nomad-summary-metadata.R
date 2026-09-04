test_that("NOMAD smooth-coefficient metadata matches the manual constructor", {
  old <- options(np.messages = FALSE)
  on.exit(options(old))
  build <- getFromNamespace(".npscoefbw_build_scbandwidth", "np")
  normalize <- getFromNamespace(".npscoefbw_normalize_nomad_scbw", "np")
  x <- data.frame(x = (1:24) / 24)
  y <- sin(1:24)
  for (kind in c("continuous", "unordered", "ordered")) for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
    mixed <- kind != "continuous"
    z <- if (mixed) data.frame(z1 = x$x, z2 = factor(rep(c("a", "b"), 12), ordered = kind == "ordered")) else
      data.frame(z1 = x$x, z2 = x$x^2)
    bw <- c(if (type == "fixed") .3 else 16,
            if (mixed) .2 else if (type == "fixed") .4 else 18)
    reg.args <- list(regtype = "lp", degree = if (mixed) 0L else c(0L, 0L),
                     bwtype = type)
    bare <- build(x, y, z, bw, FALSE, reg.args)
    normalized <- normalize(bare, z)
    manual <- do.call(npscoefbw, c(list(xdat = x, ydat = y, zdat = z,
                                       bws = bw, bandwidth.compute = FALSE), reg.args))
    for (field in c("sumNum", "bandwidth", "sfactor")) {
      expect_identical(names(normalized[[field]]), names(manual[[field]]))
      # The regression-backed normalizer already retains inner variable names;
      # the manual smooth-coefficient constructor does not. Compare values.
      expect_identical(unname(unlist(normalized[[field]])),
                       unname(unlist(manual[[field]])))
    }
    expect_length(unlist(normalized$sumNum), ncol(z))
    expect_output(summary(normalized), "Var. Name:")
  }
  x2 <- data.frame(x1 = x$x, x2 = x$x^2)
  implicit <- normalize(build(x2, y, NULL, c(.3, .4), FALSE,
                              list(regtype = "lp", degree = c(0L, 0L))), x2)
  expect_identical(names(implicit$sumNum), "x")
  expect_length(unlist(implicit$sumNum), 2L)
  expect_true(all(is.finite(unlist(implicit$sumNum))))
  expect_output(summary(implicit), "Var. Name:")
})

test_that("fixed NOMAD publishes complete scalar and multivariate metadata", {
  skip_if_not_installed("crs")
  old <- options(np.messages = FALSE)
  on.exit(options(old))
  set.seed(902L)
  n <- 32L
  x <- data.frame(x = runif(n))
  z <- data.frame(z1 = runif(n), z2 = runif(n))
  y <- x$x * (1 + z$z1) + .1 * cos(1:n)
  for (p in 1:2) {
    bw <- npscoefbw(xdat = x, ydat = y, zdat = z[seq_len(p)],
                    nomad = TRUE, search.engine = "nomad",
                    degree.min = 0L, degree.max = 1L, nmulti = 1L,
                    nomad.opts = list(MAX_BB_EVAL = 8L))
    expect_length(unlist(bw$bandwidth), p)
    expect_length(unlist(bw$sfactor), p)
    expect_length(unlist(bw$sumNum), p)
    expect_true(all(is.finite(unlist(bw$sumNum))))
    expect_equal(unname(unlist(bw$bandwidth)), unname(bw$bw), tolerance = 0)
    expect_output(summary(bw), "Var. Name:")
    fit <- npscoef(bws = bw, txdat = x, tydat = y, tzdat = z[seq_len(p)])
    expect_true(all(is.finite(fitted(fit))))
    expect_true(all(is.finite(predict(fit, exdat = x[1:3, , drop = FALSE],
                                      ezdat = z[1:3, seq_len(p), drop = FALSE]))))
  }
})
