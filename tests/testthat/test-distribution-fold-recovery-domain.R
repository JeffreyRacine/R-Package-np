test_that("automatic fold recovery separates incumbent and effective domains", {
  schedule <- np:::.np_nn_ordinary_schedule(c(23, .3, 5), c(1L, 3L),
    caps = c(22, 22), incumbent.caps = 23)
  expect_identical(schedule, list(c(22, .3, 5), c(22, .3, 10),
                                 c(22, .3, 20), c(22, .3, 22)))
  expect_length(np:::.np_nn_ordinary_schedule(c(24, .3, 5), c(1L, 3L),
    caps = c(22, 22), incumbent.caps = 23), 0L)
  expect_length(np:::.np_nn_ordinary_schedule(23, 1L, caps = 22), 0L)
  seen <- list()
  z <- np:::.np_nn_find_raw_valid_start(c(23, .3), 1L, 22,
    raw.eval = function(x) { seen[[length(seen)+1L]] <<- x; .25 },
    incumbent.caps = 23)
  expect_identical(seen, list(c(22, .3)))
  expect_true(z$found)
  expect_identical(z$evaluations, 1L)
  expect_error(np:::.np_nn_ordinary_schedule(2, 1L, 4, incumbent.caps = 3),
               "invalid incumbent caps", fixed = TRUE)
})

test_that("distribution recovery respects the deleted-fit cap and raw objective", {
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = c(rep(0, 20), 1:4)); n <- nrow(x)
  oracle <- function(k) mean(vapply(seq_len(n), function(i) {
    d <- x$x[-i]
    f <- vapply(seq_along(d), function(j) {
      h <- sort(abs(d[-j] - d[j]))[k]
      if (is.na(h) || h <= 0) return(NA_real_)
      # Preserve the package's historical Gaussian-CDF coefficient. The
      # deleted-sample geometry/objective remains independently assembled.
      mean(pnorm(sqrt(2) * 0.7071067810 * (d[j] - d)/h))
    }, 0.)
    mean(((x$x[i] <= d) - f)^2)
  }, 0.))
  raw <- function(b) np:::npudistbw.dbandwidth(dat = x, bws = b,
    eval.only = TRUE, invalid.penalty = "dbmax", nmulti = 1L,
    do.full.integral = TRUE)$fval
  for (k in 20:23) {
    b <- npudistbw(x, bws = k, bwtype = "generalized_nn", bandwidth.compute = FALSE)
    if (k == 23) expect_identical(as.numeric(raw(b)), .Machine$double.xmax) else
      expect_equal(as.numeric(raw(b)), oracle(k), tolerance = 2e-12)
  }
  for (type in c("generalized_nn", "adaptive_nn")) {
    for (solver in c("powell", "mads", "mads+powell")) {
      set.seed(42)
      b <- npudistbw(x, bwtype = type, bwsolver = solver, nmulti = 1L,
        itmax = 20L, powell.remin = FALSE, do.full.integral = TRUE,
        nomad.opts = list(MAX_BB_EVAL = 30L))
      expect_true(b$bw >= 20 && b$bw <= 22)
      expect_true(is.finite(b$fval) && b$fval < .Machine$double.xmax)
      expect_equal(as.numeric(b$fval), as.numeric(raw(b)), tolerance = 2e-12)
      if (type == "generalized_nn")
        expect_equal(as.numeric(b$fval), oracle(b$bw), tolerance = 2e-12)
    }
  }
  b <- npudistbw(x, bws = 23, bwtype = "generalized_nn", bandwidth.compute = FALSE)
  options(np.extendednn = TRUE)
  expect_true(is.finite(raw(b)) && raw(b) < .Machine$double.xmax)
})
