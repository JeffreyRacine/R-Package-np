test_that("all-categorical conditional bandwidth search respects compression contract", {
  skip_on_cran()

  old_opts <- options(np.messages = FALSE, np.tree = FALSE,
                      np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  make_data <- function(seed, ordered_y = FALSE, ordered_x = FALSE) {
    set.seed(seed)
    n <- 300L
    y_raw <- sample(0:3, n, TRUE)
    x_raw <- sample(0:4, n, TRUE)
    z_raw <- sample(0:2, n, TRUE)
    data.frame(
      y = if (ordered_y) ordered(y_raw) else factor(y_raw),
      x = if (ordered_x) ordered(x_raw) else factor(x_raw),
      z = factor(z_raw)
    )
  }

  compare_bw <- function(bwfun, dat, strict_objective = TRUE, ...) {
    options(np.tree = FALSE, np.categorical.compress = FALSE)
    dense <- bwfun(y ~ x + z, data = dat, nmulti = 1, ...)

    options(np.tree = FALSE, np.categorical.compress = TRUE)
    profile <- bwfun(y ~ x + z, data = dat, nmulti = 1, ...)

    if (strict_objective)
      expect_equal(profile$fval, dense$fval, tolerance = 1e-8)
    else {
      expect_true(is.finite(profile$fval))
      expect_true(is.finite(dense$fval))
    }
    expect_true(all(is.finite(c(profile$ybw, profile$xbw))))
  }

  compare_bw(npcdensbw, make_data(20260621L), bwmethod = "cv.ml",
             strict_objective = FALSE)
  compare_bw(npcdensbw, make_data(20260622L, ordered_y = TRUE,
                                  ordered_x = TRUE),
             bwmethod = "cv.ls")
  compare_bw(npcdistbw, make_data(20260623L, ordered_y = TRUE,
                                  ordered_x = TRUE))
})

test_that("npcdist generated ordered grids preserve interval-labelled levels", {
  skip_on_cran()

  old_opts <- options(np.messages = FALSE, np.tree = FALSE,
                      np.categorical.compress = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260624L)
  n <- 180L
  x.raw <- rnorm(n)
  z.raw <- rnorm(n)
  s.raw <- x.raw + z.raw + rnorm(n, sd = 0.5)
  dat <- data.frame(
    y = ordered(cut(s.raw, quantile(s.raw, seq(0, 1, length.out = 5)),
                    include.lowest = TRUE)),
    x = ordered(cut(x.raw, quantile(x.raw, seq(0, 1, length.out = 5)),
                    include.lowest = TRUE)),
    z = ordered(cut(z.raw, quantile(z.raw, seq(0, 1, length.out = 5)),
                    include.lowest = TRUE))
  )

  expect_warning(
    bw <- npcdistbw(y ~ x + z, data = dat, nmulti = 1,
                    oxkertype = "liracine", oykertype = "liracine"),
    NA
  )
  expect_true(is.finite(bw$fval))
})
