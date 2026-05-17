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
