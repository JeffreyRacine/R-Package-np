test_that("partially linear bootstrap compose chunk matches scalar loop", {
  finish_chunk <- getFromNamespace(".np_inid_boot_from_plreg_finish_chunk", "npRmpi")

  x.train.num <- matrix(c(1, 2, 3, 4), nrow = 4L, ncol = 1L)
  x.eval.num <- matrix(c(1.5, 2.5), nrow = 2L, ncol = 1L)
  y.num <- c(2, 3, 4, 5)
  y.train.t <- matrix(
    c(
      0.2, 0.3, 0.4, 0.5,
      0.1, 0.2, 0.3, 0.4
    ),
    nrow = 2L,
    byrow = TRUE
  )
  y.eval.t <- matrix(
    c(
      1.1, 1.2,
      1.3, 1.4
    ),
    nrow = 2L,
    byrow = TRUE
  )
  x.train.t.list <- list(matrix(
    c(
      0.05, 0.10, 0.15, 0.20,
      0.10, 0.15, 0.20, 0.25
    ),
    nrow = 2L,
    byrow = TRUE
  ))
  x.eval.t.list <- list(matrix(
    c(
      0.20, 0.25,
      0.30, 0.35
    ),
    nrow = 2L,
    byrow = TRUE
  ))
  counts.mat <- matrix(
    c(
      1, 0,
      1, 2,
      1, 1,
      1, 1
    ),
    nrow = 4L,
    ncol = 2L
  )
  storage.mode(counts.mat) <- "double"

  out <- finish_chunk(
    rows = 1:2,
    x.train.num = x.train.num,
    x.eval.num = x.eval.num,
    y.num = y.num,
    y.train.t = y.train.t,
    y.eval.t = y.eval.t,
    x.train.t.list = x.train.t.list,
    x.eval.t.list = x.eval.t.list,
    counts.mat = counts.mat
  )

  expected <- matrix(NA_real_, nrow = 2L, ncol = 2L)
  for (b in 1:2) {
    xres.train <- x.train.num[, 1L] - x.train.t.list[[1L]][b, ]
    xres.eval <- x.eval.num[, 1L] - x.eval.t.list[[1L]][b, ]
    yres <- y.num - y.train.t[b, ]
    beta <- sum(counts.mat[, b] * xres.train * yres) /
      sum(counts.mat[, b] * xres.train * xres.train)
    expected[b, ] <- y.eval.t[b, ] + xres.eval * beta
  }

  expect_equal(out, expected, tolerance = 1e-12)
})
