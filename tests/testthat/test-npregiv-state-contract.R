test_that("npregiv ridging loops avoid superassignment", {
  src_path <- testthat::test_path("..", "..", "R", "npregiv.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")
  lines <- grep("<<-", strsplit(src, "\n", fixed = TRUE)[[1L]], value = TRUE, fixed = TRUE)
  expect_false(any(grepl(
    "ridge|doridge|coef\\.mat|mean\\.loo|ghat|WzkWz",
    lines
  )))
})

test_that("npregiv uses shared seed enter helper", {
  src_path <- testthat::test_path("..", "..", "R", "npregiv.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")
  expect_true(grepl("\\.np_seed_enter\\(random\\.seed\\)", src))
  expect_false(grepl("exists\\(\"\\.Random\\.seed\"", src))
  expect_true(grepl("\\.np_seed_exit\\(seed\\.state, remove_if_absent = TRUE\\)", src))
})

test_that("npregiv Tikhonov refinement preserves the T-star-r right-hand side", {
  src_path <- testthat::test_path("..", "..", "R", "npregiv.R")
  skip_if_not(file.exists(src_path), "source R files unavailable in installed test context")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")
  expect_match(src,
               "T\\.star\\.r <- .*glpreg\\(tydat=E\\.y\\.w",
               perl = TRUE)
  expect_match(src, "Cr\\.r = T\\.star\\.r", perl = TRUE)
  expect_false(any(grepl("^\\s*E\\.E\\.phi\\.w\\.z <-",
                         readLines(src_path, warn = FALSE), perl = TRUE)))
})

test_that("npregiv stop selection is boundary safe and deterministic", {
  select <- npRmpi:::.npregiv_select_stop_index

  expect_identical(select(1)$index, 1L)
  expect_false(select(1)$monotone.failure)
  expect_identical(select(c(1, 2, 3))$index, 1L)
  expect_true(select(c(1, 2, 3))$monotone.failure)
  expect_identical(select(c(3, 2, 1))$index, 3L)
  expect_identical(select(c(1, 3, 2, 1, 2))$index, 4L)
  expect_identical(select(c(1, 3, 2, 2, 3))$index, 3L)
})

test_that("npregiv stopping rows describe their same-index states", {
  set.seed(20260721)
  n <- 28L
  w <- rnorm(n)
  v <- rnorm(n, sd = 0.22)
  z <- 0.45 * w + v
  y <- z^2 - 0.4 * v + rnorm(n, sd = 0.055)

  fit <- suppressWarnings(suppressMessages(npregiv(
    y = y,
    z = z,
    w = w,
    regtype = "lc",
    nmulti = 1L,
    iterate.max = 4L,
    stop.on.increase = FALSE,
    iterate.diff.tol = 0
  )))

  E.y.w <- fitted(npreg(
    tydat = y,
    txdat = w,
    exdat = w,
    bws = fit$bw.E.y.w,
    regtype = "lc"
  ))
  denominator <- sum(E.y.w^2)

  reconstructed <- vapply(seq_along(fit$norm.stop), function(N) {
    state.residual <- fitted(npreg(
      tydat = y - fit$phi.mat[, N],
      txdat = w,
      exdat = w,
      bws = fit$bw.resid.w[N, ],
      regtype = "lc"
    ))
    N * sum(state.residual^2) / denominator
  }, numeric(1))

  expect_equal(fit$norm.stop, reconstructed, tolerance = 1e-12)
  expect_identical(as.numeric(fit$phi),
                   as.numeric(fit$phi.mat[, fit$norm.index]))
})
