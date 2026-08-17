adaptive_gaussian2 <- function(u) stats::dnorm(u)

adaptive_epanechnikov2 <- function(u) {
  edge <- sqrt(5)
  value <- numeric(length(u))
  inside <- abs(u) < edge
  value[inside] <- 3 / (4 * edge) * (1 - u[inside]^2 / 5)
  value
}

adaptive_literal_radius <- function(x, held_out, donor, k) {
  sort(abs(x[-c(held_out, donor)] - x[donor]), method = "quick")[k]
}

adaptive_literal_lc_cvls <- function(x, y, k, kernel,
                                     unordered = NULL,
                                     lambda = numeric()) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (!is.null(unordered))
    unordered <- as.data.frame(unordered)
  fitted <- numeric(n)

  for (held_out in seq_len(n)) {
    numerator <- denominator <- 0
    for (donor in setdiff(seq_len(n), held_out)) {
      weight <- 1
      for (coordinate in seq_len(ncol(x))) {
        bandwidth <- adaptive_literal_radius(
          x[, coordinate], held_out, donor, k[coordinate])
        weight <- weight * kernel(
          (x[held_out, coordinate] - x[donor, coordinate]) /
            bandwidth) / bandwidth
      }
      if (!is.null(unordered)) {
        for (coordinate in seq_len(ncol(unordered))) {
          categories <- length(unique(unordered[[coordinate]]))
          categorical_weight <- if (
            unordered[[coordinate]][held_out] ==
              unordered[[coordinate]][donor]) {
            1 - lambda[coordinate]
          } else {
            lambda[coordinate] / (categories - 1)
          }
          weight <- weight * categorical_weight
        }
      }
      numerator <- numerator + weight * y[donor]
      denominator <- denominator + weight
    }
    fitted[held_out] <- numerator / denominator
  }
  mean((y - fitted)^2)
}

adaptive_package_lc_cvls <- function(x, y, k, ckertype,
                                     unordered = NULL,
                                     lambda = numeric(),
                                     tree = FALSE,
                                     accelerate = FALSE,
                                     invalid.penalty = "baseline") {
  options(np.tree = tree, np.macMseries.accelerate = accelerate,
          np.messages = FALSE)
  x <- as.data.frame(x)
  names(x) <- paste0("x", seq_len(ncol(x)))
  if (!is.null(unordered)) {
    unordered <- as.data.frame(unordered)
    names(unordered) <- paste0("u", seq_len(ncol(unordered)))
    x <- cbind(x, unordered)
  }
  bw <- npregbw(
    xdat = x, ydat = y, bws = c(as.double(k), as.double(lambda)),
    regtype = "lc", bwmethod = "cv.ls", bwtype = "adaptive_nn",
    bwscaling = FALSE, ckertype = ckertype, ckerorder = 2L,
    ukertype = "aitchisonaitken", bandwidth.compute = FALSE)
  as.numeric(np:::.npregbw_eval_only(
    x, y, bw, invalid.penalty = invalid.penalty)$objective)
}

test_that("adaptive regression CVLS uses exact delete-one donor radii", {
  old_options <- options()
  on.exit(options(old_options), add = TRUE)

  fixtures <- list(
    p1 = list(
      x = matrix(c(-2.8, -1.5, -0.7, -0.1, 0.4, 1.2, 2.7, 4.1),
                 ncol = 1L),
      k = 3L),
    p2 = list(
      x = cbind(seq(-1.7, 2.3, length.out = 13L),
                sin(seq_len(13L) * 0.73) + seq_len(13L) / 17),
      k = c(4L, 5L))
  )

  for (fixture in fixtures) {
    x <- fixture$x
    y <- 0.3 + 0.7 * x[, 1L] +
      (if (ncol(x) > 1L) -0.4 * x[, 2L] else 0) +
      sin(seq_len(nrow(x)) / 2.7) / 6
    for (kernel_name in c("gaussian", "epanechnikov")) {
      kernel <- if (kernel_name == "gaussian") {
        adaptive_gaussian2
      } else {
        adaptive_epanechnikov2
      }
      expected <- adaptive_literal_lc_cvls(
        x, y, fixture$k, kernel)
      observed <- adaptive_package_lc_cvls(
        x, y, fixture$k, kernel_name)
      expect_equal(observed, expected,
                   tolerance = 4096 * .Machine$double.eps)
    }
  }
})

test_that("adaptive exact CVLS preserves categorical and owner-axis semantics", {
  old_options <- options()
  on.exit(options(old_options), add = TRUE)

  x <- matrix(c(-2.4, -1.5, -0.8, -0.1, 0.5, 1.1, 1.9, 3.0),
              ncol = 1L)
  unordered <- data.frame(
    group = factor(c("a", "b", "a", "c", "b", "a", "c", "b")))
  y <- 0.4 + 0.6 * x[, 1L] +
    c(a = -0.2, b = 0.1, c = 0.35)[as.character(unordered$group)] +
    sin(seq_len(nrow(x))) / 20
  expected <- adaptive_literal_lc_cvls(
    x, y, 3L, adaptive_gaussian2,
    unordered = unordered, lambda = 0.23)
  observed <- numeric()

  for (tree in c(FALSE, TRUE))
    for (accelerate in c(FALSE, TRUE))
      observed <- c(observed, adaptive_package_lc_cvls(
        x, y, 3L, "gaussian", unordered, 0.23,
        tree = tree, accelerate = accelerate))

  expect_equal(observed, rep(expected, length(observed)),
               tolerance = 4096 * .Machine$double.eps)
})

test_that("adaptive exact CVLS rejects impossible and zero-radius folds", {
  old_options <- options()
  on.exit(options(old_options), add = TRUE)

  zero_radius <- adaptive_package_lc_cvls(
    matrix(c(0, 0, 0, 0, 1, 2), ncol = 1L),
    seq_len(6L) / 10, 2L, "gaussian",
    invalid.penalty = "dbmax")
  impossible_rank <- adaptive_package_lc_cvls(
    matrix(seq_len(6L), ncol = 1L), seq_len(6L),
    5L, "gaussian", invalid.penalty = "dbmax")

  expect_true(is.finite(zero_radius))
  expect_gt(zero_radius, 1e300)
  expect_identical(impossible_rank, .Machine$double.xmax)
})
