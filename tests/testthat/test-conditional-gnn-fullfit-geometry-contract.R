conditional_gnn_fullfit_radius <- function(train, evaluation, k, excluded = integer()) {
  keep <- if (length(excluded)) setdiff(seq_along(train), excluded) else seq_along(train)
  sort(abs(train[keep] - evaluation), method = "radix")[[k]]
}

conditional_gnn_fullfit_oracle <- function(x, y, kx, ky, degree = 0L,
                                            distribution = FALSE,
                                            evaluation = NULL) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  ex <- if (is.null(evaluation)) x else as.matrix(evaluation$x)
  ey <- if (is.null(evaluation)) y else as.numeric(evaluation$y)
  eval_to_train <- if (is.null(evaluation)) seq_len(nrow(x)) else
    rep.int(NA_integer_, nrow(ex))

  vapply(seq_len(nrow(ex)), function(i) {
    excluded <- if (is.na(eval_to_train[[i]])) integer() else eval_to_train[[i]]
    hx <- vapply(seq_len(ncol(x)), function(l) {
      conditional_gnn_fullfit_radius(x[, l], ex[i, l], kx[[l]], excluded)
    }, numeric(1L))
    hy <- conditional_gnn_fullfit_radius(y, ey[[i]], ky, excluded)
    wx <- rep.int(1, nrow(x))
    for (l in seq_len(ncol(x)))
      wx <- wx * dnorm((ex[i, l] - x[, l]) / hx[[l]]) / hx[[l]]
    yrow <- if (distribution) {
      pnorm((ey[[i]] - y) / hy)
    } else {
      dnorm((ey[[i]] - y) / hy) / hy
    }

    if (degree == 0L) {
      influence <- wx / sum(wx)
    } else {
      design <- cbind(1, x[, 1L] - ex[i, 1L])
      influence <- drop(c(1, 0) %*%
        solve(crossprod(design, wx * design), t(design) * rep(wx, each = 2L)))
    }
    sum(influence * yrow)
  }, numeric(1L))
}

conditional_gnn_fullfit_bw <- function(x, y, kx, ky, degree = 0L,
                                        distribution = FALSE) {
  arguments <- list(
    xdat = as.data.frame(x), ydat = data.frame(y = y),
    bws = c(ky, kx), bandwidth.compute = FALSE,
    bwmethod = if (distribution) "cv.ls" else "cv.ml",
    bwtype = "generalized_nn", bwscaling = FALSE,
    cxkertype = "gaussian", cykertype = "gaussian",
    regtype = if (degree == 0L) "lc" else "lp"
  )
  if (distribution)
    arguments$ngrid <- 7L
  if (degree > 0L) {
    arguments$basis <- "glp"
    arguments$degree <- rep.int(degree, ncol(as.matrix(x)))
  }
  do.call(if (distribution) npcdistbw else npcdensbw, arguments)
}

conditional_gnn_fullfit_native <- function(x, y, bw, distribution = FALSE,
                                            evaluation = NULL) {
  arguments <- list(
    bws = bw, txdat = as.data.frame(x), tydat = data.frame(y = y)
  )
  if (!is.null(evaluation)) {
    arguments$exdat <- as.data.frame(evaluation$x)
    arguments$eydat <- data.frame(y = evaluation$y)
  }
  fit <- do.call(if (distribution) npcdist else npcdens, arguments)
  if (distribution) fit$condist else fit$condens
}

test_that("conditional generalized-NN fits delete mapped radius occurrences", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(20260818)
  n <- 13L
  x1 <- sort(runif(n, -0.9, 1.1))
  x2 <- runif(n, -0.6, 0.8)
  y <- 0.8 * x1 - 0.35 * x2 + rnorm(n, sd = 0.17)
  cases <- list(
    list(x = data.frame(x1 = x1), kx = 4L, ky = 5L, degree = 0L),
    list(x = data.frame(x1 = x1, x2 = x2), kx = c(4L, 6L), ky = 5L, degree = 0L),
    list(x = data.frame(x1 = x1), kx = 6L, ky = 5L, degree = 1L)
  )

  for (case in cases) {
    for (distribution in c(FALSE, TRUE)) {
      expected <- conditional_gnn_fullfit_oracle(
        case$x, y, case$kx, case$ky, case$degree, distribution
      )
      bw <- conditional_gnn_fullfit_bw(
        case$x, y, case$kx, case$ky, case$degree, distribution
      )
      for (tree in c(FALSE, TRUE)) {
        options(np.tree = tree)
        observed <- conditional_gnn_fullfit_native(
          case$x, y, bw, distribution
        )
        expect_equal(observed, expected, tolerance = 2e-10)
      }
    }
  }
})

test_that("equal-valued external queries retain external generalized-NN geometry", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(20260818)
  x <- data.frame(x = sort(runif(13L, -0.9, 1.1)))
  y <- 0.8 * x$x + rnorm(13L, sd = 0.17)
  evaluation <- list(x = x[c(2L, 7L, 12L), , drop = FALSE],
                     y = y[c(2L, 7L, 12L)])
  bw <- conditional_gnn_fullfit_bw(x, y, 4L, 5L)
  expected <- conditional_gnn_fullfit_oracle(
    x, y, 4L, 5L, evaluation = evaluation
  )

  observed <- conditional_gnn_fullfit_native(x, y, bw, evaluation = evaluation)
  expect_equal(observed, expected, tolerance = 2e-12)
})

test_that("conditional generalized-NN fit uncertainty surfaces retain tree parity", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(20260820)
  x <- data.frame(x = sort(runif(13L, -0.9, 1.1)))
  y <- 0.8 * x$x + rnorm(13L, sd = 0.17)

  for (distribution in c(FALSE, TRUE)) {
    bw <- conditional_gnn_fullfit_bw(x, y, 6L, 5L, 1L, distribution)
    fits <- lapply(c(FALSE, TRUE), function(tree) {
      options(np.tree = tree)
      do.call(if (distribution) npcdist else npcdens, list(
        bws = bw, txdat = x, tydat = data.frame(y = y), gradients = TRUE
      ))
    })
    fields <- c(if (distribution) "condist" else "condens",
                "conderr", "congrad", "congerr")
    for (field in fields) {
      expect_true(all(is.finite(fits[[1L]][[field]])))
      expect_equal(fits[[2L]][[field]], fits[[1L]][[field]], tolerance = 2e-10)
    }
  }
})

test_that("conditional generalized-NN fits reject zero literal radii", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = c(0, 0, 0, 1, 2))
  y <- c(0, 0, 0, 0.5, 1)

  for (degree in 0:1) {
    bw <- conditional_gnn_fullfit_bw(x, y, 2L, 2L, degree)
    expect_error(
      conditional_gnn_fullfit_native(x, y, bw),
      "zero literal radius",
      fixed = TRUE
    )
  }
})
