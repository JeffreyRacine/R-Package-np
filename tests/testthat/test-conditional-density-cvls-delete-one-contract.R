conditional_cvls_nn_generalized_radius <- function(train, evaluation, k) {
  vapply(evaluation, function(target) {
    distance <- abs(train - target)
    exact <- sum(distance == 0)
    positive <- sort(distance[distance > 0])
    if (exact >= k) positive[[1L]] else positive[[k - exact]]
  }, numeric(1L))
}

conditional_cvls_nn_adaptive_radius <- function(train, k) {
  vapply(seq_along(train), function(index) {
    distance <- abs(train - train[[index]])
    duplicate <- sum(distance == 0) - 1L
    positive <- sort(distance[distance > 0])
    if (duplicate >= k) positive[[1L]] else
      positive[[max(1L, k - duplicate)]]
  }, numeric(1L))
}

conditional_cvls_bandwidth_matrices <- function(data, bws, bwtype,
                                                convolution = FALSE) {
  n <- nrow(data)
  ordinary <- matrix(1, nrow = n, ncol = n)
  companion <- matrix(1, nrow = n, ncol = n)

  for (coordinate in seq_along(data)) {
    train <- data[[coordinate]]
    if (identical(bwtype, "fixed")) {
      eval_bandwidth <- rep.int(bws[[coordinate]], n)
      train_bandwidth <- eval_bandwidth
    } else if (identical(bwtype, "generalized_nn")) {
      eval_bandwidth <- conditional_cvls_nn_generalized_radius(
        train, train, bws[[coordinate]]
      )
      train_bandwidth <- eval_bandwidth
    } else {
      eval_bandwidth <- conditional_cvls_nn_generalized_radius(
        train, train, bws[[coordinate]]
      )
      train_bandwidth <- conditional_cvls_nn_adaptive_radius(
        train, bws[[coordinate]]
      )
    }

    if (identical(bwtype, "adaptive_nn")) {
      ordinary <- ordinary * matrix(train_bandwidth, n, n)
    } else {
      ordinary <- ordinary * matrix(eval_bandwidth, n, n, byrow = TRUE)
    }
    if (convolution) {
      companion <- companion * outer(train_bandwidth, eval_bandwidth)
    }
  }

  if (convolution) companion else ordinary
}

conditional_cvls_gaussian_response_matrices <- function(data, bws, bwtype) {
  stopifnot(ncol(data) == 1L)
  ordinary_bandwidth <- conditional_cvls_bandwidth_matrices(
    data, bws, bwtype, convolution = FALSE
  )
  if (identical(bwtype, "adaptive_nn")) {
    train_bandwidth <- conditional_cvls_nn_adaptive_radius(
      data[[1L]], bws[[1L]]
    )
    eval_bandwidth <- train_bandwidth
  } else if (identical(bwtype, "generalized_nn")) {
    eval_bandwidth <- conditional_cvls_nn_generalized_radius(
      data[[1L]], data[[1L]], bws[[1L]]
    )
    train_bandwidth <- eval_bandwidth
  } else {
    train_bandwidth <- rep.int(bws[[1L]], nrow(data))
    eval_bandwidth <- train_bandwidth
  }
  difference <- outer(data[[1L]], data[[1L]], function(train, evaluation) {
    evaluation - train
  })
  ordinary <- dnorm(difference / ordinary_bandwidth) / ordinary_bandwidth
  convolution_scale <- sqrt(
    outer(train_bandwidth^2, eval_bandwidth^2, `+`)
  )
  convolution <- dnorm(difference / convolution_scale) / convolution_scale
  list(normal = ordinary, convolution = convolution)
}

conditional_cvls_delete_one_oracle <- function(xdat, ydat, bw) {
  n <- nrow(xdat)
  basis <- npRmpi:::W.lp(
    xdat = xdat,
    degree = bw$degree.engine,
    basis = bw$basis.engine,
    bernstein.basis = bw$bernstein.basis.engine
  )
  xraw <- npksum(
    bws = bw$xbw,
    txdat = xdat,
    exdat = xdat,
    bwtype = bw$type,
    ckertype = bw$cxkertype,
    ckerorder = bw$cxkerorder,
    operator = "normal",
    return.kernel.weights = TRUE,
    bandwidth.divide = FALSE
  )$kw
  if (identical(bw$type, "adaptive_nn")) {
    xraw <- xraw / conditional_cvls_bandwidth_matrices(
      xdat, bw$xbw, bw$type, convolution = FALSE
    )
  }
  stopifnot(identical(bw$cykertype, "gaussian"), bw$cykerorder == 2L)
  ykernel <- conditional_cvls_gaussian_response_matrices(
    ydat, bw$ybw, bw$type
  )
  ynormal <- ykernel$normal
  yconvolution <- ykernel$convolution
  delete_one <- matrix(0, nrow = n, ncol = n)

  for (evaluation in seq_len(n)) {
    weight <- xraw[, evaluation]
    gram <- crossprod(basis, basis * weight)
    coefficient <- solve(gram, basis[evaluation, ])
    full <- as.numeric(weight * (basis %*% coefficient))
    denominator <- 1 - full[[evaluation]]
    stopifnot(is.finite(denominator), denominator != 0)
    full[[evaluation]] <- 0
    delete_one[, evaluation] <- full / denominator
  }

  linear <- colSums(delete_one * ynormal)
  quadratic <- vapply(seq_len(n), function(evaluation) {
    smoother <- delete_one[, evaluation]
    sum((smoother %o% smoother) * yconvolution)
  }, numeric(1L))
  -mean(quadratic - 2 * linear)
}

.conditional_delete_one_test_env <- new.env(parent = emptyenv())

.ensure_conditional_delete_one_pool <- function() {
  if (!isTRUE(.conditional_delete_one_test_env$started)) {
    npRmpi.init(nslaves = 1L, quiet = TRUE)
    .conditional_delete_one_test_env$started <- TRUE
    withr::defer({
      if (isTRUE(.conditional_delete_one_test_env$started)) {
        try(npRmpi.quit(force = TRUE), silent = TRUE)
        .conditional_delete_one_test_env$started <- FALSE
      }
    }, envir = testthat::teardown_env())
  }
}

test_that("fixed Gaussian conditional CVLS uses one delete-one smoother", {
  skip_on_cran()
  .ensure_conditional_delete_one_pool()
  x <- c(-1.22, -0.83, -0.31, 0.04, 0.37, 0.79, 1.18)
  y <- c(-1.41, -0.88, -0.46, -0.03, 0.34, 0.82, 1.29)
  hx <- 0.43
  hy <- 0.37
  wx <- outer(x, x, function(train, evaluation) {
    dnorm((evaluation - train) / hx) / hx
  })
  wy <- outer(y, y, function(train, evaluation) {
    dnorm((evaluation - train) / hy) / hy
  })
  overlap <- outer(y, y, function(first, second) {
    dnorm((first - second) / (sqrt(2) * hy)) / (sqrt(2) * hy)
  })
  oracle <- 0

  for (evaluation in seq_along(x)) {
    keep <- setdiff(seq_along(x), evaluation)
    xweight <- wx[keep, evaluation]
    xweight <- xweight / sum(xweight)
    integrated.square <- sum(
      (xweight %o% xweight) * overlap[keep, keep, drop = FALSE]
    )
    cross.fit <- sum(xweight * wy[keep, evaluation])
    oracle <- oracle + integrated.square - 2 * cross.fit
  }
  oracle <- oracle / length(x)

  bw <- npcdensbw(
    xdat = data.frame(x = x), ydat = data.frame(y = y),
    bws = c(hy, hx), bandwidth.compute = FALSE,
    bwmethod = "cv.ls", bwtype = "fixed", bwscaling = FALSE,
    regtype = "lc", cxkertype = "gaussian", cxkerorder = 2,
    cykertype = "gaussian", cykerorder = 2
  )
  native <- npRmpi:::.npcdensbw_eval_only(
    data.frame(x = x), data.frame(y = y), bw
  )$objective

  expect_equal(native, -oracle, tolerance = 5e-10)
})

test_that("conditional LP CVLS matches independent delete-one WLS", {
  skip_on_cran()
  .ensure_conditional_delete_one_pool()
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = FALSE, np.largeh = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(2026080101L)
  n <- 54L
  xdat <- data.frame(
    x1 = runif(n, 0.05, 0.95),
    x2 = runif(n, 0.04, 0.96)
  )
  ydat <- data.frame(
    y = 0.4 + 0.25 * sin(2 * pi * xdat$x1) - 0.12 * xdat$x2 +
      rnorm(n, sd = 0.07)
  )
  cases <- list(
    list(type = "fixed", kernel = "epanechnikov", order = 6L,
         degree = c(1L, 1L), bernstein = FALSE,
         bws = c(0.31, 0.44, 0.47)),
    list(type = "fixed", kernel = "gaussian", order = 4L,
         degree = c(2L, 2L), bernstein = TRUE,
         bws = c(0.26, 0.29, 0.32)),
    list(type = "generalized_nn", kernel = "gaussian", order = 8L,
         degree = c(2L, 2L), bernstein = FALSE,
         bws = c(21L, 24L, 25L)),
    list(type = "adaptive_nn", kernel = "gaussian", order = 2L,
         degree = c(1L, 1L), bernstein = TRUE,
         bws = c(22L, 24L, 25L))
  )

  for (case in cases) {
    bw <- npcdensbw(
      xdat = xdat,
      ydat = ydat,
      bws = case$bws,
      bandwidth.compute = FALSE,
      bwmethod = "cv.ls",
      bwtype = case$type,
      regtype = "lp",
      degree = case$degree,
      bernstein.basis = case$bernstein,
      cxkertype = case$kernel,
      cxkerorder = case$order,
      cykertype = "gaussian",
      cykerorder = 2L
    )
    native <- npRmpi:::.npcdensbw_eval_only(xdat, ydat, bw)$objective
    oracle <- conditional_cvls_delete_one_oracle(xdat, ydat, bw)

    expect_equal(
      as.numeric(native), oracle, tolerance = 2e-8,
      info = paste(case$type, case$kernel, case$order,
                   if (case$bernstein) "bernstein" else "raw")
    )
  }
})

test_that("all-large CVLS is the same delete-one estimator", {
  skip_on_cran()
  .ensure_conditional_delete_one_pool()
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = FALSE, np.largeh = TRUE)
  on.exit(options(old), add = TRUE)

  set.seed(2026080102L)
  n <- 54L
  xdat <- data.frame(x1 = runif(n), x2 = runif(n))
  ydat <- data.frame(
    y = 0.3 + 0.2 * xdat$x1 - 0.15 * xdat$x2 + rnorm(n, sd = 0.12)
  )

  for (bernstein in c(FALSE, TRUE)) {
    bw <- npcdensbw(
      xdat = xdat,
      ydat = ydat,
      bws = rep.int(1e8, 3L),
      bandwidth.compute = FALSE,
      bwmethod = "cv.ls",
      bwtype = "fixed",
      regtype = "lp",
      degree = c(1L, 1L),
      bernstein.basis = bernstein,
      cxkertype = "gaussian",
      cxkerorder = 2L,
      cykertype = "gaussian",
      cykerorder = 2L
    )
    options(np.largeh = TRUE)
    accelerated <- npRmpi:::.npcdensbw_eval_only(xdat, ydat, bw)
    options(np.largeh = FALSE)
    ordinary <- npRmpi:::.npcdensbw_eval_only(xdat, ydat, bw)
    oracle <- conditional_cvls_delete_one_oracle(xdat, ydat, bw)

    expect_lt(abs(accelerated$objective - oracle), 1e-17)
    expect_lt(abs(ordinary$objective - oracle), 1e-17)
    expect_equal(as.numeric(accelerated$num.feval.fast[[1L]]), 1)
    expect_equal(as.numeric(ordinary$num.feval.fast[[1L]]), 0)
  }
})

test_that("scalar fixed-row shortcut is normal-operator only", {
  roots <- unique(c(
    Sys.getenv("NP_SOURCE_ROOT", unset = ""),
    normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE),
    normalizePath(getwd(), mustWork = FALSE)
  ))
  path <- NULL
  for (root in roots[nzchar(roots)]) {
    candidate <- file.path(root, "src", "jksum.c")
    if (file.exists(candidate)) {
      path <- candidate
      break
    }
  }
  skip_if(is.null(path), "package source unavailable")
  engine <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_match(
    engine,
    paste0(
      "(ctx->operator_y != NULL) &&\n",
      "     (ctx->operator_y[0] == OP_NORMAL))\n",
      "    return np_conditional_y_scalar_fixed_row_direct("
    ),
    fixed = TRUE
  )
  expect_match(
    engine,
    "if(ctx->operator_y == NULL || ctx->operator_y[0] != OP_NORMAL)",
    fixed = TRUE
  )
  expect_false(grepl("NP_BOUNDED_CVLS_I1_MODE_FULL", engine, fixed = TRUE))
  expect_false(grepl("NP_BOUNDED_CVLS_I1_MODE_BOOK", engine, fixed = TRUE))
})
