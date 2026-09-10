beta_cvml_guard <- function(fit) {
  vapply(fit, function(value) {
    if (value > .Machine$double.xmin)
      return(-log(value))
    if (value < -.Machine$double.xmin)
      return(log(-value) - 2 * log(.Machine$double.xmin))
    -log(.Machine$double.xmin)
  }, numeric(1L))
}

beta_cvml_weights <- function(data, bandwidth, type, kernel, order,
                              density = FALSE, evaluation = NULL) {
  arguments <- list(
    bws = bandwidth,
    txdat = data,
    bwtype = type,
    ckertype = kernel,
    ckerorder = order
  )
  if (!is.null(evaluation))
    arguments$exdat <- evaluation
  if (identical(kernel, "beta")) {
    arguments$ckerbound <- "fixed"
    arguments$ckerlb <- 0
    arguments$ckerub <- 1
  }
  if (!density) {
    arguments$return.kernel.weights <- TRUE
    return(do.call(npksum, arguments)$kw)
  }
  # Identity response columns expose the exact response-density operator.
  # For peer kernels public returned `kw` deliberately omits this division.
  arguments$tydat <- diag(nrow(data))
  arguments$bandwidth.divide <- TRUE
  out <- do.call(npksum, arguments)$ksum
  if (is.null(dim(out))) as.numeric(out) else t(out)
}

beta_cvml_adaptive_full_sample_radius <- function(train, k) {
  vapply(seq_along(train), function(index) {
    distance <- abs(train - train[[index]])
    duplicates <- sum(distance == 0) - 1L
    positive <- sort(distance[distance > 0])
    if (duplicates >= k) positive[[1L]] else
      positive[[max(1L, k - duplicates)]]
  }, numeric(1L))
}

beta_cvml_ratio_oracle <- function(x, y, bandwidth, type,
                                    xkernel, ykernel, order) {
  xweights <- beta_cvml_weights(
    x, bandwidth[[2L]], type, xkernel, order, density = FALSE
  )
  if (identical(type, "adaptive_nn") &&
      !identical(xkernel, "beta")) {
    xbandwidth <- beta_cvml_adaptive_full_sample_radius(
      x[[1L]], bandwidth[[2L]]
    )
    xweights <- xweights /
      matrix(xbandwidth, nrow = nrow(x), ncol = nrow(x))
  }
  yweights <- beta_cvml_weights(
    y, bandwidth[[1L]], type, ykernel, order, density = TRUE
  )
  fit <- vapply(seq_len(nrow(x)), function(evaluation) {
    xrow <- xweights[, evaluation]
    yrow <- yweights[, evaluation]
    xrow[[evaluation]] <- 0
    yrow[[evaluation]] <- 0
    sum(xrow * yrow) / sum(xrow)
  }, numeric(1L))

  # The public eval-only wrapper reports the maximized sign of the native
  # negative-log CVML loss.
  -sum(beta_cvml_guard(fit))
}

beta_cvml_native_objective <- function(x, y, bandwidth, type,
                                       xkernel, ykernel, order, compress) {
  old <- options(np.categorical.compress = compress)
  on.exit(options(old), add = TRUE)
  arguments <- list(
    xdat = x,
    ydat = y,
    bws = bandwidth,
    bandwidth.compute = FALSE,
    bwtype = type,
    bwmethod = "cv.ml",
    regtype = "lc",
    cxkertype = xkernel,
    cxkerorder = order,
    cykertype = ykernel,
    cykerorder = order
  )
  if (identical(xkernel, "beta")) {
    arguments$cxkerbound <- "fixed"
    arguments$cxkerlb <- 0
    arguments$cxkerub <- 1
  }
  if (identical(ykernel, "beta")) {
    arguments$cykerbound <- "fixed"
    arguments$cykerlb <- 0
    arguments$cykerub <- 1
  }
  bandwidth_object <- do.call(npcdensbw, arguments)
  evaluator <- getFromNamespace(".npcdensbw_eval_only", "npRmpi")
  as.numeric(evaluator(x, y, bandwidth_object)$objective[[1L]])
}

.beta_cvml_test_env <- new.env(parent = emptyenv())

.ensure_beta_cvml_pool <- function() {
  if (!isTRUE(.beta_cvml_test_env$started)) {
    skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
    .beta_cvml_test_env$started <- TRUE
    withr::defer({
      if (isTRUE(.beta_cvml_test_env$started)) {
        try(close_mpi_slaves(force = TRUE), silent = TRUE)
        .beta_cvml_test_env$started <- FALSE
      }
    }, envir = testthat::teardown_env())
  }
}

beta_cvml_mpi_domains <- function() {
  list(
    type = c("fixed", "generalized_nn", "adaptive_nn"),
    order = c(2L, 4L, 6L, 8L),
    sides = c("beta_beta", "beta_legacy", "legacy_beta"),
    compress = c(FALSE, TRUE)
  )
}

beta_cvml_mpi_pairwise_cases <- function() {
  domains <- beta_cvml_mpi_domains()
  types <- domains$type
  orders <- domains$order
  sides <- domains$sides

  do.call(rbind, lapply(seq_along(orders), function(order_index) {
    data.frame(
      type = types,
      order = orders[[order_index]],
      sides = sides[
        (seq_along(types) + order_index - 2L) %% length(sides) + 1L
      ],
      compress = (seq_along(types) + order_index) %% 2L == 0L,
      stringsAsFactors = FALSE
    )
  }))
}

test_that("MPI conditional beta CVML covers every route-axis pair", {
  domains <- beta_cvml_mpi_domains()
  cases <- beta_cvml_mpi_pairwise_cases()
  expect_identical(anyDuplicated(cases), 0L)
  for (pair in combn(names(cases), 2L, simplify = FALSE)) {
    observed <- unique(cases[pair])
    complete <- expand.grid(
      domains[pair], stringsAsFactors = FALSE
    )
    expect_setequal(
      do.call(paste, observed),
      do.call(paste, complete)
    )
  }
})

test_that("conditional beta CVML equals the public signed weight-ratio oracle", {
  skip_on_cran()
  .ensure_beta_cvml_pool()
  set.seed(20260801L)
  n <- 31L
  xbeta <- data.frame(x = sort(runif(n, 0.03, 0.97)))
  ybeta <- data.frame(y = sort(runif(n, 0.04, 0.96)))
  xlegacy <- data.frame(x = qnorm(xbeta$x))
  ylegacy <- data.frame(y = qnorm(ybeta$y))

  cases <- beta_cvml_mpi_pairwise_cases()
  for (case_index in seq_len(nrow(cases))) {
    type <- cases$type[[case_index]]
    order <- cases$order[[case_index]]
    sides <- cases$sides[[case_index]]
    compress <- cases$compress[[case_index]]
    bandwidth <- if (identical(type, "fixed")) c(0.19, 0.17) else c(8, 7)
    xkernel <- if (identical(sides, "legacy_beta")) "gaussian" else "beta"
    ykernel <- if (identical(sides, "beta_legacy")) "gaussian" else "beta"
    x <- if (identical(xkernel, "beta")) xbeta else xlegacy
    y <- if (identical(ykernel, "beta")) ybeta else ylegacy
    expected <- beta_cvml_ratio_oracle(
      x, y, bandwidth, type, xkernel, ykernel, order
    )
    actual <- beta_cvml_native_objective(
      x, y, bandwidth, type, xkernel, ykernel, order, compress
    )
    expect_equal(
      actual,
      expected,
      tolerance = 8e-12,
      info = paste(type, order, sides, "compress", compress)
    )
  }
})
