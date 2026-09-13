# Ordered, uncached scalar-tile oracle: retain the pre-reuse arithmetic.
j1_uncached <- function(bws, xdat, index, response.matrix, pivotal,
                        structural = NULL, context = "bootstrap statistic") {
  tile <- getFromNamespace(".np_npsig_streamed_iid_tile", "npRmpi")
  value <- numeric(ncol(response.matrix))
  for (j in seq_along(index))
    value <- value + tile(bws, xdat, index[[j]],
      response.matrix = response.matrix, null.mean = response.matrix[, 1L],
      residual.pool = response.matrix[, 1L], pivotal = pivotal,
      structural = if (is.null(structural)) NULL else structural[, j, drop = FALSE],
      context = context) / length(index)
  value
}

test_that("joint reuse preserves pure reducer ordering, zeros and diagnostics", {
  ns <- asNamespace("npRmpi")
  joint <- get(".np_npsig_streamed_response_statistic", ns)
  tile <- get(".np_npsig_streamed_iid_tile", ns)
  sandbox <- new.env(parent = ns)
  environment(joint) <- environment(tile) <- sandbox
  sandbox$.np_npsig_streamed_iid_tile <- tile
  state <- new.env(parent = emptyenv())
  state$calls <- integer()
  state$fits <- list()
  sandbox$.npRmpi_npsig_npreg_local <- function(bws, txdat, tydat, ...) {
    key <- as.integer(tydat[[1L]])
    state$calls <- c(state$calls, key)
    state$fits[[key]]
  }
  # No native execution in this pure fixture.
  sandbox$.npreghat_exact_lp_apply_from_regression_core <- function(...) {
    stop("unexpected native tile")
  }
  x <- data.frame(u = factor(c("a", "b", "b")),
                  o = ordered(c(1, 2, 3)), x = c(-1, 0, 1))
  bw <- list(icon = c(FALSE, FALSE, TRUE), iuno = c(TRUE, FALSE, FALSE),
              iord = c(FALSE, TRUE, FALSE), ncon = 1L)
  responses <- matrix(rep(1:2, each = 3L), 3L)
  masks <- matrix(FALSE, 3L, 2L)
  fit <- list(grad = cbind(u = c(0, 0, 0), o = c(1, 2, 3), x = c(1, 1, 1)),
               gerr = cbind(u = rep(NA_real_, 3L), o = c(2, 2, 2), x = c(1, 1, 1)))
  state$fits <- list(fit, fit)
  set.seed(719)
  seed <- .Random.seed
  expect_identical(joint(bw, x, 1:2, responses, TRUE, masks),
                   rep(mean((c(1, 2, 3) / 2)^2) / 2, 2L))
  expect_identical(state$calls, 1:2)
  expect_identical(.Random.seed, seed)
  expect_identical(state$fits, list(fit, fit)) # reducer must not mutate cached fits
  zero <- fit
  zero$grad[,] <- 0
  zero$gerr <- NULL
  state$fits <- list(zero, zero)
  expect_identical(joint(bw, x, 1:2, responses, TRUE, masks), c(0, 0))

  # A later component's response 1 must not preempt component 1, response 2.
  first <- second <- fit
  first$gerr[1L, 2L] <- 0
  second$grad[1L, 1L] <- 1
  second$gerr[1L, 1L] <- 0
  state$fits <- list(first, second)
  state$calls <- integer()
  expect_error(joint(bw, x, 1:2, responses, TRUE, masks, context = "witness"),
    "witness \\(tile column 2\\).*zero standard error.*'u'.*row 1")
  expect_identical(state$calls, 1:2)
  second$grad[1L, 1L] <- Inf
  state$fits[[2L]] <- second
  expect_error(joint(bw, x, 1:2, responses, TRUE, masks),
               "non-finite gradient estimates")
  second$grad[1L, 1L] <- 1e-200
  state$fits[[2L]] <- second
  expect_error(joint(bw, x, 1:2, responses, TRUE, masks),
               "zero standard error.*'u'")

  # A structural component does not fit or require its nonfinite gradient.
  state$fits <- list(fit, fit)
  masks[, 1L] <- TRUE
  state$calls <- integer()
  expect_identical(joint(bw, x, 1:2, responses, TRUE, masks),
                   rep(mean((c(1, 2, 3) / 2)^2) / 2, 2L))
  expect_identical(state$calls, 1:2)
  expect_error(joint(bw, x, 1:2, responses[, FALSE, drop = FALSE], TRUE, masks),
               "subscript out of bounds|finite n-by-at-most-8")
  bad <- responses
  bad[2L, 2L] <- NA_real_
  expect_error(joint(bw, x, 1:2, bad, TRUE, masks),
               "finite n-by-at-most-8")
})

test_that("actual joint tiles fit once per response and preserve exact component order", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(120927)
  n <- 36L
  x <- data.frame(u = factor(rep(letters[1:3], 12L)),
    o = ordered(rep(1:4, 9L)), x = runif(n, -1, 1))
  y <- x$x + .4 * (x$u == "b") + rnorm(n, sd = .3)
  responses <- vapply(seq_len(8L), function(k) y + rnorm(n, sd = .2), numeric(n))
  joint <- getFromNamespace(".np_npsig_streamed_response_statistic", "npRmpi")
  tile <- getFromNamespace(".np_npsig_streamed_iid_tile", "npRmpi")
  owner <- getFromNamespace(".npRmpi_npsig_npreg_local", "npRmpi")
  count <- new.env(parent = emptyenv())
  count$n <- 0L
  wrapped <- function(...) {
    count$n <- count$n + 1L
    owner(...)
  }
  local_mocked_bindings(.npRmpi_npsig_npreg_local = wrapped, .package = "npRmpi")
  for (regtype in c("lc", "ll", "lp")) {
    bw <- npregbw(xdat = x, ydat = y, bws = c(.25, .3, .65),
      bandwidth.compute = FALSE, regtype = regtype,
      degree = if (regtype == "lp") 2L else NULL)
    for (width in c(1L, 7L, 8L, 1L)) {
      response <- responses[, seq_len(width), drop = FALSE]
      for (index in list(1:3, c(2L, 3L, 1L), c(3L, 1L, 2L), c(2L, 1L, 2L))) {
        count$n <- 0L
        expected <- j1_uncached(bw, x, index, response, TRUE)
        expect_identical(count$n, as.integer(sum(index != 3L) * width))
        count$n <- 0L
        actual <- joint(bw, x, index, response, TRUE)
        expect_identical(actual, expected)
        expect_identical(count$n, width)
      }
    }
    for (index in list(1L, c(3L, 1L), 3L)) {
      count$n <- 0L
      expected <- j1_uncached(bw, x, index, responses, TRUE)
      before <- count$n
      count$n <- 0L
      expect_identical(joint(bw, x, index, responses, TRUE), expected)
      expect_identical(count$n, before)
    }
    count$n <- 0L
    expect_identical(joint(bw, x, c(2L, 3L, 1L), responses, FALSE),
                     j1_uncached(bw, x, c(2L, 3L, 1L), responses, FALSE))
    expect_identical(count$n, 0L)
    donor <- matrix(rep(seq_len(n), 2L), nrow = n)
    ready <- y + matrix(y[donor], nrow = n)
    expect_identical(
      tile(bw, x, 1L, donor.index = donor, null.mean = y,
           residual.pool = y, pivotal = TRUE),
      tile(bw, x, 1L, response.matrix = ready, null.mean = y,
           residual.pool = y, pivotal = TRUE))
  }
})


test_that("public MPI joint smoke retains reproducibility, metadata and RNG", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(727)
  n <- 36L
  x <- data.frame(u = factor(rep(1:3, 12L)),
                   o = ordered(rep(1:4, 9L)), x = runif(n, -1, 1))
  y <- sin(x$x) + as.integer(x$u) / 3 + rnorm(n, sd = .2)
  bw <- npregbw(xdat = x, ydat = y, bws = c(.25, .3, .65),
                 bandwidth.compute = FALSE, regtype = "ll")
  expect_error(npsigtest(bw, xdat = x, ydat = y, index = c(1L, 1L), B = 9L),
               "repeated|unique")
  for (method in c("iid", "wild-rademacher")) {
    for (pivot in c(FALSE, TRUE)) {
      for (index in list(1:3, c(2L, 3L, 1L))) {
        set.seed(171)
        before <- .Random.seed
        # This is a repeatability smoke, not an uncached fanout oracle.
        # Installed pre/post worker proof is a separate qualification gate:
        # a master-only namespace mock does not prove a remote owner changed.
        expected <- npsigtest(bw, xdat = x, ydat = y, index = index, joint = TRUE,
          B = 9L, pivot = pivot, boot.method = method, random.seed = 727)
        after <- .Random.seed
        expect_identical(after, before)
        actual <- npsigtest(bw, xdat = x, ydat = y, index = index, joint = TRUE,
          B = 9L, pivot = pivot, boot.method = method, random.seed = 727)
        expect_identical(actual, expected)
        expect_identical(.Random.seed, after)
      }
    }
  }
})
