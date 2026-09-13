test_that("whole-zero effects need no SE and retain the joint denominator", {
  statistic <- getFromNamespace(".np_npsig_statistic", "npRmpi")
  zero <- list(grad = matrix(0, 3L, 2L), gerr = NULL)
  expect_identical(statistic(zero, 1:2, TRUE), 0)
  expect_identical(statistic(zero, 1:2, FALSE), mean(zero$grad^2))
  mixed <- list(grad = cbind(u = c(0, 0, 0), x = c(1, 2, 3)),
                gerr = cbind(u = rep(NA_real_, 3L), x = rep(2, 3L)))
  expect_identical(statistic(mixed, 1:2, TRUE),
                   mean(c(0, 0, 0, (c(1, 2, 3) / 2)^2)))
  expect_identical(statistic(mixed, 1L, TRUE), 0)
  expect_identical(statistic(mixed, 1:2, FALSE), mean(mixed$grad^2))
  mixed$grad[1L, 2L] <- 0
  mixed$gerr[1L, 2L] <- 0
  expect_error(statistic(mixed, 1:2, TRUE, context = "bootstrap replication 3"),
    "bootstrap replication 3.*zero standard error.*'x'.*row 1")
})

test_that("zero effects are decided before squaring and require finite gradients", {
  statistic <- getFromNamespace(".np_npsig_statistic", "npRmpi")
  zero.effects <- getFromNamespace(".np_npsig_zero_effects", "npRmpi")
  tiny <- list(grad = matrix(rep(1e-200, 2L), ncol = 1L),
               gerr = matrix(1, 2L, 1L))
  expect_identical(statistic(tiny, 1L, TRUE), 0)
  expect_identical(statistic(tiny, 1L, FALSE), 0)
  expect_identical(zero.effects(tiny, 1L), FALSE)
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    nonfinite <- list(grad = matrix(c(0, bad), ncol = 1L), gerr = NULL)
    expect_identical(zero.effects(nonfinite, 1L), FALSE)
    expect_error(statistic(nonfinite, 1L, TRUE), "non-finite gradient estimates")
    expect_error(statistic(nonfinite, 1L, FALSE), "non-finite gradient estimates")
  }
})

test_that("zero-effect draw advancement matches the MPI seed planner", {
  planner <- getFromNamespace(".npRmpi_npsig_bootstrap_seed_plan", "npRmpi")
  advance <- getFromNamespace(".np_npsig_advance_bootstrap_rng", "npRmpi")
  wild <- function(n, a, b, p.a) ifelse(runif(n) <= p.a, a, b)
  for (method in c("iid", "pairwise", "wild", "wild-rademacher")) {
    planned <- withr::with_seed(81, {
      planner(36L, 9L, method, wild, -.6180339887499, 1.6180339887499,
              if (method == "wild-rademacher") .5 else .72360679774998)
      .Random.seed
    })
    skipped <- withr::with_seed(81, {
      advance(36L, 9L, method)
      .Random.seed
    })
    expect_identical(skipped, planned, info = method)
  }
})

test_that("categorical MPI tiles use each response's paired-contrast HC0", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old), add = TRUE)
  set.seed(120926)
  n <- 48L
  x <- data.frame(u = factor(rep(letters[1:3], 16L)),
    o = ordered(rep(c("low", "mid", "high", "top"), 12L),
                levels = c("low", "mid", "high", "top")),
    x = runif(n, -1, 1))
  y <- x$x + .4 * (x$u == "b") + rnorm(n, sd = .3)
  response <- cbind(y, y + rnorm(n, sd = .2), rev(y))
  tile <- getFromNamespace(".np_npsig_streamed_iid_tile", "npRmpi")
  joint <- getFromNamespace(".np_npsig_streamed_response_statistic", "npRmpi")
  for (regtype in c("lc", "ll", "lp")) {
    bw <- npregbw(xdat = x, ydat = y, bws = c(.25, .3, .65),
      bandwidth.compute = FALSE, regtype = regtype,
      degree = if (regtype == "lp") 2L else NULL)
    direct <- matrix(0, ncol(response), ncol(x))
    for (k in seq_len(ncol(response))) {
      fit <- npreg(bws = bw, txdat = x, tydat = response[, k],
                   gradients = TRUE, se = TRUE)
      for (i in seq_len(ncol(x))) {
        # Only unordered reference rows have identical endpoints here.
        # Ordered first-level rows use a forward difference and are required.
        reference <- if (i == 1L) as.integer(x$u) == 1L else rep(FALSE, n)
        expect_true(all(is.finite(fit$gerr[!reference, i]) & fit$gerr[!reference, i] > 0))
        ratio <- numeric(n)
        ratio[!reference] <- fit$grad[!reference, i] / fit$gerr[!reference, i]
        direct[k, i] <- mean(ratio^2)
      }
    }
    for (i in 1:2) {
      result <- tile(bw, x, i, response.matrix = response,
                     null.mean = y, residual.pool = y, pivotal = TRUE)
      expect_equal(result, direct[, i], tolerance = 2e-10)
    }
    expect_equal(joint(bw, x, 1:3, response, pivotal = TRUE),
                 rowMeans(direct), tolerance = 2e-10)
  }
})

test_that("observed-zero MPI tests retain shape, raw arithmetic and later draws", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old), add = TRUE)
  n <- 36L
  x <- data.frame(u = factor(rep(1:3, length.out = n)),
    w = factor(rep(3:1, length.out = n)), x = seq(-1.5, 1.5, length.out = n))
  y <- x$x + as.integer(x$u) / 2 + sin(seq_len(n)) / 5
  bw <- npregbw(xdat = x, ydat = y, bws = c(1, 1, .6),
    bandwidth.compute = FALSE, regtype = "lc", ukertype = "liracine")
  fit <- npreg(bws = bw, txdat = x, tydat = y, gradients = TRUE, se = TRUE)
  expect_true(all(is.finite(fit$grad[, 1:2]) & fit$grad[, 1:2] == 0))
  for (pivot in c(FALSE, TRUE)) {
    set.seed(91)
    outer <- .Random.seed
    both <- npsigtest(bw, xdat = x, ydat = y, B = 9L,
      pivot = pivot, boot.method = "iid", random.seed = 81)
    expect_identical(.Random.seed, outer)
    expect_identical(both$bootstrap.executed, c(0L, 0L, 9L))
    expect_identical(both$P[1:2], c(1, 1))
    expect_true(all(is.na(both$In.bootstrap[, 1:2])))
    expect_true(all(is.finite(both$In.bootstrap[, 3L])))
    expect_identical(dim(both$In.bootstrap), c(9L, 3L))
    expect_output(print(both), "Analytic non-rejection", fixed = TRUE)

    # Recreate the later predictor's IID responses without namespace mocks:
    # MPI workers receive planned per-replication seeds, not a serial mock.
    scaled <- scale(residuals(npreg(bws = bw, txdat = x, tydat = y)))
    null.frame <- x
    null.frame$x <- uocquantile(x$x, .5)
    null.mean <- npreg(bws = bw, txdat = x, tydat = y, exdat = null.frame)$mean
    residual.pool <- as.numeric(scale(y - null.mean) * attr(scaled, "scaled:scale") +
                                attr(scaled, "scaled:center"))
    residual.pool <- residual.pool - mean(residual.pool)
    responses <- withr::with_seed(81, {
      for (k in seq_len(18L)) sample.int(n, replace = TRUE)
      vapply(seq_len(9L), function(k)
        null.mean + residual.pool[sample.int(n, replace = TRUE)], numeric(n))
    })
    oracle <- vapply(seq_len(9L), function(k) {
      boot.fit <- npreg(bws = bw, txdat = x, tydat = responses[, k],
                        gradients = TRUE, se = pivot)
      if (pivot) mean((boot.fit$grad[, 3L] / boot.fit$gerr[, 3L])^2)
      else mean(boot.fit$grad[, 3L]^2)
    }, numeric(1L))
    expect_equal(both$In.bootstrap[, 3L], oracle, tolerance = 2e-10)

    # Type II must return before bandwidth reselection when all tested
    # observed contrasts are zero. No mocked serial reselection assumptions.
    all.zero <- npsigtest(bw, xdat = x, ydat = y, index = 1:2,
      B = 9L, joint = TRUE, pivot = pivot, boot.type = "II")
    expect_identical(all.zero$In, 0)
    expect_identical(all.zero$P, 1)
    expect_identical(all.zero$bootstrap.executed, 0L)
    expect_identical(all.zero$bootstrap.reason, "all observed contrasts identically zero")
    expect_identical(dim(all.zero$In.bootstrap), c(9L, 1L))
    expect_true(all(is.na(all.zero$In.bootstrap)))

    mixed.joint <- npsigtest(bw, xdat = x, ydat = y, B = 9L,
      joint = TRUE, pivot = pivot, random.seed = 81)
    expected <- if (pivot) {
      ratio <- matrix(0, n, 3L)
      ratio[, 3L] <- fit$grad[, 3L] / fit$gerr[, 3L]
      mean(ratio^2)
    } else mean(fit$grad^2)
    expect_equal(mixed.joint$In, expected, tolerance = 2e-10)
    expect_identical(mixed.joint$bootstrap.executed, 9L)
    expect_true(all(is.finite(mixed.joint$In.bootstrap)))
  }
})
