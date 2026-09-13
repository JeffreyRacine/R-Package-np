test_that("categorical response statistics use each response's paired-contrast HC0", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(120926)
  n <- 48L
  x <- data.frame(u = factor(rep(letters[1:3], 16L)),
    o = ordered(rep(c("low", "mid", "high", "top"), 12L),
                levels = c("low", "mid", "high", "top")),
    x = runif(n, -1, 1))
  y <- x$x + .4 * (x$u == "b") + rnorm(n, sd = .3)
  response <- cbind(y, y + rnorm(n, sd = .2), rev(y))
  tile <- getFromNamespace(".np_npsig_streamed_iid_tile", "np")
  joint <- getFromNamespace(".np_npsig_streamed_response_statistic", "np")
  for (regtype in c("lc", "ll", "lp")) {
    bw <- npregbw(xdat = x, ydat = y, bws = c(.25, .3, .65),
      bandwidth.compute = FALSE, regtype = regtype,
      degree = if (regtype == "lp") 2L else NULL)
    direct <- matrix(0, ncol(response), ncol(x))
    for (k in seq_len(ncol(response))) {
      fit <- npreg(bws = bw, txdat = x, tydat = response[, k],
                   gradients = TRUE, se = TRUE)
      for (i in seq_len(ncol(x))) {
        # Independent row convention: unordered reference only; ordered
        # first-level rows have a forward difference and remain required.
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
      expect_identical(result, direct[, i])
    }
    expect_equal(joint(bw, x, 1:3, response, pivotal = TRUE),
                 rowMeans(direct), tolerance = 2e-10)
  }
})

test_that("an exactly zero observed contrast skips bootstrap but preserves later draws", {
  old <- options(np.messages = FALSE, np.largelambda = TRUE)
  on.exit(options(old), add = TRUE)
  set.seed(110926)
  n <- 48L
  x <- data.frame(u = factor(rep(c("a", "b"), 24L)), x = runif(n, -1, 1))
  y <- x$x + rnorm(n, sd = .3)
  bw <- npregbw(xdat = x, ydat = y, bws = c(.5, .6),
                regtype = "ll", bandwidth.compute = FALSE)
  seed.enter <- getFromNamespace(".np_seed_enter", "np")
  for (method in c("iid", "wild", "wild-rademacher", "pairwise")) {
    set.seed(91)
    outer <- .Random.seed
    both <- npsigtest(bw, B = 9L, boot.method = method, random.seed = 81)
    expect_identical(.Random.seed, outer)
    expect_identical(both$bootstrap.executed, c(0L, 9L))
    expect_identical(both$In[[1L]], 0)
    expect_identical(both$P[[1L]], 1)
    expect_true(all(is.na(both$In.bootstrap[, 1L])))
    # Reproduce the later stream independently, using the established donor
    # laws. random.seed=NULL is not a continuation request: set.seed(NULL)
    # initializes a new state, so it must not be used as this oracle.
    later <- testthat::with_mocked_bindings(
      npsigtest(bw, index = 2L, B = 9L, boot.method = method, random.seed = 81),
      .np_seed_enter = function(random.seed) {
        state <- seed.enter(random.seed)
        for (k in seq_len(9L)) {
          if (method %in% c("iid", "pairwise")) sample.int(n, replace = TRUE)
          else runif(n)
        }
        state
      }, .package = "np")
    expect_identical(both$In[[2L]], later$In[[1L]])
    expect_identical(both$In.bootstrap[, 2L], later$In.bootstrap[, 1L])
  }
  all.zero <- npsigtest(bw, index = 1L, B = 9L, joint = TRUE)
  expect_identical(all.zero$bootstrap.executed, 0L)
  expect_identical(all.zero$P, 1)
  expect_true(all(is.na(all.zero$In.bootstrap)))
})

test_that("whole-zero effects need no SE and retain the joint denominator", {
  statistic <- getFromNamespace(".np_npsig_statistic", "np")
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
  statistic <- getFromNamespace(".np_npsig_statistic", "np")
  zero.effects <- getFromNamespace(".np_npsig_zero_effects", "np")
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

test_that("Type II zero effects preserve later draws and reselection hot starts", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  n <- 36L
  x <- data.frame(u = factor(rep(1:3, length.out = n)),
                  x = seq(-1.5, 1.5, length.out = n))
  y <- x$x + as.integer(x$u) / 2 + sin(seq_len(n)) / 5
  bw <- npregbw(xdat = x, ydat = y, bws = c(1, .6),
    bandwidth.compute = FALSE, regtype = "lc", ukertype = "liracine")
  seed.enter <- getFromNamespace(".np_seed_enter", "np")
  run <- function(index, advance = FALSE) {
    calls <- results <- list()
    set.seed(91)
    outer <- .Random.seed
    value <- testthat::with_mocked_bindings(
      npsigtest(bw, xdat = x, ydat = y, index = index, B = 9L,
                 boot.type = "II", boot.method = "iid", random.seed = 81),
      # Exercise the real reselection helper and its nmulti policy, but do
      # not run a bandwidth search. Retain every supplied response and seed.
      npregbw = function(...) {
        args <- list(...)
        k <- length(calls) + 1L
        calls[[k]] <<- args
        selected <- args$bws
        selected$bw[2L] <- .55 + .01 * k
        results[[k]] <<- selected
        selected
      },
      .np_seed_enter = function(random.seed) {
        state <- seed.enter(random.seed)
        if (advance) for (k in seq_len(9L)) sample.int(n, replace = TRUE)
        state
      }, .package = "np")
    expect_identical(.Random.seed, outer)
    list(value = value, calls = calls, selected = results)
  }
  both <- run(1:2)
  later <- run(2L, advance = TRUE)
  expect_identical(both$value$bootstrap.executed, c(0L, 9L))
  expect_identical(both$value$P[[1L]], 1)
  expect_true(all(is.na(both$value$In.bootstrap[, 1L])))
  expect_length(both$calls, 9L)
  expect_identical(both$calls, later$calls)
  expect_identical(both$calls[[1L]]$bws, bw)
  expect_false("nmulti" %in% names(both$calls[[1L]]))
  for (k in 2:9) {
    expect_identical(both$calls[[k]]$bws, both$selected[[k - 1L]])
    expect_identical(both$calls[[k]]$nmulti, 1L)
  }
  expect_identical(both$value$In[[2L]], later$value$In[[1L]])
  expect_identical(both$value$In.bootstrap[, 2L], later$value$In.bootstrap[, 1L])
  expect_identical(both$value$P[[2L]], later$value$P[[1L]])
  expected.bw <- bw
  expected.bw$bw[2L] <- both$selected[[9L]]$bw[2L]
  expect_identical(both$value$bws, expected.bw)
})
