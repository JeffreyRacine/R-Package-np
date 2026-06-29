test_that("np.objective.cache controls npindex continuous NN R optimizer caching", {
  old <- options(np.messages = FALSE, np.objective.cache = TRUE)
  on.exit(options(old), add = TRUE)

  run_bw <- function(method, bwtype, cache) {
    set.seed(123 + match(bwtype, c("generalized_nn", "adaptive_nn")) +
               if (identical(method, "kleinspady")) 10L else 0L)
    n <- 70L
    xdat <- data.frame(x1 = runif(n), x2 = runif(n))
    eta <- (xdat$x1 + xdat$x2) / 2
    ydat <- if (identical(method, "kleinspady")) {
      as.integer(runif(n) < plogis(2 * eta - 1))
    } else {
      eta + rnorm(n, sd = 0.25)
    }
    options(np.objective.cache = cache)
    np::npindexbw(
      xdat = xdat,
      ydat = ydat,
      method = method,
      regtype = "lc",
      bwtype = bwtype,
      nmulti = 1L,
      optim.maxit = 35L,
      optim.maxattempts = 1L
    )
  }

  for (method in c("ichimura", "kleinspady")) {
    for (bwtype in c("generalized_nn", "adaptive_nn")) {
      cached <- run_bw(method, bwtype, TRUE)
      uncached <- run_bw(method, bwtype, FALSE)

      expect_equal(cached$beta, uncached$beta)
      expect_equal(cached$bw, uncached$bw)
      expect_equal(cached$fval, uncached$fval, tolerance = 0)
      expect_equal(cached$num.feval, uncached$num.feval)

      expect_equal(unname(cached$nn.cache[["enabled"]]), 1)
      expect_gte(unname(cached$nn.cache[["hits"]]), 0)
      expect_gte(unname(cached$nn.cache[["visits"]]), unname(cached$nn.cache[["raw.evals"]]))
      expect_equal(unname(cached$nn.cache[["visits"]]),
                   unname(cached$nn.cache[["raw.evals"]]) + unname(cached$nn.cache[["hits"]]))
      expect_gte(as.numeric(cached$num.feval.fast[1L]),
                 unname(cached$nn.cache[["hits"]]))
      expect_lte(as.numeric(cached$num.feval.fast[1L]),
                 as.numeric(cached$num.feval[1L]))

      expect_equal(unname(uncached$nn.cache[["enabled"]]), 0)
      expect_equal(unname(uncached$nn.cache[["hits"]]), 0)
    }
  }
})

test_that("npindex R NN cache leaves fixed bandwidth searches unmarked", {
  old <- options(np.messages = FALSE, np.objective.cache = TRUE)
  on.exit(options(old), add = TRUE)

  set.seed(20260524)
  n <- 70L
  xdat <- data.frame(x1 = runif(n), x2 = runif(n))
  ydat <- xdat$x1 + xdat$x2 + rnorm(n, sd = 0.25)

  bw <- np::npindexbw(
    xdat = xdat,
    ydat = ydat,
    method = "ichimura",
    regtype = "lc",
    bwtype = "fixed",
    nmulti = 1L,
    optim.maxit = 25L,
    optim.maxattempts = 1L
  )

  expect_null(bw$nn.cache)
})
