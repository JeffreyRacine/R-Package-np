test_that("fixed full-support tree candidates preserve canonical LP objectives", {
  old <- options(np.messages = FALSE, np.largeh = TRUE,
                 np.macMseries.accelerate = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(20260812L)
  n <- 89L
  all_x <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1),
    x3 = runif(n, -1, 1)
  )
  y <- 0.5 + 0.7 * all_x$x1 - 0.4 * all_x$x2 +
    0.3 * all_x$x1 * all_x$x2 + rnorm(n, sd = 0.25)
  core <- getFromNamespace(".npregbw_call_fixed_degree_core", "np")
  cases <- list(
    list(p = 1L, degree = 0L, bernstein = FALSE,
         kernel = "epanechnikov", order = 2L),
    list(p = 2L, degree = c(1L, 1L), bernstein = TRUE,
         kernel = "uniform", order = 2L),
    list(p = 2L, degree = c(2L, 2L), bernstein = FALSE,
         kernel = "epanechnikov", order = 4L),
    list(p = 3L, degree = rep(2L, 3L), bernstein = TRUE,
         kernel = "epanechnikov", order = 6L),
    list(p = 3L, degree = rep(3L, 3L), bernstein = TRUE,
         kernel = "epanechnikov", order = 8L)
  )

  for (case in cases) {
    x <- all_x[seq_len(case$p)]
    range_x <- vapply(x, function(value) diff(range(value)), numeric(1L))
    bandwidth <- 10 * range_x
    bw_args <- list(
      xdat = x, ydat = y, regtype = "lp", degree = case$degree,
      bernstein.basis = case$bernstein, bwmethod = "cv.ls",
      bwtype = "fixed", ckertype = case$kernel,
      bandwidth.compute = FALSE, bws = bandwidth
    )
    if (case$kernel != "uniform")
      bw_args$ckerorder <- case$order
    bw <- do.call(npregbw, bw_args)

    options(np.tree = FALSE)
    dense <- core(xdat = x, ydat = y, bws = bw, eval.only = TRUE)
    options(np.tree = TRUE)
    tree <- core(xdat = x, ydat = y, bws = bw, eval.only = TRUE)

    expect_equal(
      tree$objective,
      dense$objective,
      tolerance = 2e-11,
      info = paste(case$p, paste(case$degree, collapse = ","),
                   case$bernstein, case$kernel, case$order)
    )
    expect_identical(tree$num.feval.fast, dense$num.feval.fast)
  }

  mixed_x <- data.frame(
    x1 = all_x$x1,
    unordered = factor(rep(c("a", "b", "c"), length.out = n)),
    ordered = ordered(rep(1:4, length.out = n))
  )
  mixed_bw <- npregbw(
    xdat = mixed_x,
    ydat = y,
    regtype = "lp",
    degree = 1L,
    bernstein.basis = TRUE,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    ckertype = "epanechnikov",
    bandwidth.compute = FALSE,
    bws = c(10 * diff(range(mixed_x$x1)), 0.3, 0.3)
  )
  options(np.tree = FALSE)
  mixed_dense <- core(xdat = mixed_x, ydat = y, bws = mixed_bw,
                      eval.only = TRUE)
  options(np.tree = TRUE)
  mixed_tree <- core(xdat = mixed_x, ydat = y, bws = mixed_bw,
                     eval.only = TRUE)
  expect_equal(mixed_tree$objective, mixed_dense$objective,
               tolerance = 2e-11)
  expect_identical(mixed_tree$num.feval.fast, mixed_dense$num.feval.fast)
})

test_that("fixed LP capability reaches the canonical row owner", {
  old <- options(np.messages = FALSE, np.largeh = FALSE,
                 np.macMseries.accelerate = FALSE)
  on.exit({
    options(old)
    Sys.unsetenv("NP_LP_TREE_ORACLE")
  }, add = TRUE)

  set.seed(20260816L)
  n <- 83L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- x$x1 + x$x2 + rnorm(n, sd = 0.25)
  ranges <- vapply(x, function(value) diff(range(value)), numeric(1L))
  core <- getFromNamespace(".npregbw_call_fixed_degree_core", "np")

  evaluate <- function(degree, bandwidth) {
    bw <- npregbw(
      xdat = x, ydat = y, regtype = "lp", degree = degree,
      bwmethod = "cv.ls", bwtype = "fixed",
      ckertype = "epanechnikov", bandwidth.compute = FALSE,
      bws = bandwidth
    )
    options(np.tree = TRUE)
    Sys.setenv(NP_LP_TREE_ORACLE = "1")
    messages <- capture.output(
      tree <- core(xdat = x, ydat = y, bws = bw, eval.only = TRUE),
      type = "message"
    )
    Sys.unsetenv("NP_LP_TREE_ORACLE")
    options(np.tree = FALSE)
    dense <- core(xdat = x, ydat = y, bws = bw, eval.only = TRUE)
    list(
      sparse = any(grepl("NP_LP_TREE_ORACLE", messages, fixed = TRUE)),
      tree = tree,
      dense = dense
    )
  }

  full <- evaluate(c(1L, 1L), 1.25 * ranges / sqrt(5))
  expect_false(full$sparse)
  expect_equal(full$tree$objective, full$dense$objective,
               tolerance = 2e-11)

  prunable <- evaluate(c(1L, 1L), 0.30 * ranges / sqrt(5))
  expect_true(prunable$sparse)
  expect_equal(prunable$tree$objective, prunable$dense$objective,
               tolerance = 2e-11)

  width_two <- evaluate(c(1L, 0L), 1.25 * ranges / sqrt(5))
  expect_true(width_two$sparse)
  expect_equal(width_two$tree$objective, width_two$dense$objective,
               tolerance = 2e-11)
})

test_that("full-support owner is objective-method neutral", {
  old <- options(np.messages = FALSE, np.largeh = TRUE,
                 np.macMseries.accelerate = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(20260813L)
  n <- 97L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  y <- 0.4 + 0.7 * x$x1 - 0.5 * x$x2 + rnorm(n, sd = 0.3)
  bandwidth <- 10 * vapply(x, function(value) diff(range(value)), numeric(1L))
  core <- getFromNamespace(".npregbw_call_fixed_degree_core", "np")

  evaluate <- function(method, objective = "ls") {
    bw <- npregbw(
      xdat = x, ydat = y, regtype = "lp", degree = c(1L, 1L),
      bernstein.basis = TRUE, bwmethod = method, bwtype = "fixed",
      ckertype = "epanechnikov", bandwidth.compute = FALSE,
      bws = bandwidth
    )
    options(np.tree = FALSE)
    dense <- core(xdat = x, ydat = y, bws = bw, eval.only = TRUE,
                  objective = objective)
    options(np.tree = TRUE)
    tree <- core(xdat = x, ydat = y, bws = bw, eval.only = TRUE,
                 objective = objective)
    expect_equal(tree$objective, dense$objective, tolerance = 2e-11,
                 info = paste(method, objective))
    expect_identical(tree$num.feval.fast, dense$num.feval.fast)
  }

  evaluate("cv.ls", "ls")
  evaluate("cv.ls", "ks")
  evaluate("cv.aic", "ls")
})

test_that("full-support owner preserves least-squares quantile CV", {
  old <- options(np.messages = FALSE, np.largeh = TRUE,
                 np.macMseries.accelerate = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(20260814L)
  n <- 89L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  y <- 0.4 + 0.7 * x$x1 - 0.5 * x$x2 + rnorm(n, sd = 0.3)
  h <- 10 * vapply(x, function(value) diff(range(value)), numeric(1L))
  args <- list(
    xdat = x, ydat = y, scale = rep.int(1, n), tau = 0.5,
    bws = h, regtype = "lp", degree = c(2L, 2L),
    bernstein.basis = TRUE, ckertype = "epanechnikov",
    bandwidth.compute = FALSE
  )

  options(np.tree = FALSE)
  dense <- do.call(nplsqregbw, args)
  options(np.tree = TRUE)
  tree <- do.call(nplsqregbw, args)

  expect_equal(tree$objective, dense$objective, tolerance = 2e-11)
  expect_identical(tree$num.feval.fast, dense$num.feval.fast)
})
