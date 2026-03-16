npreghat_generalized_lp_owner_case <- function() {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- "NPREGHAT_GENERALIZED_LP_OWNER_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "run_eval_case <- function(train_df, y, eval_df, degree, basis = 'glp', bern = FALSE, tol = 1e-10) {",
    "  frame <- train_df",
    "  frame$y <- y",
    "  fml <- as.formula(paste('y ~', paste(names(train_df), collapse = ' + ')))",
    "  bw <- npregbw(fml, data = frame, regtype = 'lp', degree = degree, basis = basis, bernstein.basis = bern, bwtype = 'generalized_nn')",
    "  H <- npreghat(bws = bw, txdat = train_df, exdat = eval_df)",
    "  a <- npreghat(bws = bw, txdat = train_df, exdat = eval_df, y = y, output = 'apply')",
    "  g <- npreg(bws = bw, exdat = eval_df)",
    "  stopifnot(isTRUE(all.equal(drop(H %*% y), a, tolerance = tol)))",
    "  stopifnot(isTRUE(all.equal(drop(H %*% y), g$mean, tolerance = tol)))",
    "}",
    "set.seed(101)",
    "n <- 80L",
    "x <- runif(n, -1, 1)",
    "y <- x^2 + rnorm(n, sd = 0.15)",
    "xe <- data.frame(x = seq(-1, 1, length.out = 41L))",
    "run_eval_case(data.frame(x = x), y, xe, degree = 2L)",
    "run_eval_case(data.frame(x = x), y, xe, degree = 3L)",
    "set.seed(102)",
    "n <- 70L",
    "x1 <- runif(n, -1, 1)",
    "x2 <- runif(n, -1, 1)",
    "y <- x1^2 - x2 + x1 * x2 + rnorm(n, sd = 0.1)",
    "xe <- expand.grid(x1 = seq(-1, 1, length.out = 6L), x2 = seq(-1, 1, length.out = 6L))",
    "for (b in c('glp', 'additive', 'tensor')) {",
    "  run_eval_case(data.frame(x1 = x1, x2 = x2), y, xe, degree = c(2L, 2L), basis = b)",
    "}",
    "set.seed(103)",
    "n <- 70L",
    "x1 <- c(0, 1, runif(n - 2L))",
    "x2 <- c(1, 0, runif(n - 2L))",
    "y <- sin(pi * x1) + x2^2 + rnorm(n, sd = 0.1)",
    "xe <- expand.grid(x1 = seq(0.1, 0.9, length.out = 6L), x2 = seq(0.1, 0.9, length.out = 6L))",
    "run_eval_case(data.frame(x1 = x1, x2 = x2), y, xe, degree = c(2L, 2L), basis = 'tensor', bern = TRUE, tol = 1e-9)",
    "set.seed(104)",
    "n <- 70L",
    "x1 <- runif(n, -1, 1)",
    "x2 <- runif(n, -1, 1)",
    "u <- factor(sample(c('a', 'b', 'c'), n, replace = TRUE))",
    "o <- ordered(sample(1:3, n, replace = TRUE))",
    "y <- x1^2 + x2 - 0.5 * (u == 'b') + 0.3 * (o == 3) + rnorm(n, sd = 0.1)",
    "xe <- data.frame(x1 = seq(-1, 1, length.out = 12L), x2 = seq(-1, 1, length.out = 12L), u = factor(rep(c('a', 'b', 'c'), length.out = 12L), levels = levels(u)), o = ordered(rep(1:3, length.out = 12L), levels = levels(o)))",
    "run_eval_case(data.frame(x1 = x1, x2 = x2, u = u, o = o), y, xe, degree = c(2L, 2L), basis = 'tensor')",
    "old_tree <- getOption('np.tree')",
    "options(np.tree = TRUE)",
    "on.exit(options(np.tree = old_tree), add = TRUE)",
    "set.seed(203)",
    "n <- 80L",
    "x <- runif(n, -1, 1)",
    "y <- x^2 + rnorm(n, sd = 0.2)",
    "d <- data.frame(x = x)",
    "bw <- npregbw(y ~ x, data = data.frame(d, y = y), regtype = 'lp', degree = 2L, basis = 'glp', bwtype = 'generalized_nn')",
    "H <- npreghat(bws = bw, txdat = d)",
    "a <- npreghat(bws = bw, txdat = d, y = y, output = 'apply')",
    "g <- npreg(bws = bw)",
    "stopifnot(isTRUE(all.equal(drop(H %*% y), a, tolerance = 1e-10)))",
    "stopifnot(isTRUE(all.equal(drop(H %*% y), g$mean, tolerance = 1e-10)))",
    sprintf("cat('%s\\n')", ok_tag)
  )

  res <- npRmpi_run_rscript_subprocess(
    lines = lines,
    timeout = 180L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl(ok_tag, res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
}

test_that("npRmpi generalized higher-degree lp owner matches npreg and apply on common cells", {
  npreghat_generalized_lp_owner_case()
})

npreghat_generalized_lp_derivative_owner_case <- function() {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- "NPREGHAT_GENERALIZED_LP_DERIVATIVE_OWNER_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "run_eval_case <- function(train_df, y, eval_df, degree, s, basis = 'glp', bern = FALSE, tol = 1e-9) {",
    "  frame <- train_df",
    "  frame$y <- y",
    "  fml <- as.formula(paste('y ~', paste(names(train_df), collapse = ' + ')))",
    "  bw <- npregbw(fml, data = frame, regtype = 'lp', degree = degree, basis = basis, bernstein.basis = bern, bwtype = 'generalized_nn')",
    "  H <- npreghat(bws = bw, txdat = train_df, exdat = eval_df, s = s)",
    "  a <- npreghat(bws = bw, txdat = train_df, exdat = eval_df, y = y, output = 'apply', s = s)",
    "  g <- npreg(bws = bw, exdat = eval_df, gradients = TRUE)",
    "  target.col <- which(as.integer(s) == 1L)[1L]",
    "  stopifnot(isTRUE(all.equal(drop(H %*% y), a, tolerance = tol)))",
    "  stopifnot(isTRUE(all.equal(drop(H %*% y), g$grad[, target.col], tolerance = tol)))",
    "}",
    "set.seed(301)",
    "n <- 80L",
    "x <- runif(n, -1, 1)",
    "y <- x^2 + rnorm(n, sd = 0.15)",
    "xe <- data.frame(x = seq(-1, 1, length.out = 41L))",
    "run_eval_case(data.frame(x = x), y, xe, degree = 2L, s = 1L)",
    "run_eval_case(data.frame(x = x), y, xe, degree = 3L, s = 1L)",
    "set.seed(302)",
    "n <- 70L",
    "x1 <- runif(n, -1, 1)",
    "x2 <- runif(n, -1, 1)",
    "y <- x1^2 - x2 + x1 * x2 + rnorm(n, sd = 0.1)",
    "xe <- expand.grid(x1 = seq(-1, 1, length.out = 6L), x2 = seq(-1, 1, length.out = 6L))",
    "run_eval_case(data.frame(x1 = x1, x2 = x2), y, xe, degree = c(2L, 2L), s = c(1L, 0L), basis = 'tensor')",
    sprintf("cat('%s\\n')", ok_tag)
  )

  res <- npRmpi_run_rscript_subprocess(
    lines = lines,
    timeout = 120L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl(ok_tag, res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
}

test_that("npRmpi generalized higher-degree lp derivative owner matches npreg and apply", {
  npreghat_generalized_lp_derivative_owner_case()
})
