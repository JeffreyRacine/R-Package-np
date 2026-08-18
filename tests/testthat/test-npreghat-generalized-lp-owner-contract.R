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

npreghat_lp_gradient_apply_contract_case <- function() {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- "NPREGHAT_LP_GRADIENT_APPLY_CONTRACT_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "set.seed(20260630)",
    "n <- 18L",
    "x1 <- runif(n, -1, 1)",
    "x2 <- runif(n, -1, 1)",
    "f <- factor(sample(letters[1:3], n, replace = TRUE))",
    "tx <- data.frame(x1 = x1, x2 = x2, f = f)",
    "eta <- x1 + 0.5 * x2 + 0.2 * (as.integer(f) - 2L)",
    "y1 <- sin(eta) + 0.15 * eta^2",
    "y2 <- cos(eta) + seq_len(n) / (10 * n)",
    "Y <- cbind(y1 = y1, y2 = y2)",
    "check_apply_contract <- function(degree, s, tol = 1e-9) {",
    "  bw <- npregbw(xdat = tx, ydat = y1, regtype = 'lp', degree = degree, bwtype = 'generalized_nn', bws = c(7, 7, 0.5), bandwidth.compute = FALSE)",
    "  H <- unclass(npreghat(bws = bw, txdat = tx, s = s))",
    "  a1 <- npreghat(bws = bw, txdat = tx, y = matrix(y1, ncol = 1L), output = 'apply', s = s)",
    "  a2 <- npreghat(bws = bw, txdat = tx, y = Y, output = 'apply', s = s)",
    "  stopifnot(isTRUE(all.equal(as.vector(a1), as.vector(H %*% y1), tolerance = tol)))",
    "  stopifnot(isTRUE(all.equal(unname(as.matrix(a2)), unname(H %*% Y), tolerance = tol)))",
    "  if (identical(as.integer(degree), c(1L, 0L)) && identical(as.integer(s), c(1L, 0L))) {",
    "    fit <- npreg(bws = bw, txdat = tx, tydat = y1, gradients = TRUE, warn.glp.gradient = FALSE)",
    "    stopifnot(isTRUE(all.equal(as.vector(H %*% y1), as.vector(fit$grad[, 1L]), tolerance = 1e-8)))",
    "  }",
    "}",
    "check_apply_contract(degree = c(0L, 0L), s = c(1L, 0L))",
    "check_apply_contract(degree = c(1L, 0L), s = c(1L, 0L))",
    sprintf("cat('%s\\n')", ok_tag)
  )

  res <- npRmpi_run_rscript_subprocess(
    lines = lines,
    timeout = 90L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl(ok_tag, res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
}

test_that("npRmpi lp gradient apply routing covers degree-zero and heterogeneous degrees", {
  npreghat_lp_gradient_apply_contract_case()
})

npreghat_gennn_mixed_degree_available_derivative_case <- function() {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- "NPREGHAT_GENNN_MIXED_DEGREE_AVAILABLE_DERIVATIVE_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "set.seed(20260701)",
    "old_tree <- getOption('np.tree')",
    "on.exit(options(np.tree = old_tree), add = TRUE)",
    "n <- 24L",
    "x1 <- runif(n, -1, 1)",
    "x2 <- runif(n, -1, 1)",
    "u <- factor(sample(letters[1:3], n, replace = TRUE))",
    "o <- ordered(sample(1:3, n, replace = TRUE))",
    "y <- sin(x1) + 0.25 * x1^2 + 0.15 * x2 + 0.2 * (as.integer(u) - 2L)",
    "Y <- cbind(y = y, shifted = y + seq_len(n) / (10 * n))",
    "check_case <- function(tx, ex, bws, degree, s_available, s_unavailable) {",
    "  bw <- npregbw(xdat = tx, ydat = y, regtype = 'lp', degree = degree, bwtype = 'generalized_nn', bws = bws, bandwidth.compute = FALSE)",
    "  for (tree in c(FALSE, TRUE)) {",
    "    options(np.tree = tree)",
    "    H <- npreghat(bws = bw, txdat = tx, exdat = ex, output = 'matrix', s = s_available)",
    "    a1 <- npreghat(bws = bw, txdat = tx, exdat = ex, y = matrix(y, ncol = 1L), output = 'apply', s = s_available)",
    "    a2 <- npreghat(bws = bw, txdat = tx, exdat = ex, y = Y, output = 'apply', s = s_available)",
    "    fit <- npreg(bws = bw, txdat = tx, tydat = y, exdat = ex, gradients = TRUE, warn.glp.gradient = FALSE)",
    "    stopifnot(isTRUE(all.equal(as.vector(a1), as.vector(unclass(H) %*% y), tolerance = 1e-9)))",
    "    stopifnot(isTRUE(all.equal(unname(as.matrix(a2)), unname(unclass(H) %*% Y), tolerance = 1e-9)))",
    "    stopifnot(isTRUE(all.equal(as.vector(unclass(H) %*% y), as.vector(fit$grad[, which(s_available > 0L)]), tolerance = 1e-8)))",
    "    err <- tryCatch(npreghat(bws = bw, txdat = tx, exdat = ex, output = 'matrix', s = s_unavailable), error = conditionMessage)",
    "    stopifnot(is.character(err), grepl(\"requested derivative order in 's' exceeds local polynomial degree\", err, fixed = TRUE))",
    "  }",
    "}",
    "tx.cont <- data.frame(x1 = x1, x2 = x2)",
    "ex.cont <- tx.cont[seq(2L, n, by = 3L), , drop = FALSE]",
    "check_case(tx.cont, ex.cont, bws = c(9, 9), degree = c(2L, 0L), s_available = c(1L, 0L), s_unavailable = c(0L, 1L))",
    "check_case(tx.cont, ex.cont, bws = c(9, 9), degree = c(0L, 2L), s_available = c(0L, 1L), s_unavailable = c(1L, 0L))",
    "tx.uno <- data.frame(x1 = x1, x2 = x2, u = u)",
    "ex.uno <- tx.uno[seq(2L, n, by = 3L), , drop = FALSE]",
    "check_case(tx.uno, ex.uno, bws = c(9, 9, 0.5), degree = c(2L, 0L), s_available = c(1L, 0L), s_unavailable = c(0L, 1L))",
    "check_case(tx.uno, ex.uno, bws = c(9, 9, 0.5), degree = c(0L, 2L), s_available = c(0L, 1L), s_unavailable = c(1L, 0L))",
    "tx.ord <- data.frame(x1 = x1, x2 = x2, o = o)",
    "ex.ord <- tx.ord[seq(2L, n, by = 3L), , drop = FALSE]",
    "check_case(tx.ord, ex.ord, bws = c(9, 9, 0.5), degree = c(2L, 0L), s_available = c(1L, 0L), s_unavailable = c(0L, 1L))",
    "check_case(tx.ord, ex.ord, bws = c(9, 9, 0.5), degree = c(0L, 2L), s_available = c(0L, 1L), s_unavailable = c(1L, 0L))",
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

test_that("npRmpi generalized-nn lp mixed-degree available derivative owner is nonrecursive", {
  npreghat_gennn_mixed_degree_available_derivative_case()
})

npreghat_tree_disabled_lp_scalar_apply_guard_case <- function() {
  skip_on_cran()

  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  ok_tag <- "NPREGHAT_TREE_DISABLED_LP_SCALAR_APPLY_GUARD_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "old_tree <- getOption('np.tree')",
    "options(np.tree=FALSE)",
    "on.exit(options(np.tree=old_tree), add=TRUE)",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "set.seed(20260630)",
    "n <- 24L",
    "x1 <- runif(n, -1, 1)",
    "x2 <- runif(n, -1, 1)",
    "f <- factor(sample(letters[1:3], n, replace = TRUE))",
    "tx <- data.frame(x1 = x1, x2 = x2, f = f)",
    "eta <- x1 + 0.5 * x2 + 0.2 * (as.integer(f) - 2L)",
    "y <- sin(eta) + 0.15 * eta^2",
    "check_scalar_apply_contract <- function(degree, s, tol = 1e-9) {",
    "  bw <- npregbw(xdat = tx, ydat = y, regtype = 'lp', degree = degree, bwtype = 'generalized_nn', bws = c(9, 9, 0.5), bandwidth.compute = FALSE)",
    "  H <- unclass(npreghat(bws = bw, txdat = tx, output = 'matrix', s = s))",
    "  a <- npreghat(bws = bw, txdat = tx, y = matrix(y, ncol = 1L), output = 'apply', s = s)",
    "  stopifnot(isTRUE(all.equal(as.vector(a), as.vector(H %*% y), tolerance = tol)))",
    "}",
    "check_scalar_apply_contract(degree = c(2L, 2L), s = c(0L, 0L))",
    "check_scalar_apply_contract(degree = c(2L, 2L), s = c(1L, 0L))",
    sprintf("cat('%s\\n')", ok_tag)
  )

  res <- npRmpi_run_rscript_subprocess(
    lines = lines,
    timeout = 90L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl(ok_tag, res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
}

test_that("npRmpi tree-disabled higher-degree lp scalar apply routes through exact matrix contract", {
  npreghat_tree_disabled_lp_scalar_apply_guard_case()
})

test_that("public positive LP mean matrices use the typed native capability", {
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "installed npRmpi unavailable for subprocess proof")
  ok_tag <- "NPRMPI_PUBLIC_LP_NATIVE_MATRIX_CONTRACT_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "npRmpi.init(nslaves=1L, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "candidate <- getFromNamespace('.npreghat_native_matrix_candidate', 'npRmpi')",
    "has_extended_nn <- getFromNamespace('npRegressionHasExtendedNn', 'npRmpi')",
    "native_matrix <- getFromNamespace('.npreghat_exact_lp_matrix_from_regression_core', 'npRmpi')",
    "local_regression <- getFromNamespace('.npRmpi_with_local_regression', 'npRmpi')",
    "set.seed(20260816)",
    "n <- 42L",
    "x <- data.frame(x=sort(runif(n, -1, 1)))",
    "y <- sin(pi*x$x) + seq_len(n)/1000",
    "make_bw <- function(type, bws, degree=2L, kernel='gaussian') npregbw(xdat=x, ydat=y, bws=bws, bandwidth.compute=FALSE, bwmethod='cv.ls', bwtype=type, bwscaling=FALSE, regtype='lp', degree=degree, degree.select='manual', basis='glp', bernstein.basis=FALSE, ckertype=kernel, ckerorder=2L)",
    "bw.fixed <- make_bw('fixed', 0.32)",
    "bw.gnn <- make_bw('generalized_nn', 9L)",
    "bw.ann <- make_bw('adaptive_nn', 9L)",
    "candidate_args <- function(bw, output='matrix', s=0L) list(bws=bw, output=output, regtype.engine='lp', degree=2L, basis='glp', bernstein.basis=FALSE, s=as.integer(s), leave.one.out=FALSE)",
    "stopifnot(do.call(candidate, candidate_args(bw.fixed)))",
    "stopifnot(do.call(candidate, candidate_args(bw.gnn)))",
    "stopifnot(!do.call(candidate, candidate_args(bw.ann)))",
    "bw.extended <- bw.gnn",
    "bw.extended$bw[bw.extended$icon] <- bw.extended$nobs + 1L",
    "stopifnot(has_extended_nn(bw.extended))",
    "stopifnot(!do.call(candidate, candidate_args(bw.extended)))",
    "stopifnot(do.call(candidate, candidate_args(bw.fixed, output='constraint')))",
    "stopifnot(!do.call(candidate, candidate_args(bw.fixed, s=1L)))",
    "for (bw in list(bw.fixed, bw.gnn)) {",
    "  public <- npreghat(bws=bw, txdat=x, output='matrix')",
    "  native <- local_regression(native_matrix(bw, txdat=x, degree=2L, basis='glp', bernstein.basis=FALSE, s=0L))",
    "  stopifnot(identical(as.vector(public), as.vector(native)))",
    "  stopifnot(identical(attr(public, 'ridge.used', exact=TRUE), attr(native, 'ridge.used', exact=TRUE)))",
    "  stopifnot(inherits(public, 'npreghat'))",
    "  stopifnot(identical(attr(public, 'trainiseval', exact=TRUE), TRUE))",
    "  public.with.y <- npreghat(bws=bw, txdat=x, y=y, output='matrix')",
    "  stopifnot(identical(attr(public.with.y, 'Hy', exact=TRUE), as.vector(public.with.y %*% y)))",
    "  constraint <- npreghat(bws=bw, txdat=x, y=y, output='constraint')",
    "  stopifnot(identical(as.vector(constraint), as.vector(t(public)*y)))",
    "}",
    "bw.low.rank <- make_bw('fixed', 0.04, degree=3L, kernel='epanechnikov')",
    "H.low.rank <- npreghat(bws=bw.low.rank, txdat=x, output='matrix')",
    "ridge.used <- attr(H.low.rank, 'ridge.used', exact=TRUE)",
    "stopifnot(is.double(ridge.used), length(ridge.used)==nrow(H.low.rank))",
    "stopifnot(all(is.finite(ridge.used)), all(ridge.used >= 0))",
    "stopifnot(sum(ridge.used > 0) > 0L)",
    sprintf("cat('%s\\n')", ok_tag)
  )
  res <- npRmpi_run_rscript_subprocess(lines=lines, timeout=90L, env=env)
  info <- paste(res$output, collapse="\n")
  expect_equal(res$status, 0L, info=info)
  expect_true(any(grepl(ok_tag, res$output, fixed=TRUE)), info=info)
})

test_that("public positive LP mean multi-RHS apply uses the canonical response owner", {
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "installed npRmpi unavailable for subprocess proof")
  ok_tag <- "NPRMPI_PUBLIC_LP_MULTI_RHS_APPLY_CONTRACT_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "options(npRmpi.autodispatch=TRUE, np.messages=FALSE)",
    "npRmpi.init(nslaves=1L, quiet=TRUE)",
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "candidate <- getFromNamespace('.npreghat_native_positive_lp_mean_apply_candidate', 'npRmpi')",
    "has_extended_nn <- getFromNamespace('npRegressionHasExtendedNn', 'npRmpi')",
    "native_apply <- getFromNamespace('.npreghat_exact_lp_apply_from_regression_core', 'npRmpi')",
    "local_regression <- getFromNamespace('.npRmpi_with_local_regression', 'npRmpi')",
    "set.seed(20260817)",
    "n <- 42L",
    "x <- data.frame(x1=runif(n,-1,1), x2=runif(n,-1,1), u=factor(rep(c('a','b','c'), length.out=n)))",
    "y1 <- sin(pi*x$x1) + 0.2*x$x2 + 0.1*(x$u=='b')",
    "y2 <- cos(pi*x$x2) - 0.3*x$x1 + seq_len(n)/1000",
    "Y <- cbind(y1,y2)",
    "make_bw <- function(type,bws) npregbw(xdat=x, ydat=y1, bws=bws, bandwidth.compute=FALSE, bwmethod='cv.ls', bwtype=type, bwscaling=FALSE, regtype='lp', degree=c(2L,1L), degree.select='manual', basis='glp', bernstein.basis=FALSE, ckertype='gaussian', ckerorder=2L)",
    "bw.fixed <- make_bw('fixed', c(0.48,0.52,0.35))",
    "bw.gnn <- make_bw('generalized_nn', c(11L,11L,0.35))",
    "bw.ann <- make_bw('adaptive_nn', c(11L,11L,0.35))",
    "candidate_args <- function(bw,y=Y,s=c(0L,0L)) list(bws=bw, output='apply', y=y, regtype.engine='lp', degree=c(2L,1L), basis='glp', bernstein.basis=FALSE, s=as.integer(s), leave.one.out=FALSE)",
    "stopifnot(do.call(candidate, candidate_args(bw.fixed)))",
    "stopifnot(do.call(candidate, candidate_args(bw.gnn)))",
    "stopifnot(!do.call(candidate, candidate_args(bw.ann)))",
    "stopifnot(!do.call(candidate, candidate_args(bw.fixed, y=Y[,1L,drop=FALSE])))",
    "stopifnot(!do.call(candidate, candidate_args(bw.fixed, s=c(1L,0L))))",
    "bw.extended <- bw.gnn",
    "bw.extended$bw[bw.extended$icon] <- bw.extended$nobs + 1L",
    "stopifnot(has_extended_nn(bw.extended))",
    "stopifnot(!do.call(candidate, candidate_args(bw.extended)))",
    "for (bw in list(bw.fixed,bw.gnn)) {",
    "  options(np.npreghat.apply.memory.threshold.mb=Inf)",
    "  public.inf <- npreghat(bws=bw, txdat=x, y=Y, output='apply')",
    "  options(np.npreghat.apply.memory.threshold.mb=0)",
    "  public.zero <- npreghat(bws=bw, txdat=x, y=Y, output='apply')",
    "  native <- local_regression(native_apply(bw, txdat=x, y=Y, degree=c(2L,1L), basis='glp', bernstein.basis=FALSE, s=c(0L,0L), return.hat=FALSE))",
    "  stopifnot(identical(unname(public.inf), unname(native)))",
    "  stopifnot(identical(unname(public.zero), unname(native)))",
    "  for (column in seq_len(ncol(Y))) {",
    "    fit <- npreg(bws=bw, txdat=x, tydat=Y[,column], warn.glp.gradient=FALSE)$mean",
    "    stopifnot(isTRUE(all.equal(public.inf[,column], fit, tolerance=1e-10)))",
    "  }",
    "}",
    sprintf("cat('%s\\n')", ok_tag)
  )
  res <- npRmpi_run_rscript_subprocess(lines=lines, timeout=90L, env=env)
  info <- paste(res$output, collapse="\n")
  expect_equal(res$status, 0L, info=info)
  expect_true(any(grepl(ok_tag, res$output, fixed=TRUE)), info=info)
})
