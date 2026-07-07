suppressPackageStartupMessages(library(npRmpi))
suppressPackageStartupMessages(library(testthat))

options(npRmpi.autodispatch = TRUE, np.messages = FALSE)
npRmpi.init(nslaves = 1L, quiet = TRUE)

# This fixture owns one MPI session for the selected npreghat contract block.
# The parent test wrapper runs each block in a subprocess so MPI lifecycle
# state cannot leak into the main testthat process or adjacent blocks.
spawn_mpi_slaves <- function(n = 1L) TRUE
close_mpi_slaves <- function(force = FALSE) invisible()

expect_try_error <- function(expr, regexp) {
  out <- try(force(expr), silent = TRUE)
  expect_true(inherits(out, "try-error"))
  expect_match(conditionMessage(attr(out, "condition")), regexp)
  invisible(out)
}

npreghat_fixture_selected <- Sys.getenv("NPREGHAT_FIXTURE_BLOCK", "all")
npreghat_fixture_index <- 0L
test_that <- function(desc, code) {
  npreghat_fixture_index <<- npreghat_fixture_index + 1L
  selected <- identical(npreghat_fixture_selected, "all") ||
    identical(npreghat_fixture_selected, as.character(npreghat_fixture_index))
  if (!selected)
    return(invisible(FALSE))

  cat(sprintf("NPREGHAT_FIXTURE_BLOCK_START %d %s\n",
              npreghat_fixture_index, desc))
  flush.console()
  force(code)
  cat(sprintf("NPREGHAT_FIXTURE_BLOCK_OK %d %s\n",
              npreghat_fixture_index, desc))
  flush.console()
  invisible(TRUE)
}

test_that("npreghat reproduces npreg fitted values for mixed-data local constant", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260223)
  n <- 120
  x <- runif(n)
  u <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  o <- ordered(sample(1:3, n, replace = TRUE))
  y <- sin(2 * pi * x) + as.numeric(u) / 3 + as.numeric(o) / 5 + rnorm(n, sd = 0.05)

  tx <- data.frame(x = x, u = u, o = o)
  bw <- npregbw(
    xdat = tx,
    ydat = y,
    bws = c(0.2, 0.4, 0.4),
    regtype = "lc",
    bandwidth.compute = FALSE
  )

  fit <- npreg(txdat = tx, tydat = y, bws = bw, warn.glp.gradient = FALSE)
  H <- npreghat(bws = bw, txdat = tx)

  expect_s3_class(H, "npreghat")
  expect_equal(as.vector(H %*% y), as.vector(fit$mean), tolerance = 1e-8)

  expect_identical(predict(H, output = "matrix"), H)
  H.loo <- npreghat(bws = bw, txdat = tx, leave.one.out = TRUE)
  expect_gt(max(abs(H - H.loo)), 1e-8)
  expect_lt(max(abs(diag(H.loo))), 1e-12)
  expect_error(
    npreghat(bws = bw, txdat = tx, exdat = tx[1:10, , drop = FALSE], leave.one.out = TRUE),
    "you may not specify 'leave.one.out = TRUE' and provide evaluation data"
  )
  hy <- predict(H.loo, y = y, output = "apply")
  expect_equal(as.vector(hy), as.vector(H.loo %*% y), tolerance = 1e-10)
  Y <- cbind(y, y + 0.1)
  expect_equal(predict(H.loo, y = Y, output = "apply"), H.loo %*% Y,
               tolerance = 1e-10)
  expect_error(
    predict(H.loo, output = "apply"),
    "argument 'y' is required when output='apply'"
  )

  ex <- tx[seq_len(10L), , drop = FALSE]
  H.ex <- npreghat(bws = bw, txdat = tx, exdat = ex, output = "matrix")
  expect_equal(predict(H, newdata = ex, output = "matrix"), H.ex,
               tolerance = 0, ignore_attr = TRUE)
})

test_that("npreghat supports lp/ll derivatives and matrix apply mode", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(777)
  n <- 150
  x <- sort(runif(n))
  y <- cos(2 * pi * x) + rnorm(n, sd = 0.03)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 40))

  bw <- npregbw(
    xdat = tx,
    ydat = y,
    bws = 0.2,
    regtype = "ll",
    bandwidth.compute = FALSE
  )

  fit <- npreg(
    txdat = tx,
    tydat = y,
    exdat = ex,
    bws = bw,
    gradients = TRUE,
    gradient.order = 1L,
    warn.glp.gradient = FALSE
  )

  H0 <- npreghat(bws = bw, txdat = tx, exdat = ex)
  H1 <- npreghat(bws = bw, txdat = tx, exdat = ex, s = 1L)

  expect_equal(as.vector(H0 %*% y), as.vector(fit$mean), tolerance = 1e-8)
  expect_equal(as.vector(H1 %*% y), as.vector(fit$grad[, 1]), tolerance = 1e-6)

  yboot <- cbind(y, y + 0.1)
  hy.apply <- npreghat(
    bws = bw,
    txdat = tx,
    exdat = ex,
    y = yboot,
    output = "apply"
  )

  expect_true(isTRUE(all.equal(
    hy.apply,
    H0 %*% yboot,
    tolerance = 1e-10,
    check.attributes = FALSE
  )))
})

test_that("npreghat lc derivative matrix matches analytic npreg for fixed and generalized-nn", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260311)
  n <- 70L
  x <- sort(runif(n))
  y <- sin(2 * pi * x) + 0.3 * x + rnorm(n, sd = 0.04)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(0.1, 0.9, length.out = 20L))

  for (bwtype in c("fixed", "generalized_nn")) {
    bw <- npregbw(
      xdat = tx,
      ydat = y,
      regtype = "lc",
      bwtype = bwtype,
      bandwidth.compute = FALSE,
      bws = if (identical(bwtype, "fixed")) 0.2 else 9L
    )

    fit <- npreg(
      bws = bw,
      txdat = tx,
      tydat = y,
      exdat = ex,
      gradients = TRUE,
      warn.glp.gradient = FALSE
    )
    H <- npreghat(bws = bw, txdat = tx, exdat = ex, output = "matrix", s = 1L)
    grad.apply <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = "apply", s = 1L)

    expect_equal(as.vector(H %*% y), as.vector(fit$grad[, 1]), tolerance = 1e-6)
    expect_equal(as.vector(H %*% y), as.vector(grad.apply), tolerance = 1e-8)
  }
})

test_that("npreghat lp bernstein path matches predict semantics", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260305)
  n <- 220
  x <- sort(runif(n))
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.08)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 60))

  bw.raw <- npregbw(
    xdat = tx,
    ydat = y,
    bws = 0.18,
    regtype = "lp",
    degree = 3L,
    bernstein = FALSE,
    bandwidth.compute = FALSE
  )
  bw.bern <- npregbw(
    xdat = tx,
    ydat = y,
    bws = 0.18,
    regtype = "lp",
    degree = 3L,
    bernstein = TRUE,
    bandwidth.compute = FALSE
  )

  fit.raw <- npreg(txdat = tx, tydat = y, exdat = ex, bws = bw.raw, gradients = FALSE, warn.glp.gradient = FALSE)
  fit.bern <- npreg(txdat = tx, tydat = y, exdat = ex, bws = bw.bern, gradients = FALSE, warn.glp.gradient = FALSE)
  grad.raw <- npreg(txdat = tx, tydat = y, exdat = ex, bws = bw.raw, gradients = TRUE, warn.glp.gradient = FALSE)
  grad.bern <- npreg(txdat = tx, tydat = y, exdat = ex, bws = bw.bern, gradients = TRUE, warn.glp.gradient = FALSE)

  hat.raw <- npreghat(bws = bw.raw, txdat = tx, exdat = ex, y = y, output = "apply", s = 0L)
  hat.bern <- npreghat(bws = bw.bern, txdat = tx, exdat = ex, y = y, output = "apply", s = 0L)
  dhat.raw <- npreghat(bws = bw.raw, txdat = tx, exdat = ex, y = y, output = "apply", s = 1L)
  dhat.bern <- npreghat(bws = bw.bern, txdat = tx, exdat = ex, y = y, output = "apply", s = 1L)

  expect_equal(as.vector(hat.raw), as.vector(fit.raw$mean), tolerance = 1e-8)
  expect_equal(as.vector(hat.bern), as.vector(fit.bern$mean), tolerance = 1e-8)
  expect_equal(as.vector(fit.bern$mean), as.vector(fit.raw$mean), tolerance = 1e-8)
  expect_equal(as.vector(dhat.raw), as.vector(grad.raw$grad[, 1L]), tolerance = 1e-8)
  expect_equal(as.vector(dhat.bern), as.vector(grad.bern$grad[, 1L]), tolerance = 1e-8)
  expect_equal(as.vector(grad.bern$grad[, 1L]), as.vector(grad.raw$grad[, 1L]), tolerance = 1e-8)
})

test_that("npreghat apply mode matches npreg across bwtypes", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260308)
  n <- 90
  x <- sort(runif(n))
  y <- sin(2 * pi * x) + 0.25 * x + rnorm(n, sd = 0.04)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(0.05, 0.95, length.out = 25))

  make_bw <- function(regtype, bwtype) {
    bw_args <- list(
      xdat = tx,
      ydat = y,
      regtype = regtype,
      bwtype = bwtype,
      bandwidth.compute = FALSE,
      bws = if (identical(bwtype, "fixed")) 0.18 else 9
    )
    if (identical(regtype, "lp")) {
      bw_args$degree <- 1L
      bw_args$basis <- "glp"
      bw_args$bernstein.basis <- FALSE
    }
    do.call(npregbw, bw_args)
  }

  for (regtype in c("lc", "ll", "lp")) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      bw <- make_bw(regtype, bwtype)

      fit.in <- npreg(
        bws = bw,
        txdat = tx,
        tydat = y,
        gradients = TRUE,
        warn.glp.gradient = FALSE
      )
      fit.ex <- npreg(
        bws = bw,
        txdat = tx,
        tydat = y,
        exdat = ex,
        gradients = TRUE,
        warn.glp.gradient = FALSE
      )
      pred.ex <- predict(fit.in, newdata = ex)

      hat.in.mean <- npreghat(bws = bw, txdat = tx, y = y, output = "apply")
      hat.ex.mean <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = "apply")

      expect_equal(
        as.vector(hat.in.mean),
        as.vector(fit.in$mean),
        tolerance = 1e-8,
        info = paste("mean in-sample", regtype, bwtype)
      )
      expect_equal(
        as.vector(hat.ex.mean),
        as.vector(fit.ex$mean),
        tolerance = 1e-8,
        info = paste("mean exdat", regtype, bwtype)
      )
      expect_equal(
        as.vector(hat.ex.mean),
        as.vector(pred.ex),
        tolerance = 1e-8,
        info = paste("mean predict", regtype, bwtype)
      )

      hat.in.grad <- npreghat(bws = bw, txdat = tx, y = y, output = "apply", s = 1L)
      hat.ex.grad <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = "apply", s = 1L)

      expect_equal(
        as.vector(hat.in.grad),
        as.vector(fit.in$grad[, 1]),
        tolerance = 1e-6,
        info = paste("grad in-sample", regtype, bwtype)
      )
      expect_equal(
        as.vector(hat.ex.grad),
        as.vector(fit.ex$grad[, 1]),
        tolerance = 1e-6,
        info = paste("grad exdat", regtype, bwtype)
      )
    }
  }
})

test_that("npreghat nonfixed higher-order lp operator matches npreg and matrix apply semantics", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  make_case <- function(seed, bwtype) {
    set.seed(seed)
    n <- 28
    x <- data.frame(x = runif(n))
    y <- cos(2 * pi * x$x) + rnorm(n, sd = if (identical(bwtype, "generalized_nn")) 0.05 else 0.03)
    Y <- cbind(y, y^2)
    bw <- npregbw(
      xdat = x,
      ydat = y,
      regtype = "lp",
      degree = 2,
      basis = "tensor",
      bernstein.basis = TRUE,
      bwmethod = "cv.ls",
      bwtype = bwtype
    )
    fit <- npreg(bws = bw, txdat = x, tydat = y)
    H <- npreghat(bws = bw, txdat = x, output = "matrix")
    a.vec <- npreghat(bws = bw, txdat = x, y = y, output = "apply")
    a.mat <- npreghat(bws = bw, txdat = x, y = Y, output = "apply")
    list(fit = fit, H = H, a.vec = a.vec, a.mat = a.mat, y = y, Y = Y)
  }

  for (case in list(
    make_case(3, "generalized_nn"),
    make_case(5, "adaptive_nn")
  )) {
    expect_equal(as.vector(case$H %*% case$y), as.vector(case$fit$mean), tolerance = 1e-8)
    expect_equal(as.vector(case$a.vec), as.vector(case$fit$mean), tolerance = 1e-8)
    expect_equal(as.vector(case$H %*% case$y), as.vector(case$a.vec), tolerance = 1e-10)
    expect_equal(case$H %*% case$Y, case$a.mat, tolerance = 1e-10, ignore_attr = TRUE)
  }
})

test_that("npreghat constraint output is exact row-weighted transpose", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260513)
  n <- 24L
  x <- sort(runif(n))
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.05)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(0.05, 0.95, length.out = 6L))

  make_bw <- function(regtype, bwtype) {
    args <- list(
      xdat = tx,
      ydat = y,
      regtype = regtype,
      bwtype = bwtype,
      bandwidth.compute = FALSE,
      bws = if (identical(bwtype, "fixed")) 0.25 else 7L
    )
    if (identical(regtype, "lp")) {
      args$degree <- 2L
      args$basis <- "glp"
      args$bernstein.basis <- FALSE
    }
    do.call(npregbw, args)
  }

  for (regtype in c("lc", "ll", "lp")) {
    bwtype <- "fixed"
    bw <- make_bw(regtype, bwtype)
    H <- npreghat(bws = bw, txdat = tx, exdat = ex, output = "matrix")
    A <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = "constraint")
    expect_equal(
      A,
      t(H) * y,
      tolerance = 0,
      ignore_attr = TRUE,
      info = paste("mean", regtype, bwtype)
    )

    if (!identical(regtype, "lc")) {
      H.grad <- npreghat(bws = bw, txdat = tx, exdat = ex, output = "matrix", s = 1L)
      A.grad <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y,
                         output = "constraint", s = 1L)
      expect_equal(
        A.grad,
        t(H.grad) * y,
        tolerance = 1e-14,
        ignore_attr = TRUE,
        info = paste("gradient", regtype, bwtype)
      )
    }
  }

  bw <- make_bw("ll", "fixed")
  H <- npreghat(bws = bw, txdat = tx, exdat = ex, output = "matrix")
  H.obj <- npreghat(bws = bw, txdat = tx, output = "matrix")
  A.stored <- predict(H, y = y, output = "constraint")
  expect_equal(A.stored, t(H) * y, tolerance = 1e-14, ignore_attr = TRUE)
  expect_try_error(
    predict(H, output = "constraint"),
    "argument 'y' is required"
  )
  A.pred <- predict(H.obj, newdata = ex, y = y, output = "constraint")
  expect_equal(A.pred, t(H) * y, tolerance = 1e-14, ignore_attr = TRUE)
})

run_npreghat_constraint_error_case <- function(case, regexp) {
  selected <- identical(npreghat_fixture_selected, "all") ||
    identical(npreghat_fixture_selected, as.character(case))
  if (!selected)
    return(invisible(FALSE))

  set.seed(20260513)
  n <- 24L
  x <- sort(runif(n))
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.05)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(0.05, 0.95, length.out = 6L))
  bw <- npregbw(
    xdat = tx,
    ydat = y,
    regtype = "ll",
    bwtype = "fixed",
    bandwidth.compute = FALSE,
    bws = 0.25
  )

  cat(sprintf("NPREGHAT_FIXTURE_BLOCK_START %d constraint error\n", case))
  flush.console()
  if (case == 8L) {
    expect_try_error(
      npreghat(bws = bw, txdat = tx, exdat = ex, output = "constraint"),
      regexp
    )
  } else if (case == 9L) {
    expect_try_error(
      npreghat(bws = bw, txdat = tx, exdat = ex, y = cbind(y, y),
               output = "constraint"),
      regexp
    )
  } else if (case == 10L) {
    expect_try_error(
      npreghat(bws = bw, txdat = tx, exdat = ex, y = y[-1L],
               output = "constraint"),
      regexp
    )
  }
  cat(sprintf("NPREGHAT_FIXTURE_BLOCK_OK %d constraint error\n", case))
  flush.console()
  invisible(TRUE)
}

run_npreghat_constraint_error_case(8L, "argument 'y' is required")
run_npreghat_constraint_error_case(9L, "one-column")
run_npreghat_constraint_error_case(10L, "must match")

cat("NPREGHAT_FULL_FILE_SUBPROCESS_OK\n")
try(npRmpi.quit(force = TRUE), silent = TRUE)
