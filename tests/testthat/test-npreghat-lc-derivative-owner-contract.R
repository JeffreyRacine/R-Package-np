library(np)

test_that("lc first-derivative owner matches npreg and apply on common cells", {
  run_case <- function(bwtype, bwval, tol = 1e-10) {
    train <- data.frame(x = x)
    frame <- data.frame(train, y = y)
    bw <- npregbw(y ~ x,
                  data = frame,
                  regtype = "lc",
                  bwtype = bwtype,
                  bandwidth.compute = FALSE,
                  bws = bwval)

    H.train <- npreghat(bws = bw, txdat = train, s = 1L)
    H.eval <- npreghat(bws = bw, txdat = train, exdat = eval, s = 1L)
    a.train <- npreghat(bws = bw, txdat = train, y = y, s = 1L, output = "apply")
    a.eval <- npreghat(bws = bw, txdat = train, exdat = eval, y = y, s = 1L, output = "apply")
    g.train <- npreg(bws = bw, gradients = TRUE)
    g.eval <- npreg(bws = bw, exdat = eval, gradients = TRUE)

    expect_equal(drop(H.train %*% y), a.train, tolerance = tol)
    expect_equal(drop(H.eval %*% y), a.eval, tolerance = tol)
    expect_equal(drop(H.train %*% y), g.train$grad[, 1L], tolerance = tol)
    expect_equal(drop(H.eval %*% y), g.eval$grad[, 1L], tolerance = tol)
  }

  set.seed(410)
  n <- 60L
  x <- rnorm(n)
  y <- x^2 + rnorm(n, sd = 0.2)
  eval <- data.frame(x = seq(min(x), max(x), length.out = 31L))

  run_case("fixed", 0.45)
  run_case("generalized_nn", 9)
  run_case("adaptive_nn", 9)
})

test_that("mixed lc first-derivative owner supports matrix, apply, and constraint", {
  set.seed(20260630)
  n <- 36L
  idx <- seq_len(n)
  tx <- data.frame(
    x1 = seq(0.08, 0.92, length.out = n),
    x2 = sin(seq(0.2, 2.8, length.out = n)) + idx / 200,
    u = factor(rep(c("a", "b", "c"), length.out = n), levels = c("a", "b", "c")),
    o = ordered(rep(c("low", "mid", "high"), length.out = n),
                levels = c("low", "mid", "high"))
  )
  ex <- tx[seq_len(14L), , drop = FALSE]
  ex$x1 <- pmin(0.98, pmax(0.02, ex$x1 + 0.015))
  ex$x2 <- pmin(0.98, pmax(0.02, ex$x2 + 0.015))

  y <- 1.0 + sin(2 * pi * tx$x1) + 0.25 * tx$x2^2 +
    0.12 * as.integer(tx$u) - 0.08 * as.integer(tx$o)
  Y <- cbind(y = y, y2 = y^2)
  s <- c(1L, 0L)

  run_case <- function(bwtype, bws, eval_data = NULL, tol = 1e-8) {
    bw <- npregbw(
      xdat = tx,
      ydat = y,
      bws = bws,
      bandwidth.compute = FALSE,
      bwtype = bwtype,
      regtype = "lc"
    )

    base_args <- list(bws = bw, txdat = tx, s = s)
    if (!is.null(eval_data))
      base_args$exdat <- eval_data

    H <- do.call(npreghat, c(base_args, list(output = "matrix")))
    A1 <- do.call(npreghat, c(base_args, list(y = y, output = "apply")))
    A2 <- do.call(npreghat, c(base_args, list(y = Y, output = "apply")))
    C <- do.call(npreghat, c(base_args, list(y = y, output = "constraint")))

    expect_true(is.matrix(H))
    expect_true(all(is.finite(as.matrix(H))))
    expect_equal(ncol(H), nrow(tx))
    expect_equal(nrow(H), if (is.null(eval_data)) nrow(tx) else nrow(eval_data))
    Hmat <- matrix(as.numeric(H), nrow = nrow(H), ncol = ncol(H))
    expect_equal(as.vector(A1), as.vector(Hmat %*% y), tolerance = tol)
    expect_equal(unname(as.matrix(A2)), unname(Hmat %*% Y), tolerance = tol)
    expect_equal(unname(as.matrix(C)), unname(t(Hmat) * as.vector(y)), tolerance = tol)
  }

  for (eval_data in list(NULL, ex)) {
    run_case("fixed", c(0.55, 0.55, 0.45, 0.45), eval_data = eval_data)
    run_case("generalized_nn", c(7, 7, 0.45, 0.45), eval_data = eval_data)
    run_case("adaptive_nn", c(7, 7, 0.45, 0.45), eval_data = eval_data)
  }

  bw <- npregbw(
    xdat = tx[, c("x1", "x2")],
    ydat = y,
    bws = c(0.55, 0.55),
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    regtype = "lc"
  )
  expect_error(
    npreghat(bws = bw, txdat = tx[, c("x1", "x2")], s = c(2L, 0L)),
    "requested derivative order in 's' exceeds local polynomial degree"
  )
})
