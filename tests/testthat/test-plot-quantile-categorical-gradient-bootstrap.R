test_that("quantile categorical gradient bootstrap helper matches explicit local refits", {
  skip_if_not_installed("np")

  library(np)

  helper <- getFromNamespace(".np_inid_boot_from_quantile_gradient_local", "np")
  eval_fun <- getFromNamespace(".np_plot_quantile_eval", "np")

  set.seed(20260312)
  n <- 24L
  x1 <- factor(sample(c("a", "b"), n, replace = TRUE))
  x2 <- rnorm(n)
  y <- rnorm(n)
  xdat <- data.frame(x1 = x1, x2 = x2)
  ydat <- data.frame(y = y)
  bw <- suppressWarnings(npcdistbw(xdat = xdat, ydat = ydat, nmulti = 1L))
  tau <- 0.5

  slice <- suppressWarnings(plot(
    bw,
    xdat = xdat,
    ydat = ydat,
    tau = tau,
    gradients = TRUE,
    perspective = FALSE,
    output = "data"
  ))[[1L]]

  counts <- cbind(
    rep(1L, n),
    c(rep(2L, 4), rep(0L, 4), rep(1L, n - 8)),
    c(rep(0L, 3), rep(3L, 3), rep(1L, n - 6))
  )
  if (!is.matrix(counts))
    counts <- matrix(counts, nrow = n)

  boot <- helper(
    xdat = xdat,
    ydat = y,
    exdat = slice$xeval,
    bws = bw,
    B = ncol(counts),
    tau = tau,
    gradient.index = 1L,
    counts = counts
  )

  explicit_fit <- function(idx) {
    as.vector(eval_fun(
      bws = bw,
      txdat = xdat[idx, , drop = FALSE],
      tydat = y[idx],
      exdat = slice$xeval,
      tau = tau,
      gradients = TRUE
    )$quantgrad[, 1L])
  }

  oracle.t0 <- explicit_fit(seq_len(n))
  oracle.t <- vapply(
    seq_len(ncol(counts)),
    function(j) explicit_fit(rep.int(seq_len(n), counts[, j])),
    numeric(length(oracle.t0))
  )

  expect_equal(as.vector(boot$t0), oracle.t0, tolerance = 1e-8)
  expect_equal(boot$t, t(oracle.t), tolerance = 1e-8)
})

test_that("quantile categorical bootstrap gradients work again", {
  skip_if_not_installed("np")

  library(np)

  set.seed(20260312)
  n <- 24L
  x1 <- factor(sample(c("a", "b"), n, replace = TRUE))
  x2 <- rnorm(n)
  y <- rnorm(n)
  xdat <- data.frame(x1 = x1, x2 = x2)
  ydat <- data.frame(y = y)
  bw <- suppressWarnings(npcdistbw(xdat = xdat, ydat = ydat, nmulti = 1L))
  fit <- npqreg(bws = bw, txdat = xdat, tydat = y, tau = 0.5)

  for (case in list(
    list(label = "bw", obj = bw, args = list(xdat = xdat, ydat = ydat)),
    list(label = "fit", obj = fit, args = list())
  )) {
    for (boot.method in c("inid", "fixed", "geom")) {
      plot.args <- c(
        list(
          case$obj,
          tau = 0.5,
          gradients = TRUE,
          perspective = FALSE,
          output = "data",
          errors = "bootstrap",
          bootstrap = boot.method,
          B = 5L
        ),
        case$args
      )
      out <- suppressWarnings(do.call(plot, plot.args))
      expect_true(is.list(out), info = paste(case$label, boot.method))
      expect_true(length(out) >= 1L, info = paste(case$label, boot.method))
      expect_true(length(out[[1L]]$bxp) > 0L, info = paste(case$label, boot.method))
      expect_equal(length(out[[1L]]$bxp$names), 2L, info = paste(case$label, boot.method))
    }
  }
})
