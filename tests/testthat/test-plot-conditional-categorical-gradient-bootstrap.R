test_that("conditional categorical gradient bootstrap helper matches explicit local refits", {
  skip_if_not_installed("np")

  library(np)

  helper <- getFromNamespace(".np_inid_boot_from_conditional_gradient_local", "np")
  eval_fun <- getFromNamespace(".np_plot_conditional_eval", "np")

  set.seed(20260312)
  n <- 24L
  x1 <- factor(sample(c("a", "b"), n, replace = TRUE))
  x2 <- rnorm(n)
  y <- rnorm(n)
  xdat <- data.frame(x1 = x1, x2 = x2)
  ydat <- data.frame(y = y)

  fixtures <- list(
    list(
      label = "dens",
      cdf = FALSE,
      bw = suppressWarnings(npcdensbw(xdat = xdat, ydat = ydat, nmulti = 1L))
    ),
    list(
      label = "dist",
      cdf = TRUE,
      bw = suppressWarnings(npcdistbw(xdat = xdat, ydat = ydat, nmulti = 1L))
    )
  )

  counts <- cbind(
    rep(1L, n),
    c(rep(2L, 4), rep(0L, 4), rep(1L, n - 8)),
    c(rep(0L, 3), rep(3L, 3), rep(1L, n - 6))
  )
  if (!is.matrix(counts))
    counts <- matrix(counts, nrow = n)

  for (fixture in fixtures) {
    slice <- suppressWarnings(plot(
      fixture$bw,
      xdat = xdat,
      ydat = ydat,
      gradients = TRUE,
      perspective = FALSE,
      output = "data"
    ))[[1L]]

    boot <- helper(
      xdat = xdat,
      ydat = ydat,
      exdat = slice$xeval,
      eydat = slice$yeval,
      bws = fixture$bw,
      B = ncol(counts),
      cdf = fixture$cdf,
      gradient.index = 1L,
      counts = counts
    )

    explicit_fit <- function(idx) {
      as.vector(eval_fun(
        bws = fixture$bw,
        xdat = xdat[idx, , drop = FALSE],
        ydat = ydat[idx, , drop = FALSE],
        exdat = slice$xeval,
        eydat = slice$yeval,
        cdf = fixture$cdf,
        gradients = TRUE
      )$congrad[, 1L])
    }

    oracle.t0 <- explicit_fit(seq_len(n))
    oracle.t <- vapply(
      seq_len(ncol(counts)),
      function(j) explicit_fit(rep.int(seq_len(n), counts[, j])),
      numeric(length(oracle.t0))
    )

    expect_equal(as.vector(boot$t0), oracle.t0, tolerance = 1e-8, info = fixture$label)
    expect_equal(boot$t, t(oracle.t), tolerance = 1e-8, info = fixture$label)
  }
})

test_that("conditional density/distribution categorical bootstrap gradients work again", {
  skip_if_not_installed("np")

  library(np)

  set.seed(20260312)
  n <- 24L
  x1 <- factor(sample(c("a", "b"), n, replace = TRUE))
  x2 <- rnorm(n)
  y <- rnorm(n)
  xdat <- data.frame(x1 = x1, x2 = x2)
  ydat <- data.frame(y = y)

  cd.bw <- suppressWarnings(npcdensbw(xdat = xdat, ydat = ydat, nmulti = 1L))
  cdist.bw <- suppressWarnings(npcdistbw(xdat = xdat, ydat = ydat, nmulti = 1L))
  cases <- list(
    list(label = "npcdens-bw", obj = cd.bw, args = list(xdat = xdat, ydat = ydat)),
    list(label = "npcdens-fit", obj = npcdens(bws = cd.bw), args = list()),
    list(label = "npcdist-bw", obj = cdist.bw, args = list(xdat = xdat, ydat = ydat)),
    list(label = "npcdist-fit", obj = npcdist(bws = cdist.bw), args = list())
  )

  for (case in cases) {
    for (boot.method in c("inid", "fixed", "geom")) {
      plot.args <- c(
        list(
          case$obj,
          perspective = FALSE,
          output = "data",
          gradients = TRUE,
          errors = "bootstrap",
          bootstrap = boot.method,
          B = 7L
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
