library(npRmpi)

with_master_only <- function(expr) {
  if (isTRUE(getOption("npRmpi.mpi.initialized", FALSE))) {
    try(npRmpi.quit(mode = "spawn", force = TRUE), silent = TRUE)
  }
  suppressWarnings(npRmpi.init(nslaves = 0, quiet = TRUE))
  on.exit(try(npRmpi.quit(mode = "spawn", force = TRUE), silent = TRUE), add = TRUE)
  force(expr)
}

test_that("direct kernel-weight helper matches npksum.default with evaluation data", {
  skip_on_cran()
  with_master_only({
    kw.fun <- getFromNamespace(".np_kernel_weights_direct", "npRmpi")
    ksum.fun <- getFromNamespace("npksum.default", "npRmpi")
    set.seed(4021)
    n <- 60
    x <- runif(n)
    z <- rnorm(n)
    tx <- data.frame(x = x, z = z)
    ex <- data.frame(
      x = seq(min(x), max(x), length.out = 21),
      z = seq(min(z), max(z), length.out = 21)
    )
    y <- sin(x) + 0.1 * z + rnorm(n, sd = 0.05)

    bw <- npregbw(
      xdat = tx,
      ydat = y,
      regtype = "ll",
      bws = c(0.3, 0.4),
      bandwidth.compute = FALSE
    )

    kw.direct <- kw.fun(
      bws = bw,
      txdat = tx,
      exdat = ex,
      leave.one.out = FALSE,
      bandwidth.divide = TRUE
    )

    kw.ref <- ksum.fun(
      bws = bw,
      txdat = tx,
      exdat = ex,
      return.kernel.weights = TRUE,
      bandwidth.divide = TRUE
    )$kw

    expect_equal(kw.direct, kw.ref, tolerance = 1e-12)
  })
})

test_that("direct kernel-weight helper matches npksum.default leave-one-out", {
  skip_on_cran()
  with_master_only({
    kw.fun <- getFromNamespace(".np_kernel_weights_direct", "npRmpi")
    ksum.fun <- getFromNamespace("npksum.default", "npRmpi")
    set.seed(4022)
    n <- 50
    x <- runif(n)
    y <- cos(2 * pi * x) + rnorm(n, sd = 0.04)
    tx <- data.frame(x = x)

    bw <- npregbw(
      xdat = tx,
      ydat = y,
      regtype = "lc",
      bws = 0.25,
      bandwidth.compute = FALSE
    )

    kw.direct <- kw.fun(
      bws = bw,
      txdat = tx,
      exdat = NULL,
      leave.one.out = TRUE,
      bandwidth.divide = TRUE
    )

    kw.ref <- ksum.fun(
      bws = bw,
      txdat = tx,
      leave.one.out = TRUE,
      return.kernel.weights = TRUE,
      bandwidth.divide = TRUE
    )$kw
    if (nrow(kw.ref) == ncol(kw.ref))
      diag(kw.ref) <- 0

    expect_equal(kw.direct, kw.ref, tolerance = 1e-12)
  })
})
