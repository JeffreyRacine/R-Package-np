library(npRmpi)

test_that("generalized-nn degree-1 raw and Bernstein fits share canonical sandwich SEs", {
  skip_on_cran()
  if (!spawn_mpi_slaves())
    skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260723)

  n <- 193L
  tx <- data.frame(
    x1 = runif(n, -0.75, 1.15),
    x2 = runif(n, -0.60, 1.30)
  )
  y <- 0.4 + sin(1.7 * tx$x1) - 0.5 * cos(2.1 * tx$x2) +
    0.25 * tx$x1 * tx$x2 + rnorm(n, sd = 0.16)
  ex <- tx[seq.int(3L, n, by = 4L), , drop = FALSE]

  make_bw <- function(regtype, bernstein = FALSE) {
    args <- list(
      xdat = tx,
      ydat = y,
      regtype = regtype,
      bwmethod = "cv.ls",
      bwtype = "generalized_nn",
      ckertype = "gaussian",
      ckerorder = 2L,
      bws = c(37L, 41L),
      bandwidth.compute = FALSE
    )
    if (identical(regtype, "lp")) {
      args$degree <- c(1L, 1L)
      args$degree.select <- "manual"
      args$basis <- "glp"
      args$bernstein.basis <- bernstein
    }
    do.call(npregbw, args)
  }

  fit <- function(bw) {
    npreg(
      bws = bw,
      txdat = tx,
      tydat = y,
      exdat = ex,
      gradients = TRUE,
      warn.glp.gradient = FALSE
    )
  }

  bw.ll <- make_bw("ll")
  bw.raw <- make_bw("lp", bernstein = FALSE)
  bw.bernstein <- make_bw("lp", bernstein = TRUE)

  fit.ll <- fit(bw.ll)
  fit.raw <- fit(bw.raw)
  fit.bernstein <- fit(bw.bernstein)

  for (surface in c("mean", "merr", "grad", "gerr")) {
    expect_equal(
      fit.ll[[surface]],
      fit.raw[[surface]],
      tolerance = 1e-10,
      info = paste("ll versus raw", surface)
    )
    expect_equal(
      fit.raw[[surface]],
      fit.bernstein[[surface]],
      tolerance = 1e-10,
      info = paste("raw versus Bernstein", surface)
    )
  }

  hat.raw <- npreghat(
    bws = bw.raw, txdat = tx, exdat = ex, y = y, output = "apply"
  )
  hat.bernstein <- npreghat(
    bws = bw.bernstein, txdat = tx, exdat = ex, y = y, output = "apply"
  )
  expect_equal(hat.raw, hat.bernstein, tolerance = 1e-10)
})
