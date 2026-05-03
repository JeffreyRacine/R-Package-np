library(npRmpi)

npreg_ll_lp1_payload <- function(fit) {
  list(
    mean = as.numeric(fit$mean),
    merr = as.numeric(fit$merr),
    grad = as.numeric(fit$grad),
    gerr = as.numeric(fit$gerr)
  )
}

npreg_expect_payload_equal <- function(a, b, tolerance = 1e-10) {
  pa <- npreg_ll_lp1_payload(a)
  pb <- npreg_ll_lp1_payload(b)
  for (nm in names(pa)) {
    expect_equal(pa[[nm]], pb[[nm]], tolerance = tolerance, info = nm)
  }
}

test_that("npreg ll uses the lp degree-one engine without public metadata drift", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old.opts <- options(npRmpi.autodispatch = TRUE,
                      np.messages = FALSE,
                      np.tree = FALSE)
  on.exit(options(old.opts), add = TRUE)

  set.seed(20260503)
  n <- 90L
  d <- data.frame(
    y = rnorm(n),
    x1 = runif(n),
    x2 = runif(n),
    u = factor(sample(letters[1:3], n, replace = TRUE)),
    o = ordered(sample(1:3, n, replace = TRUE))
  )
  d$y <- sin(2 * pi * d$x1) + d$x2 + as.numeric(d$u) / 10 + rnorm(n, sd = 0.04)
  ex <- data.frame(
    x1 = runif(25L),
    x2 = runif(25L),
    u = factor(sample(levels(d$u), 25L, replace = TRUE), levels = levels(d$u)),
    o = ordered(sample(levels(d$o), 25L, replace = TRUE), levels = levels(d$o))
  )

  specs <- list(
    fixed = list(bwtype = "fixed", bws = c(0.24, 0.24, 0.35, 0.35)),
    generalized_nn = list(bwtype = "generalized_nn", bws = c(12, 12, 0.35, 0.35))
  )

  for (spec in specs) {
    common <- list(
      formula = y ~ x1 + x2 + u + o,
      data = d,
      bwtype = spec$bwtype,
      bwmethod = "cv.ls",
      bws = spec$bws,
      bandwidth.compute = FALSE
    )
    bw_ll <- do.call(npregbw, c(common, list(regtype = "ll")))
    bw_lp <- do.call(npregbw, c(common, list(
      regtype = "lp",
      basis = "glp",
      degree = c(1L, 1L),
      bernstein.basis = FALSE
    )))

    expect_identical(bw_ll$regtype, "ll")
    expect_null(bw_ll$degree)
    expect_identical(bw_ll$regtype.engine, "lp")
    expect_identical(bw_ll$basis.engine, "glp")
    expect_equal(bw_ll$degree.engine, c(1L, 1L))
    expect_false(isTRUE(bw_ll$bernstein.basis.engine))

    fit_ll <- npreg(bws = bw_ll, gradients = TRUE)
    fit_lp <- npreg(bws = bw_lp, gradients = TRUE)
    ex_ll <- npreg(bws = bw_ll, exdat = ex, gradients = TRUE)
    ex_lp <- npreg(bws = bw_lp, exdat = ex, gradients = TRUE)

    npreg_expect_payload_equal(fit_ll, fit_lp)
    npreg_expect_payload_equal(ex_ll, ex_lp)
  }
})
