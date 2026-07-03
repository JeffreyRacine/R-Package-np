local_npcdens_audit_bw <- function(x, y, ybw = 0.28, xbw = 0.22,
                                   regtype = "lp", degree = 1L) {
  npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(ybw, xbw),
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    regtype = regtype,
    degree = degree
  )
}

local_npcdens_audit_grid <- function(xgrid, ygrid, jitter = 0) {
  out <- do.call(rbind, lapply(xgrid, function(xx) {
    jj <- if (jitter > 0) seq(-jitter, jitter, length.out = length(ygrid)) else 0
    data.frame(y = ygrid, x = xx + jj)
  }))
  row.names(out) <- NULL
  out
}

test_that("npcdens formula reentry honors explicit data argument", {
  set.seed(1101)
  d1 <- data.frame(x = seq(-1, 1, length.out = 44L))
  d1$y <- 0.25 + d1$x + rnorm(nrow(d1), sd = 0.2)
  d2 <- data.frame(x = seq(2, 4, length.out = 31L))
  d2$y <- -1 + 0.5 * d2$x + rnorm(nrow(d2), sd = 0.2)

  bw <- npcdensbw(
    y ~ x,
    data = d1,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    bwtype = "fixed"
  )

  fit <- npcdens(bws = bw, data = d2)
  direct <- npcdens(bws = bw, txdat = d2["x"], tydat = d2["y"])

  expect_equal(fit$ntrain, nrow(d2))
  expect_equal(fitted(fit), fitted(direct), tolerance = 1e-12)
})

test_that("npcdens proper repair warns on under-covering grids and tolerates numeric x jitter", {
  set.seed(2202)
  n <- 80L
  x <- runif(n, -1, 1)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.20)
  bw <- local_npcdens_audit_bw(x, y, ybw = 0.28, xbw = 0.22, regtype = "lp", degree = 3L)
  xgrid <- unname(quantile(x, probs = c(0.25, 0.75)))
  yfull <- seq(min(y) - 0.5, max(y) + 0.5, length.out = 80L)
  ynarrow <- seq(quantile(y, 0.40), quantile(y, 0.60), length.out = 50L)

  full.grid <- local_npcdens_audit_grid(xgrid, yfull)
  jitter.grid <- local_npcdens_audit_grid(xgrid, yfull, jitter = 1e-12)
  narrow.grid <- local_npcdens_audit_grid(xgrid, ynarrow)

  fit.full <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = full.grid["x"],
    eydat = full.grid["y"],
    proper = TRUE,
    proper.control = list(tol = 1e-8, mass.warn.tol = 0)
  )
  fit.jitter <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = jitter.grid["x"],
    eydat = jitter.grid["y"],
    proper = TRUE,
    proper.control = list(tol = 1e-8, mass.warn.tol = 0)
  )

  expect_true(isTRUE(fit.full$proper.applied))
  expect_true(isTRUE(fit.jitter$proper.applied))
  expect_equal(fit.jitter$proper.info$slice.count, fit.full$proper.info$slice.count)
  expect_equal(fitted(fit.jitter), fitted(fit.full), tolerance = 1e-7)

  expect_warning(
    fit.narrow <- npcdens(
      bws = bw,
      txdat = data.frame(x = x),
      tydat = data.frame(y = y),
      exdat = narrow.grid["x"],
      eydat = narrow.grid["y"],
      proper = TRUE
    ),
    "under-cover response support"
  )

  expect_true(isTRUE(fit.narrow$proper.applied))
  expect_true(any(fit.narrow$proper.info$integral.raw < 0.95))
})

test_that("npcdens records train and evaluation omit metadata without changing legacy rows.omit", {
  tx <- data.frame(x = c(0.0, 0.2, NA, 0.5, 0.9, 1.1))
  ty <- data.frame(y = c(1, 2, 3, 4, 5, NA))
  ex <- data.frame(x = c(0.1, NA, 0.7, 0.8))
  ey <- data.frame(y = c(1, 2, NA, 4))

  bw <- npcdensbw(
    xdat = tx,
    ydat = ty,
    bws = c(0.4, 0.4),
    bandwidth.compute = FALSE,
    bwtype = "fixed"
  )
  fit <- npcdens(bws = bw, txdat = tx, tydat = ty, exdat = ex, eydat = ey)

  expect_equal(unname(fit$rows.omit), c(2L, 3L))
  expect_equal(unname(fit$train.rows.omit), c(3L, 6L))
  expect_equal(unname(fit$eval.rows.omit), c(2L, 3L))
  expect_equal(fit$train.nobs.omit, 2L)
  expect_equal(fit$eval.nobs.omit, 2L)
})

test_that("npcdens bounded cv.ls method metadata and infinite-bound surrogate are explicit", {
  set.seed(3303)
  x <- data.frame(x = runif(32L))
  y <- data.frame(y = 0.25 + x$x + rchisq(32L, df = 3))

  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    regtype = "lp",
    degree = 0,
    cxkerbound = "fixed",
    cxkerlb = 0,
    cxkerub = 1,
    cykerbound = "fixed",
    cykerlb = 0,
    cykerub = Inf
  )
  bw$method <- NULL

  guard <- getFromNamespace(".npcdensbw_assert_bounded_cvls_supported", "np")
  expect_error(guard(bw, where = "probe"), "method metadata")

  marshal <- getFromNamespace(".npcdensbw_marshal_y_bounds", "np")
  bounded <- marshal(-Inf, Inf, "fixed")
  expect_equal(bounded$lb, -1e300)
  expect_equal(bounded$ub, 1e300)
})
