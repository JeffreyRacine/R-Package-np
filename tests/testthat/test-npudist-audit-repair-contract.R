with_npudist_audit_runtime <- function(code) {
  code <- substitute(code)
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  if (exists("spawn_mpi_slaves", mode = "function")) {
    skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }

  eval(code, envir = parent.frame())
}

expect_omit <- function(x, value) {
  expect_identical(as.integer(x), as.integer(value))
}

test_that("npudistbw CDF grid helper preserves public grid modes", {
  with_npudist_audit_runtime({
    dat <- data.frame(
      x1 = seq(-1, 1, length.out = 9L),
      x2 = c(-2, -1.5, -0.5, 0, 0.25, 0.5, 1.2, 1.8, 2.5)
    )
    bw <- npudistbw(dat = dat, bws = c(0.4, 0.5),
                    bandwidth.compute = FALSE)
    probs <- seq(0, 1, length.out = 7L)

    grid <- .npudistbw_prepare_cdf_grid(dat = dat, bws = bw, ngrid = 7L)
    expect_false(grid$cdf_on_train)
    expect_identical(grid$nog, 7L)
    expect_equal(grid$gcon[, 1], as.numeric(uocquantile(dat$x1, probs)))
    expect_equal(grid$gcon[, 2], as.numeric(uocquantile(dat$x2, probs)))

    explicit <- .npudistbw_prepare_cdf_grid(
      dat = dat,
      bws = bw,
      gdat = dat[3:5, , drop = FALSE]
    )
    expect_false(explicit$cdf_on_train)
    expect_identical(explicit$nog, 3L)
    expect_equal(unname(explicit$gcon),
                 unname(as.matrix(dat[3:5, , drop = FALSE])))

    full <- .npudistbw_prepare_cdf_grid(
      dat = dat,
      bws = bw,
      do.full.integral = TRUE
    )
    expect_true(full$cdf_on_train)
    expect_identical(full$nog, 0L)
    expect_identical(ncol(full$gcon), 0L)
  })
})

test_that("npudistbw default grid matches explicit marginal-quantile gdat", {
  with_npudist_audit_runtime({
    dat <- data.frame(
      x1 = seq(-1, 1, length.out = 12L),
      x2 = rev(seq(-2, 2, length.out = 12L))
    )
    bw <- npudistbw(dat = dat, bws = c(0.35, 0.45),
                    bandwidth.compute = FALSE)
    grid <- .npudistbw_prepare_cdf_grid(dat = dat, bws = bw, ngrid = 6L)
    gdat <- data.frame(grid$gcon)
    names(gdat) <- names(dat)

    default.path <- npudistbw(dat = dat, bws = bw, eval.only = TRUE,
                              ngrid = 6L)
    explicit.path <- npudistbw(dat = dat, bws = bw, eval.only = TRUE,
                               gdat = gdat)
    full.path <- npudistbw(dat = dat, bws = bw, eval.only = TRUE,
                           do.full.integral = TRUE)

    expect_equal(default.path$fval, explicit.path$fval, tolerance = 1e-14)
    expect_true(is.finite(full.path$fval))
  })
})

test_that("npudist keeps train and evaluation omitted rows separately", {
  with_npudist_audit_runtime({
    train <- data.frame(x = c(0, NA, 1, 2, NA, 3))
    eval.clean <- data.frame(x = c(-0.5, 0.5, 1.5))
    eval.dirty <- data.frame(x = c(-0.5, NA, 1.5, NA))
    bw <- npudistbw(dat = data.frame(x = c(0, 1, 2, 3)),
                    bws = 0.5, bandwidth.compute = FALSE)

    clean.fit <- npudist(bws = bw, tdat = train, edat = eval.clean)
    dirty.fit <- npudist(bws = bw, tdat = train, edat = eval.dirty)

    expect_omit(clean.fit$train.rows.omit, c(2L, 5L))
    expect_identical(clean.fit$ntrain.omit, 2L)
    expect_identical(clean.fit$eval.rows.omit, NA)
    expect_identical(clean.fit$neval.omit, 0)
    expect_identical(clean.fit$rows.omit, NA)
    expect_identical(clean.fit$nobs.omit, 0)

    expect_omit(dirty.fit$train.rows.omit, c(2L, 5L))
    expect_omit(dirty.fit$eval.rows.omit, c(2L, 4L))
    expect_omit(dirty.fit$rows.omit, c(2L, 4L))
    expect_identical(dirty.fit$ntrain.omit, 2L)
    expect_identical(dirty.fit$neval.omit, 2L)
    expect_identical(dirty.fit$nobs.omit, 2L)
  })
})

test_that("npudist zero-row inputs fail before low-level shape errors", {
  with_npudist_audit_runtime({
    bw <- npudistbw(dat = data.frame(x = c(0, 1, 2, 3)),
                    bws = 0.5, bandwidth.compute = FALSE)

    expect_error(
      npudist(bws = bw, tdat = data.frame(x = c(NA_real_, NA_real_))),
      "no rows without NAs"
    )
    expect_error(
      npudist(bws = bw, tdat = data.frame(x = c(0, 1, 2, 3)),
              edat = data.frame(x = c(NA_real_, NA_real_))),
      "no rows without NAs"
    )
    expect_error(
      npudistbw(dat = data.frame(x = c(NA_real_, NA_real_)),
                nmulti = 1L, itmax = 1L),
      "no rows without NAs"
    )
  })
})

test_that("npudist fit mirrors bandwidth-side unordered rejection", {
  with_npudist_audit_runtime({
    dat <- data.frame(u = factor(c("a", "a", "b", "b")))
    bw <- dbandwidth(bw = 0, nobs = nrow(dat), xdati = untangle(dat),
                     xnames = names(dat), bandwidth.compute = FALSE,
                     bwmethod = "cv.cdf")

    expect_error(
      npudist(bws = bw, tdat = dat),
      "does not support unordered data types"
    )
  })
})

test_that("npudistbw validates method metadata without rejecting normal-reference", {
  with_npudist_audit_runtime({
    x <- data.frame(x = seq(0, 1, length.out = 8L))
    bw <- npudistbw(dat = x, bws = 0.3, bandwidth.compute = FALSE)
    bad <- bw
    bad$method <- "bogus"

    expect_error(
      npudistbw(dat = x, bws = bad, nmulti = 1L, itmax = 1L),
      "does not support bwmethod"
    )

    nr <- npudistbw(dat = x, bwmethod = "normal-reference")
    expect_identical(nr$method, "normal-reference")
    expect_true(is.finite(nr$bw[1]))
  })
})
