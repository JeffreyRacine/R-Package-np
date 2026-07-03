npcdist_audit_pkg <- if ("package:npRmpi" %in% search()) "npRmpi" else "np"
npcdist_audit_ns <- asNamespace(npcdist_audit_pkg)

npcdist_audit_call <- function(name, ...) {
  do.call(get(name, envir = npcdist_audit_ns), list(...))
}

npcdist_audit_expect_omit <- function(x, expected) {
  if (!length(expected)) {
    expect_true(is.null(x) || length(x) == 0L || identical(x, NA))
  } else {
    expect_equal(unname(as.integer(x)), expected)
  }
}

npcdist_audit_with_runtime <- function(expr) {
  if (identical(npcdist_audit_pkg, "npRmpi")) {
    if (exists("spawn_mpi_slaves", mode = "function") &&
        exists("close_mpi_slaves", mode = "function")) {
      skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
      on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
    } else {
      npRmpi_init <- get("npRmpi.init", envir = npcdist_audit_ns)
      npRmpi_quit <- get("npRmpi.quit", envir = npcdist_audit_ns)
      ok <- try(npRmpi_init(nslaves = 1L, quiet = TRUE), silent = TRUE)
      if (inherits(ok, "try-error"))
        skip("Could not initialize MPI slaves")
      on.exit(try(npRmpi_quit(force = TRUE), silent = TRUE), add = TRUE)
    }
  }
  force(expr)
}

npcdist_audit_grid <- function(xgrid, ygrid) {
  out <- do.call(rbind, lapply(xgrid, function(xx) {
    data.frame(y = ygrid, x = xx)
  }))
  row.names(out) <- NULL
  out
}

test_that("npcdist formula reentry honors explicit data argument", {
  npcdist_audit_with_runtime({
    set.seed(1101)
    d1 <- data.frame(x = seq(-1, 1, length.out = 44L))
    d1$y <- 0.25 + d1$x + rnorm(nrow(d1), sd = 0.2)
    d2 <- data.frame(x = seq(2, 4, length.out = 31L))
    d2$y <- -1 + 0.5 * d2$x + rnorm(nrow(d2), sd = 0.2)

    bw <- npcdist_audit_call(
      "npcdistbw",
      formula = y ~ x,
      data = d1,
      bws = c(0.35, 0.35),
      bandwidth.compute = FALSE,
      bwtype = "fixed"
    )

    fit <- npcdist_audit_call("npcdist", bws = bw, data = d2)
    direct <- npcdist_audit_call(
      "npcdist",
      bws = bw,
      txdat = d2["x"],
      tydat = d2["y"]
    )

    expect_equal(fit$ntrain, nrow(d2))
    expect_equal(fitted(fit), fitted(direct), tolerance = 1e-12)
  })
})

test_that("npcdist records train and evaluation omit metadata without changing legacy rows.omit", {
  npcdist_audit_with_runtime({
    tx <- data.frame(x = c(0.0, 0.2, NA, 0.5, 0.9, 1.1))
    ty <- data.frame(y = c(1, 2, 3, 4, 5, NA))
    ex <- data.frame(x = c(0.1, NA, 0.7, 0.8))
    ey <- data.frame(y = c(1, 2, NA, 4))

    bw <- npcdist_audit_call(
      "npcdistbw",
      xdat = tx,
      ydat = ty,
      bws = c(0.4, 0.4),
      bandwidth.compute = FALSE,
      bwtype = "fixed"
    )
    fit <- npcdist_audit_call(
      "npcdist",
      bws = bw,
      txdat = tx,
      tydat = ty,
      exdat = ex,
      eydat = ey
    )
    fit.train <- npcdist_audit_call("npcdist", bws = bw, txdat = tx, tydat = ty)

    expect_equal(unname(as.integer(fit$rows.omit)), c(2L, 3L))
    expect_equal(unname(as.integer(fit$train.rows.omit)), c(3L, 6L))
    expect_equal(unname(as.integer(fit$eval.rows.omit)), c(2L, 3L))
    expect_equal(fit$train.nobs.omit, 2L)
    expect_equal(fit$eval.nobs.omit, 2L)
    expect_equal(unname(as.integer(fit.train$rows.omit)), c(3L, 6L))
    expect_equal(unname(as.integer(fit.train$train.rows.omit)), c(3L, 6L))
    npcdist_audit_expect_omit(fit.train$eval.rows.omit, integer(0))
  })
})

test_that("npcdist fit uses the same generalized-NN tree predicate as npcdistbw", {
  npcdist_audit_with_runtime({
    set.seed(2202)
    x <- data.frame(x = seq(-0.8, 0.8, length.out = 18L))
    y <- data.frame(y = sin(2 * x$x))
    bw <- npcdist_audit_call(
      "npcdistbw",
      xdat = x,
      ydat = y,
      bws = c(5, 5),
      bandwidth.compute = FALSE,
      bwmethod = "cv.ls",
      bwtype = "generalized_nn",
      cxkertype = "epanechnikov",
      cykertype = "epanechnikov"
    )

    tree_code <- get(".npcdistbw_tree_code", envir = npcdist_audit_ns)
    plain_tree <- get("npDoTreeOrCategoricalCompress", envir = npcdist_audit_ns)
    do_tree_no <- get("DO_TREE_NO", envir = npcdist_audit_ns)
    ncon <- bw$yncon + bw$xncon
    ncat <- bw$ynuno + bw$ynord + bw$xnuno + bw$xnord

    old_opts <- options(np.tree = TRUE, np.categorical.compress = FALSE)
    on.exit(options(old_opts), add = TRUE)

    expect_identical(tree_code(bw, ncon = ncon, ncat = ncat), do_tree_no)
    expect_false(identical(plain_tree(ncon = ncon, ncat = ncat, bws = bw), do_tree_no))

    fit.tree <- npcdist_audit_call("npcdist", bws = bw, txdat = x, tydat = y)
    options(np.tree = FALSE)
    fit.no.tree <- npcdist_audit_call("npcdist", bws = bw, txdat = x, tydat = y)
    expect_equal(fitted(fit.tree), fitted(fit.no.tree), tolerance = 1e-12)

    bad <- bw
    bad$method <- NULL
    options(np.tree = TRUE)
    expect_error(
      npcdist_audit_call("npcdist", bws = bad, txdat = x, tydat = y),
      "valid bwmethod metadata"
    )
  })
})

test_that("npcdistbw method guards reject incoherent method metadata", {
  npcdist_audit_with_runtime({
    x <- data.frame(x = seq(-0.8, 0.8, length.out = 18L))
    y <- data.frame(y = sin(2 * x$x))
    bw <- npcdist_audit_call(
      "npcdistbw",
      xdat = x,
      ydat = y,
      bws = c(5, 5),
      bandwidth.compute = FALSE,
      bwmethod = "cv.ls",
      bwtype = "generalized_nn",
      cxkertype = "epanechnikov",
      cykertype = "epanechnikov"
    )
    ncon <- bw$yncon + bw$xncon
    ncat <- bw$ynuno + bw$ynord + bw$xnuno + bw$xnord
    tree_code <- get(".npcdistbw_tree_code", envir = npcdist_audit_ns)

    old_opts <- options(np.tree = TRUE, np.categorical.compress = FALSE)
    on.exit(options(old_opts), add = TRUE)

    bad <- bw
    bad$method <- "cv.ml"
    expect_error(tree_code(bad, ncon = ncon, ncat = ncat), "does not support bwmethod")

    missing <- bw
    missing$method <- NULL
    expect_error(tree_code(missing, ncon = ncon, ncat = ncat), "valid bwmethod metadata")
  })
})

test_that("npcdist proper CDF repair remains bounded and monotone on explicit slices", {
  npcdist_audit_with_runtime({
    set.seed(3303)
    n <- 55L
    x <- data.frame(x = runif(n, -1, 1))
    y <- data.frame(y = sin(2 * pi * x$x) + rnorm(n, sd = 0.25))
    ygrid <- seq(-1.2, 1.2, length.out = 9L)
    nd <- npcdist_audit_grid(c(-0.35, 0.45), ygrid)

    bw <- npcdist_audit_call(
      "npcdistbw",
      xdat = x,
      ydat = y,
      bws = c(0.30, 0.24),
      bandwidth.compute = FALSE,
      regtype = "lp",
      degree = 3L
    )
    fit <- npcdist_audit_call(
      "npcdist",
      bws = bw,
      txdat = x,
      tydat = y,
      exdat = nd["x"],
      eydat = nd["y"],
      proper = TRUE
    )
    pred <- predict(fit, newdata = nd, proper = TRUE)

    expect_true(isTRUE(fit$proper.requested))
    expect_true(isTRUE(fit$proper.applied))
    expect_false(is.null(fit$condist.raw))
    expect_true(all(fit$condist >= -1e-10))
    expect_true(all(fit$condist <= 1 + 1e-10))
    for (ii in split(seq_along(fit$condist), rep(seq_len(2L), each = length(ygrid))))
      expect_true(all(diff(fit$condist[ii]) >= -1e-9))
    expect_equal(as.numeric(pred), as.numeric(fitted(fit)), tolerance = 1e-10)
  })
})
