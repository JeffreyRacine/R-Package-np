test_that("npindex prepared regression leaf refresh matches clean reconstruction", {
  set.seed(20260806)
  n <- 48L
  xmat <- cbind(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  y <- sin(xmat[, 1L] - 0.4 * xmat[, 2L]) + rnorm(n, sd = 0.05)
  bws <- sibandwidth(
    beta = c(1, -0.4), h = 0.3,
    method = "ichimura", regtype = "lp", degree = 0L,
    bwtype = "fixed", ckertype = "gaussian", ckerbound = "range",
    nobs = n, xdati = untangle(as.data.frame(xmat)),
    ydati = untangle(data.frame(y = y)),
    xnames = colnames(xmat), ynames = "y", bandwidth.compute = FALSE
  )
  resolve <- getFromNamespace(".npindex_resolve_spec", "npRmpi")
  policy <- getFromNamespace(".npindex_objective_policy", "npRmpi")
  prepare <- getFromNamespace(
    ".npindexbw_prepare_lp_regression_leaf_descriptor", "npRmpi"
  )
  refresh <- getFromNamespace(".npindexbw_refresh_lp_regression_leaf", "npRmpi")

  degree0 <- policy(bws, resolve(bws), TRUE)$objective.spec
  degree2 <- degree0
  degree2$regtype <- "lp"
  degree2$degree <- 2L
  degree2$regtype.engine <- "lp"
  degree2$degree.engine <- 2L
  degree2$basis.engine <- "glp"
  degree2$bernstein.basis.engine <- TRUE

  index0 <- drop(xmat %*% c(1, -0.4))
  index2 <- drop(xmat %*% c(1, 0.2))
  descriptor <- prepare(index0, y, 0.3, bws, degree0)
  descriptor.bytes <- serialize(descriptor, NULL, version = 3L)
  retained <- refresh(descriptor, index2, 0.22, degree2)
  rebuilt <- refresh(
    prepare(index2, y, 0.22, bws, degree2),
    index2, 0.22, degree2
  )

  expect_identical(retained$xdat, rebuilt$xdat)
  expect_identical(retained$bws, rebuilt$bws)
  expect_identical(serialize(descriptor, NULL, version = 3L), descriptor.bytes)
  expect_identical(retained$bws$xdati, untangle(retained$xdat))
  expect_identical(retained$bws$dati$x, retained$bws$xdati)
  expect_identical(retained$bws$ckerlb, min(index2))
  expect_identical(retained$bws$ckerub, max(index2))
  expect_identical(retained$bws$regtype.engine, "lp")
  expect_identical(retained$bws$degree.engine, 2L)
  expect_true(retained$bws$bernstein.basis.engine)
})

test_that("npindex serial and service owners pass one prepared leaf descriptor", {
  nomad.source <- paste(
    deparse(getFromNamespace(".npindexbw_nomad_search", "npRmpi")),
    collapse = "\n"
  )
  traced.source <- paste(
    deparse(getFromNamespace(".npindexbw_eval_objective_service_traced", "npRmpi")),
    collapse = "\n"
  )
  expect_match(nomad.source,
               "leaf\\.descriptor <- \\.npindexbw_prepare_lp_regression_leaf_owner")
  expect_match(nomad.source, "service\\.ctx\\$leaf\\.descriptor <- leaf\\.descriptor")
  expect_match(traced.source, "leaf\\.descriptor = ctx\\$leaf\\.descriptor")
  expect_false(exists(".npindexbw_build_lp_regression_leaf",
                      envir = asNamespace("npRmpi"), inherits = FALSE))
})
