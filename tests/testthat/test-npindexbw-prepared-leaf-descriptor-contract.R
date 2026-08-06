test_that("npindex prepared regression leaf refresh matches clean reconstruction", {
  set.seed(20260806)
  n <- 48L
  xmat <- cbind(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  y <- sin(xmat[, 1L] - 0.4 * xmat[, 2L]) + rnorm(n, sd = 0.05)
  bws <- npindexbw(
    xdat = as.data.frame(xmat), ydat = y,
    bws = c(1, -0.4, 0.3), bandwidth.compute = FALSE,
    method = "ichimura", regtype = "lp", degree = 0L,
    bwtype = "fixed", ckertype = "gaussian", ckerbound = "range"
  )
  resolve <- getFromNamespace(".npindex_resolve_spec", "np")
  policy <- getFromNamespace(".npindex_objective_policy", "np")
  prepare <- getFromNamespace(
    ".npindexbw_prepare_lp_regression_leaf_descriptor", "np"
  )
  refresh <- getFromNamespace(".npindexbw_refresh_lp_regression_leaf", "np")

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

test_that("npindex search owners pass one prepared leaf descriptor", {
  source <- paste(deparse(getFromNamespace(".npindexbw_nomad_search", "np")),
                  collapse = "\n")
  expect_match(source, "leaf\\.descriptor <- \\.npindexbw_prepare_lp_regression_leaf_owner", fixed = FALSE)
  expect_match(source, "leaf\\.descriptor = leaf\\.descriptor", fixed = FALSE)
  expect_false(exists(".npindexbw_build_lp_regression_leaf",
                      envir = asNamespace("np"), inherits = FALSE))
})
