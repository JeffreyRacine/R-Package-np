.npcopula_audit_data <- function(n = 42L) {
  set.seed(20260703)
  data.frame(
    x = seq(-0.2, 1.2, length.out = n),
    y = seq(1.2, -0.2, length.out = n) + rnorm(n, sd = 0.04)
  )
}

.npcopula_audit_density_bw <- function(d) {
  npudensbw(
    dat = d,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE,
    ckerbound = "fixed",
    ckerlb = c(-0.2, -0.2),
    ckerub = c(1.2, 1.2)
  )
}

.npcopula_manual_density <- function(bws, d) {
  marginal.bw <- getFromNamespace(".npcopula_marginal_bw", "npRmpi")
  marginal.data <- getFromNamespace(".npcopula_marginal_data", "npRmpi")
  nzd <- getFromNamespace("NZD", "npRmpi")

  joint <- fitted(npudens(bws = bws, tdat = d))
  denom <- rep.int(1.0, nrow(d))
  for (j in seq_along(bws$xnames)) {
    mbw <- marginal.bw(bws, d, j, target = "density")
    tdat <- marginal.data(bws, d, j)
    denom <- denom * fitted(npudens(bws = mbw, tdat = tdat))
  }
  joint / nzd(denom)
}

test_that("npcopula marginal helper maps bounded continuous slots by position", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  d <- .npcopula_audit_data()
  bw <- .npcopula_audit_density_bw(d)
  marginal.bw <- getFromNamespace(".npcopula_marginal_bw", "npRmpi")

  m1 <- marginal.bw(bw, d, 1L, target = "density")
  m2 <- marginal.bw(bw, d, 2L, target = "density")

  expect_equal(m1$ckerlb, bw$ckerlb[1L])
  expect_equal(m1$ckerub, bw$ckerub[1L])
  expect_equal(m2$ckerlb, bw$ckerlb[2L])
  expect_equal(m2$ckerub, bw$ckerub[2L])

  k2 <- marginal.bw(bw, d, 2L, target = "density", kbandwidth = TRUE)
  expect_equal(k2$ckerlb, bw$ckerlb[2L])
  expect_equal(k2$ckerub, bw$ckerub[2L])
})

test_that("npcopula marginal helper keeps continuous bounds off ordered margins", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  d <- data.frame(
    z = ordered(rep(c("low", "mid", "high"), length.out = 42L)),
    x = seq(-0.2, 1.2, length.out = 42L)
  )
  bw <- list(
    xnames = names(d),
    icon = c(FALSE, TRUE),
    bw = c(0.4, 0.35),
    type = "fixed",
    ckerorder = 2L,
    ckertype = "gaussian",
    okertype = "liracine",
    ckerbound = "fixed",
    ckerlb = -0.2,
    ckerub = 1.2
  )
  marginal.args <- getFromNamespace(".npcopula_marginal_bw_args", "npRmpi")

  ordered.margin <- marginal.args(bw, d, 1L, target = "distribution")
  continuous.margin <- marginal.args(bw, d, 2L, target = "distribution")

  expect_null(ordered.margin$ckerbound)
  expect_null(ordered.margin$ckerlb)
  expect_null(ordered.margin$ckerub)
  expect_equal(continuous.margin$ckerlb, bw$ckerlb[1L])
  expect_equal(continuous.margin$ckerub, bw$ckerub[1L])
})

test_that("npcopula bounded density sample path uses canonical marginal denominators", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  d <- .npcopula_audit_data()
  bw <- .npcopula_audit_density_bw(d)

  fit <- npcopula(data = d, bws = bw)
  manual <- .npcopula_manual_density(bw, d)

  expect_equal(fit$copula, manual, tolerance = 1e-10)
})

test_that("npcopula handles non-syntactic names in fixed sample routes", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  d <- .npcopula_audit_data()
  names(d) <- c("x one", "y two")

  dbw <- npudistbw(dat = d, bws = c(0.35, 0.35), bandwidth.compute = FALSE)
  dfit <- npcopula(data = d, bws = dbw)
  expect_s3_class(dfit, "npcopula")
  expect_equal(length(dfit$copula), nrow(d))

  fbw <- npudensbw(dat = d, bws = c(0.35, 0.35), bandwidth.compute = FALSE)
  ffit <- npcopula(data = d, bws = fbw)
  expect_s3_class(ffit, "npcopula")
  expect_equal(length(ffit$copula), nrow(d))
})

test_that("predict.npcopula keeps u explicit and rejects ambiguous newdata names", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  d <- .npcopula_audit_data()
  bw <- npudistbw(dat = d, bws = c(0.35, 0.35), bandwidth.compute = FALSE)
  fit <- npcopula(
    data = d,
    bws = bw,
    u = data.frame(x = c(0.25, 0.75), y = c(0.25, 0.75)),
    n.quasi.inv = 40
  )

  u.named <- data.frame(x = c(0.2, 0.4), y = c(0.3, 0.6))
  u.alias <- data.frame(u1 = c(0.2, 0.4), u2 = c(0.3, 0.6))

  expect_equal(
    predict(fit, u = u.named, n.quasi.inv = 40),
    predict(fit, newdata = u.alias, n.quasi.inv = 40),
    tolerance = 1e-12
  )
  expect_error(
    predict(fit, newdata = u.named, n.quasi.inv = 40),
    "newdata with original variable names is ambiguous",
    fixed = TRUE
  )
})

test_that("npRmpi npcopula bootstrap chunk helpers preserve count order and shape", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  d <- .npcopula_audit_data(12L)
  bw <- npudistbw(dat = d, bws = c(0.35, 0.35), bandwidth.compute = FALSE)
  counts <- matrix(
    c(2L, 0L, 1L, 0L,
      0L, 2L, 1L, 1L,
      rep(1L, 40L)),
    nrow = nrow(d),
    ncol = 4L
  )
  count.tasks <- getFromNamespace(".npcopula_bootstrap_count_tasks", "npRmpi")
  chunk.values <- getFromNamespace(".npcopula_bootstrap_chunk_values", "npRmpi")
  tasks <- count.tasks(counts, chunk.size = 2L)

  expect_equal(vapply(tasks, `[[`, integer(1L), "bsz"), c(2L, 2L))
  expect_equal(do.call(cbind, lapply(tasks, `[[`, "counts")), counts)

  out <- chunk.values(
    task = tasks[[1L]],
    data = d,
    xgrid = d[1:3, , drop = FALSE],
    bws = bw,
    density = FALSE,
    method = "inid",
    blocklen = NULL,
    marginal.bws = NULL,
    nout = 3L
  )

  expect_equal(dim(out), c(2L, 3L))
  expect_true(all(is.finite(out)))
})
