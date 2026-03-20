shadow_conditional_xweights_lp_nprmpi <- function(bw, xdat, ydat, row_index) {
  make_mat <- function(z) {
    z <- as.matrix(z)
    storage.mode(z) <- "double"
    z
  }

  xmat <- npRmpi:::toMatrix(xdat)
  ymat <- npRmpi:::toMatrix(npRmpi:::toFrame(ydat))
  ycon <- make_mat(ymat[, bw$iycon, drop = FALSE])
  yuno <- make_mat(ymat[, bw$iyuno, drop = FALSE])
  yord <- make_mat(ymat[, bw$iyord, drop = FALSE])
  xcon <- make_mat(xmat[, bw$ixcon, drop = FALSE])
  xuno <- make_mat(xmat[, bw$ixuno, drop = FALSE])
  xord <- make_mat(xmat[, bw$ixord, drop = FALSE])

  bw_vec <- c(
    bw$xbw[bw$ixcon],
    bw$ybw[bw$iycon],
    bw$ybw[bw$iyuno],
    bw$ybw[bw$iyord],
    bw$xbw[bw$ixuno],
    bw$xbw[bw$ixord]
  )

  .Call(
    "C_np_shadow_cv_xweights_conditional",
    yuno,
    yord,
    ycon,
    xuno,
    xord,
    xcon,
    as.double(bw_vec),
    as.integer(npRmpi:::BW_FIXED),
    as.integer(switch(
      bw$cxkertype,
      gaussian = npRmpi:::CKER_GAUSS + bw$cxkerorder / 2 - 1,
      epanechnikov = npRmpi:::CKER_EPAN + bw$cxkerorder / 2 - 1,
      uniform = npRmpi:::CKER_UNI,
      "truncated gaussian" = npRmpi:::CKER_TGAUSS
    )),
    as.integer(switch(
      bw$uxkertype,
      aitchisonaitken = npRmpi:::UKER_AIT,
      liracine = npRmpi:::UKER_LR
    )),
    as.integer(switch(
      bw$oxkertype,
      wangvanryzin = npRmpi:::OKER_WANG,
      liracine = npRmpi:::OKER_LR,
      racineliyan = npRmpi:::OKER_RLY
    )),
    as.integer(FALSE),
    as.integer(npRmpi:::REGTYPE_LP),
    as.integer(bw$degree.engine),
    as.integer(isTRUE(bw$bernstein.basis.engine)),
    as.integer(npRmpi:::npLpBasisCode(bw$basis.engine)),
    as.integer(row_index),
    PACKAGE = "npRmpi"
  )$streamed
}

bounded_gaussian_kernel_nprmpi <- function(x0, X, h, lower, upper) {
  denom <- h * (pnorm((upper - x0) / h) - pnorm((lower - x0) / h))
  dnorm((x0 - X) / h) / denom
}

smooth_cv_loglik_nprmpi <- function(fit_values, cutoff = .Machine$double.xmin) {
  out <- numeric(length(fit_values))
  log_cutoff <- log(cutoff)

  pos <- fit_values > cutoff
  neg <- fit_values < -cutoff
  mid <- !(pos | neg)

  out[pos] <- log(fit_values[pos])
  out[neg] <- -log(abs(fit_values[neg])) + 2 * log_cutoff
  out[mid] <- log_cutoff
  out
}

test_that("npcdensbw cv.ml LP degree-0 bounded objective matches delete-one reconstruction", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(1)
  n <- 80L
  xdat <- data.frame(x = runif(n))
  ydat <- rbeta(n, 1, 1)

  bw <- npRmpi::npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bwmethod = "cv.ml",
    regtype = "lp",
    degree = 0L,
    bws = c(0.15, 0.12),
    bandwidth.compute = FALSE,
    cxkerbound = "range",
    cykerbound = "range"
  )

  manual_rows <- vapply(seq_len(n), function(i) {
    weights <- shadow_conditional_xweights_lp_nprmpi(bw, xdat, ydat, i)
    ky <- bounded_gaussian_kernel_nprmpi(
      ydat[i],
      ydat,
      bw$ybw[1L],
      bw$cykerlb[bw$iycon][1L],
      bw$cykerub[bw$iycon][1L]
    )
    sum(weights * ky)
  }, numeric(1))

  expect_true(all(manual_rows > .Machine$double.xmin))

  manual_objective <- sum(log(manual_rows))
  np_objective <- npRmpi:::.npcdensbw_eval_only(xdat, ydat, bw)$objective

  expect_equal(np_objective, manual_objective, tolerance = 1e-5)
})

test_that("npcdensbw cv.ml LP objective uses smooth penalty for negative delete-one fits", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(1)
  n <- 80L
  xdat <- data.frame(x = runif(n))
  ydat <- rbeta(n, 1, 1)

  bw <- npRmpi::npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bwmethod = "cv.ml",
    regtype = "lp",
    degree = 3L,
    bws = c(0.15, 0.12),
    bandwidth.compute = FALSE,
    cxkerbound = "range",
    cykerbound = "range"
  )

  manual_rows <- vapply(seq_len(n), function(i) {
    weights <- shadow_conditional_xweights_lp_nprmpi(bw, xdat, ydat, i)
    ky <- bounded_gaussian_kernel_nprmpi(
      ydat[i],
      ydat,
      bw$ybw[1L],
      bw$cykerlb[bw$iycon][1L],
      bw$cykerub[bw$iycon][1L]
    )
    sum(weights * ky)
  }, numeric(1))

  expect_gt(sum(manual_rows < 0), 0L)

  smooth_objective <- sum(smooth_cv_loglik_nprmpi(manual_rows))
  constant_terms <- rep.int(log(.Machine$double.xmin), length(manual_rows))
  pos <- manual_rows > .Machine$double.xmin
  constant_terms[pos] <- log(manual_rows[pos])
  constant_objective <- sum(constant_terms)
  np_objective <- npRmpi:::.npcdensbw_eval_only(xdat, ydat, bw)$objective

  expect_lt(abs(np_objective - smooth_objective), 20)
  expect_gt(abs(np_objective - constant_objective), 1000)
})
