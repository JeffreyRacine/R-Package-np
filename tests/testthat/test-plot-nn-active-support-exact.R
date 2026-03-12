test_that("nonfixed unconditional exact bootstrap matches duplicate-row oracle", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  xdat <- data.frame(x = c(0, 0, 1, 3, 3, 10))
  exdat <- data.frame(x = c(0, 2, 5, 9))
  counts <- cbind(
    c(1, 0, 1, 2, 0, 2),
    c(0, 2, 1, 0, 1, 2)
  )
  storage.mode(counts) <- "double"

  bw <- npudensbw(
    dat = xdat,
    bwtype = "generalized_nn",
    bws = 1,
    bandwidth.compute = FALSE
  )

  helper <- npRmpi:::.np_inid_boot_from_ksum_unconditional_exact(
    xdat = xdat,
    exdat = exdat,
    bws = bw,
    B = ncol(counts),
    operator = "normal",
    counts = counts
  )

  manual <- vapply(seq_len(ncol(counts)), function(j) {
    idx <- npRmpi:::.np_counts_to_indices(counts[, j])
    npRmpi:::.np_ksum_unconditional_eval_exact(
      xdat = xdat[idx, , drop = FALSE],
      exdat = exdat,
      bws = bw,
      operator = "normal"
    )
  }, numeric(nrow(exdat)))

  expect_equal(helper$t0, npRmpi:::.np_ksum_unconditional_eval_exact(
    xdat = xdat,
    exdat = exdat,
    bws = bw,
    operator = "normal"
  ))
  expect_equal(helper$t, t(manual), tolerance = 1e-14)
})

test_that("nonfixed conditional exact bootstrap matches duplicate-row oracle", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  xdat <- data.frame(x = c(0, 0, 1, 3, 3, 10))
  ydat <- data.frame(y = c(0, 0, 1, 2, 2, 4))
  exdat <- data.frame(x = c(0, 2, 5, 9))
  eydat <- data.frame(y = c(0, 1, 2, 4))
  counts <- cbind(
    c(1, 0, 1, 2, 0, 2),
    c(0, 2, 1, 0, 1, 2)
  )
  storage.mode(counts) <- "double"

  bw <- npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bwtype = "generalized_nn",
    bws = c(1, 1),
    bandwidth.compute = FALSE
  )

  helper <- npRmpi:::.np_inid_boot_from_ksum_conditional_exact(
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    eydat = eydat,
    bws = bw,
    B = ncol(counts),
    cdf = FALSE,
    counts = counts
  )

  kbx <- npRmpi:::.np_con_make_kbandwidth_x(bws = bw, xdat = xdat)
  kbxy <- npRmpi:::.np_con_make_kbandwidth_xy(bws = bw, xdat = xdat, ydat = ydat)

  manual <- vapply(seq_len(ncol(counts)), function(j) {
    idx <- npRmpi:::.np_counts_to_indices(counts[, j])
    npRmpi:::.np_ksum_conditional_eval_exact(
      xdat = xdat[idx, , drop = FALSE],
      ydat = ydat[idx, , drop = FALSE],
      exdat = exdat,
      eydat = eydat,
      kbx = kbx,
      kbxy = kbxy,
      cdf = FALSE
    )
  }, numeric(nrow(exdat)))

  expect_equal(helper$t0, npRmpi:::.np_ksum_conditional_eval_exact(
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    eydat = eydat,
    kbx = kbx,
    kbxy = kbxy,
    cdf = FALSE
  ))
  expect_equal(helper$t, t(manual), tolerance = 1e-14)
})
