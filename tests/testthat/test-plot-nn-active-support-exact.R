test_that("nonfixed unconditional exact bootstrap matches duplicate-row oracle", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  xdat <- data.frame(x = c(0, 0, 1, 3, 3, 10))
  exdat <- data.frame(x = c(0, 2, 5, 9))
  counts <- cbind(
    c(1, 0, 1, 2, 0, 2),
    c(0, 2, 1, 0, 1, 2)
  )
  storage.mode(counts) <- "double"

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw.val <- if (identical(bt, "adaptive_nn")) 2 else 1
    bw <- npudensbw(
      dat = xdat,
      bwtype = bt,
      bws = bw.val,
      bandwidth.compute = FALSE
    )

    helper <- np:::.np_inid_boot_from_ksum_unconditional_exact(
      xdat = xdat,
      exdat = exdat,
      bws = bw,
      B = ncol(counts),
      operator = "normal",
      counts = counts
    )

    manual <- vapply(seq_len(ncol(counts)), function(j) {
      idx <- np:::.np_counts_to_indices(counts[, j])
      np:::.np_ksum_unconditional_eval_exact(
        xdat = xdat[idx, , drop = FALSE],
        exdat = exdat,
        bws = bw,
        operator = "normal"
      )
    }, numeric(nrow(exdat)))

    expect_equal(helper$t0, np:::.np_ksum_unconditional_eval_exact(
      xdat = xdat,
      exdat = exdat,
      bws = bw,
      operator = "normal"
    ), info = bt)
    expect_equal(helper$t, t(manual), tolerance = 1e-14, info = bt)
  }
})

test_that("active bootstrap sample preserves exact row-sliced training payload", {
  xdat <- data.frame(
    x = c(0, 1, 3, 7),
    xf = factor(c("a", "b", "b", "c"))
  )
  ydat <- data.frame(y = c(2, 4, 6, 8))
  counts <- c(0, 2, 0, 1)

  out <- np:::.np_active_boot_sample(xdat = xdat, ydat = ydat, counts.col = counts)

  expect_identical(out$xdat, xdat[c(2, 4), , drop = FALSE])
  expect_identical(out$ydat, ydat[c(2, 4), , drop = FALSE])
  expect_identical(out$weights, matrix(c(2, 1), ncol = 1L))
  expect_identical(out$n.total, sum(counts))
})

test_that("fast exact frame bind matches data.frame semantics", {
  xdat <- data.frame("a b" = c(1, 2), xf = factor(c("a", "b")), check.names = FALSE)
  ydat <- data.frame("a b" = c(3, 4), y = c(5, 6), check.names = FALSE)

  expect_identical(
    np:::.np_bind_data_frames_fast(xdat, ydat),
    data.frame(xdat, ydat)
  )
})

test_that("nonfixed unconditional exact helper matches direct kbandwidth precompute", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  xdat <- data.frame(x = c(0, 0, 1, 3, 3, 10))
  exdat <- data.frame(x = c(0, 2, 5, 9))
  counts <- c(1, 1, 1, 2, 1, 2)
  active <- counts > 0

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw.val <- if (identical(bt, "adaptive_nn")) 2 else 1
    bw <- npudensbw(
      dat = xdat,
      bwtype = bt,
      bws = bw.val,
      bandwidth.compute = FALSE
    )

    kb <- np:::.np_make_kbandwidth_unconditional(bws = bw, xdat = xdat)
    weighted.active <- matrix(as.double(counts[active]), ncol = 1L)

    direct.bw <- np:::.np_ksum_unconditional_eval_exact(
      xdat = xdat[active, , drop = FALSE],
      exdat = exdat,
      bws = bw,
      operator = "normal",
      weights = weighted.active,
      n.total = sum(counts)
    )

    direct.kb <- np:::.np_ksum_unconditional_eval_exact(
      xdat = xdat[active, , drop = FALSE],
      exdat = exdat,
      bws = kb,
      operator = "normal",
      weights = weighted.active,
      n.total = sum(counts)
    )

    expect_equal(direct.kb, direct.bw, tolerance = 0, info = bt)
  }
})

test_that("nonfixed conditional exact bootstrap matches duplicate-row oracle", {
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

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw.val <- if (identical(bt, "adaptive_nn")) c(2, 2) else c(1, 1)
    bw <- npcdensbw(
      xdat = xdat,
      ydat = ydat,
      bwtype = bt,
      bws = bw.val,
      bandwidth.compute = FALSE
    )

    helper <- np:::.np_inid_boot_from_ksum_conditional_exact(
      xdat = xdat,
      ydat = ydat,
      exdat = exdat,
      eydat = eydat,
      bws = bw,
      B = ncol(counts),
      cdf = FALSE,
      counts = counts
    )

    kbx <- np:::.np_con_make_kbandwidth_x(bws = bw, xdat = xdat)
    kbxy <- np:::.np_con_make_kbandwidth_xy(bws = bw, xdat = xdat, ydat = ydat)

    manual <- vapply(seq_len(ncol(counts)), function(j) {
      idx <- np:::.np_counts_to_indices(counts[, j])
      np:::.np_ksum_conditional_eval_exact(
        xdat = xdat[idx, , drop = FALSE],
        ydat = ydat[idx, , drop = FALSE],
        exdat = exdat,
        eydat = eydat,
        kbx = kbx,
        kbxy = kbxy,
        cdf = FALSE
      )
    }, numeric(nrow(exdat)))

    expect_equal(helper$t0, np:::.np_ksum_conditional_eval_exact(
      xdat = xdat,
      ydat = ydat,
      exdat = exdat,
      eydat = eydat,
      kbx = kbx,
      kbxy = kbxy,
      cdf = FALSE
    ), info = bt)
    expect_equal(helper$t, t(manual), tolerance = 1e-14, info = bt)
  }
})

test_that("generalized conditional exact state apply matches weighted active-support oracle", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  xdat <- data.frame(x = c(0, 0, 1, 3, 3, 10))
  ydat <- data.frame(y = c(0, 0, 1, 2, 2, 4))
  exdat <- data.frame(x = c(0, 2, 5, 9))
  eydat <- data.frame(y = c(0, 1, 2, 4))
  counts <- c(1, 0, 1, 2, 0, 2)
  active <- counts > 0

  bw <- npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bwtype = "generalized_nn",
    bws = c(1, 1),
    bandwidth.compute = FALSE
  )

  den.state <- np:::.np_ksum_exact_state_build(
    bws = np:::.np_con_make_kbandwidth_x(bws = bw, xdat = xdat),
    exdat = exdat,
    operator = "normal"
  )
  num.state <- np:::.np_ksum_exact_state_build(
    bws = np:::.np_con_make_kbandwidth_xy(bws = bw, xdat = xdat, ydat = ydat),
    exdat = np:::.np_bind_data_frames_fast(exdat, eydat),
    operator = c("normal", "normal")
  )

  sample <- np:::.np_active_boot_sample_matrix(
    xmat = data.matrix(xdat),
    ymat = cbind(data.matrix(xdat), data.matrix(ydat)),
    counts.col = counts
  )

  den <- as.numeric(np:::.np_ksum_eval_exact_state(
    state = den.state,
    txdat = sample$xmat,
    weights = sample$weights
  )) / sample$n.total

  num <- as.numeric(np:::.np_ksum_eval_exact_state(
    state = num.state,
    txdat = sample$ymat,
    weights = sample$weights
  )) / sample$n.total

  expect_equal(
    num / pmax(den, .Machine$double.eps),
    np:::.np_ksum_conditional_eval_exact(
      xdat = xdat[active, , drop = FALSE],
      ydat = ydat[active, , drop = FALSE],
      exdat = exdat,
      eydat = eydat,
      kbx = np:::.np_con_make_kbandwidth_x(bws = bw, xdat = xdat),
      kbxy = np:::.np_con_make_kbandwidth_xy(bws = bw, xdat = xdat, ydat = ydat),
      cdf = FALSE,
      weights = counts[active],
      n.total = sum(counts)
    ),
    tolerance = 0
  )
})

test_that("adaptive conditional exact handles tiny-support resamples consistently", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  xdat <- data.frame(x = c(0, 0, 1, 3, 3, 10))
  ydat <- data.frame(y = c(0, 0, 1, 2, 2, 4))
  exdat <- data.frame(x = c(0, 2, 5, 9))
  eydat <- data.frame(y = c(0, 1, 2, 4))
  counts <- matrix(c(0, 2, 1, 0, 3, 0), ncol = 1L)
  storage.mode(counts) <- "double"

  bw <- npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bwtype = "adaptive_nn",
    bws = c(2, 2),
    bandwidth.compute = FALSE
  )

  idx <- np:::.np_counts_to_indices(counts[, 1L])
  explicit <- npcdens(
    txdat = xdat[idx, , drop = FALSE],
    tydat = ydat[idx, , drop = FALSE],
    exdat = exdat,
    eydat = eydat,
    bws = npcdensbw(
      xdat = xdat[idx, , drop = FALSE],
      ydat = ydat[idx, , drop = FALSE],
      bwtype = "adaptive_nn",
      bws = c(2, 2),
      bandwidth.compute = FALSE
    )
  )

  helper <- np:::.np_inid_boot_from_ksum_conditional_exact(
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    eydat = eydat,
    bws = bw,
    B = 1L,
    cdf = FALSE,
    counts = counts
  )

  expect_equal(as.numeric(helper$t[1L, ]), as.numeric(explicit$condens), tolerance = 1e-15)
})
