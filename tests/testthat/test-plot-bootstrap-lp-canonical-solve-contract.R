test_that("fixed-LP bootstrap batch projection uses one canonical response solve", {
  project <- getFromNamespace(".np_inid_lp_batch_project", "npRmpi")
  set.seed(81701)

  for (p in c(1L, 2L, 3L, 9L)) {
    x <- matrix(rnorm((p + 4L) * p), ncol = p)
    gram <- crossprod(x) + diag(0.5, p)
    projection <- rnorm(p)
    response <- rnorm(p)
    packed <- unlist(lapply(seq_len(p), function(a) gram[a, a:p]))

    scalar <- project(
      Mvals = matrix(packed, nrow = 1L),
      Zvals = matrix(response, nrow = 1L),
      rhs = projection,
      represented.mass = p + 4L,
      diagnostics = TRUE
    )
    expect_identical(scalar$status, 0L)
    expect_identical(scalar$ridge_steps, 0L)
    expect_identical(scalar$ridge_total, 0)
    expect_equal(
      scalar$values[1L, 1L],
      drop(projection %*% solve(gram, response)),
      tolerance = 1e-12,
      info = sprintf("single response width %d", p)
    )

    responses <- matrix(rnorm(p * 3L), nrow = p)
    moments <- lapply(seq_len(p), function(a) {
      matrix(responses[a, ], nrow = 1L)
    })
    multi <- project(
      Mvals = matrix(packed, nrow = 1L),
      Zvals = moments,
      rhs = projection,
      represented.mass = p + 4L,
      diagnostics = TRUE
    )
    expect_identical(multi$status, 0L)
    expect_identical(multi$ridge_steps, 0L)
    expect_identical(multi$ridge_total, 0)
    expect_equal(
      drop(multi$values[1L, ]),
      drop(projection %*% solve(gram, responses)),
      tolerance = 1e-12,
      info = sprintf("multiple response width %d", p)
    )
  }
})

test_that("fixed-LP bootstrap batch projection has typed transactional failures", {
  project <- getFromNamespace(".np_inid_lp_batch_project", "npRmpi")
  registered <- getNativeSymbolInfo("C_np_lp_batch_project", PACKAGE = "npRmpi")
  expect_identical(registered$dll[["name"]], "npRmpi")

  expect_error(
    project(
      Mvals = matrix(c(1, 0), nrow = 1L),
      Zvals = matrix(c(1, 1), nrow = 1L),
      rhs = c(1, 1),
      represented.mass = 2
    ),
    "packed_gram width does not match projection length"
  )
  expect_error(
    project(
      Mvals = matrix(c(1, 0, 1), nrow = 1L),
      Zvals = list(matrix(1, nrow = 1L)),
      rhs = c(1, 1),
      represented.mass = 2
    ),
    "multi-response moments must be a length-p list"
  )
  expect_error(
    project(
      Mvals = matrix(c(1, 0, 1), nrow = 1L),
      Zvals = matrix(c(1, 1), nrow = 1L),
      rhs = c(1, 1),
      represented.mass = 0
    ),
    "represented_mass must contain positive finite values"
  )

  nonfinite <- .Call(
    "C_np_lp_batch_project",
    matrix(c(1, 0, NaN), nrow = 1L),
    matrix(c(1, 1), nrow = 1L),
    c(1, 1),
    2,
    TRUE,
    PACKAGE = "npRmpi"
  )
  expect_identical(nonfinite$status, 2L)
  expect_null(nonfinite$values)
  expect_identical(nonfinite$failed_system, 1L)

  ridged <- .Call(
    "C_np_lp_batch_project",
    matrix(c(1, 0, 0), nrow = 1L),
    matrix(c(1, 1), nrow = 1L),
    c(1, 1),
    2,
    TRUE,
    PACKAGE = "npRmpi"
  )
  expect_identical(ridged$status, 0L)
  expect_identical(ridged$failed_system, 0L)
  expect_identical(ridged$ridge_steps, 1L)
  expect_identical(ridged$ridge_total, 0.5)

  final.failed <- .Call(
    "C_np_lp_batch_project",
    matrix(c(1, 0, 0), nrow = 1L),
    matrix(c(1, 1), nrow = 1L),
    c(1, 2),
    1e308,
    TRUE,
    PACKAGE = "npRmpi"
  )
  expect_identical(final.failed$status, 4L)
  expect_null(final.failed$values)
  expect_identical(final.failed$failed_system, 1L)
  expect_identical(final.failed$ridge_steps, 1L)

  ordinary <- .Call(
    "C_np_lp_batch_project",
    matrix(c(2, 0, 3), nrow = 1L),
    matrix(c(1, 1), nrow = 1L),
    c(1, 1),
    2,
    FALSE,
    PACKAGE = "npRmpi"
  )
  expect_identical(ordinary$status, 0L)
  expect_null(ordinary$ridge_steps)
  expect_null(ordinary$ridge_total)
})

test_that("migrated fixed-LP index and smooth-coefficient consumers match refits", {
  withr::local_options(npRmpi.local.regression.mode = TRUE)
  local_eval <- getFromNamespace(".np_plot_with_local_compiled_eval", "npRmpi")
  set.seed(81702)
  n <- 24L
  counts <- rmultinom(n = 3L, size = n, prob = rep.int(1 / n, n))

  x1 <- runif(n)
  x2 <- runif(n)
  index.x <- data.frame(x1 = x1, x2 = x2)
  index.y <- sin(x1 + x2) + rnorm(n, sd = 0.05)
  index.bw <- npindexbw(
    xdat = index.x,
    ydat = index.y,
    bws = c(1, 1, 0.32),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 2L,
    basis = "tensor"
  )
  index.fast <- getFromNamespace(".np_inid_boot_from_index", "npRmpi")(
    xdat = index.x,
    ydat = index.y,
    bws = index.bw,
    B = ncol(counts),
    counts = counts
  )
  index.refit <- t(vapply(seq_len(ncol(counts)), function(b) {
    ii <- rep.int(seq_len(n), counts[, b])
    local_eval(npindex(
      txdat = index.x[ii, , drop = FALSE],
      tydat = index.y[ii],
      exdat = index.x,
      bws = index.bw,
      gradients = FALSE
    ))$mean
  }, numeric(n)))
  expect_equal(index.fast$t, index.refit, tolerance = 1e-8)

  sx <- data.frame(x = runif(n, -1, 1))
  sz <- data.frame(z = rnorm(n))
  sy <- with(sx, (0.5 + x) * (1 + sz$z) + rnorm(n, sd = 0.05))
  sex <- data.frame(x = seq(-0.8, 0.8, length.out = 4L))
  sez <- data.frame(z = seq(-0.7, 0.7, length.out = 4L))
  scoef.bw <- npscoefbw(
    xdat = sx,
    zdat = sz,
    ydat = sy,
    bws = 0.7,
    bwtype = "fixed",
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 2L,
    basis = "tensor",
    bernstein.basis = FALSE
  )
  scoef.fast <- getFromNamespace(".np_inid_boot_from_scoef", "npRmpi")(
    txdat = sx,
    ydat = sy,
    tzdat = sz,
    exdat = sex,
    ezdat = sez,
    bws = scoef.bw,
    B = ncol(counts),
    counts = counts,
    mode = "exact"
  )
  scoef.refit <- t(vapply(seq_len(ncol(counts)), function(b) {
    ii <- rep.int(seq_len(n), counts[, b])
    as.vector(local_eval(npscoef(
      bws = scoef.bw,
      txdat = sx[ii, , drop = FALSE],
      tzdat = sz[ii, , drop = FALSE],
      tydat = sy[ii],
      exdat = sex,
      ezdat = sez,
      iterate = FALSE,
      se = FALSE
    ))$mean)
  }, numeric(nrow(sex))))
  expect_equal(scoef.fast$t, scoef.refit, tolerance = 1e-8)
})

test_that("migrated fixed-LP partial and conditional consumers match refits", {
  withr::local_options(npRmpi.local.regression.mode = TRUE)
  local_eval <- getFromNamespace(".np_plot_with_local_compiled_eval", "npRmpi")
  set.seed(81703)
  n <- 24L
  counts <- rmultinom(n = 3L, size = n, prob = rep.int(1 / n, n))

  px <- data.frame(x = runif(n, -1, 1))
  pz <- data.frame(z = rnorm(n))
  py <- 0.7 * px$x + sin(pz$z) + rnorm(n, sd = 0.05)
  pl.bw <- npplregbw(
    xdat = px,
    zdat = pz,
    ydat = py,
    bws = matrix(c(0.55, 0.55), nrow = 2L),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 2L,
    basis = "tensor",
    bernstein.basis = FALSE
  )
  pl.fast <- getFromNamespace(".np_inid_boot_from_plreg", "npRmpi")(
    txdat = px,
    ydat = py,
    tzdat = pz,
    exdat = px,
    ezdat = pz,
    bws = pl.bw,
    B = ncol(counts),
    counts = counts
  )
  pl.refit <- t(vapply(seq_len(ncol(counts)), function(b) {
    ii <- rep.int(seq_len(n), counts[, b])
    local_eval(npplreg(
      bws = pl.bw,
      txdat = px[ii, , drop = FALSE],
      tydat = py[ii],
      tzdat = pz[ii, , drop = FALSE],
      exdat = px,
      ezdat = pz
    ))$mean
  }, numeric(n)))
  expect_equal(pl.fast$t, pl.refit, tolerance = 1e-8)

  cx <- data.frame(x = rnorm(n))
  cy <- data.frame(y = 0.5 * cx$x + rnorm(n))
  grid <- expand.grid(y = seq(-1, 1, length.out = 3L), x = c(-0.5, 0.5))
  cex <- data.frame(x = grid$x)
  cey <- data.frame(y = grid$y)
  cond.bw <- npcdensbw(
    xdat = cx,
    ydat = cy,
    bws = c(0.7, 0.7),
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    regtype = "lp",
    degree = 2L,
    basis = "tensor",
    bernstein.basis = FALSE
  )
  cond.fast <- getFromNamespace(".np_inid_boot_from_ksum_conditional", "npRmpi")(
    xdat = cx,
    ydat = cy,
    exdat = cex,
    eydat = cey,
    bws = cond.bw,
    B = ncol(counts),
    cdf = FALSE,
    counts = counts
  )
  cond.refit <- t(vapply(seq_len(ncol(counts)), function(b) {
    ii <- rep.int(seq_len(n), counts[, b])
    local_eval(npcdens(
      txdat = cx[ii, , drop = FALSE],
      tydat = cy[ii, , drop = FALSE],
      exdat = cex,
      eydat = cey,
      bws = cond.bw
    ))$condens
  }, numeric(nrow(cex))))
  expect_equal(cond.fast$t, cond.refit, tolerance = 1e-8)
})
