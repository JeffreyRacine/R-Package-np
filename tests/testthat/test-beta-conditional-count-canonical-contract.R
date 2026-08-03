locate_beta_conditional_count_sources <- function() {
  candidates <- unique(c(
    normalizePath(file.path("..", ".."), mustWork = FALSE),
    normalizePath(".", mustWork = FALSE)
  ))
  for (root in candidates) {
    if (file.exists(file.path(root, "src", "np.c")) &&
        file.exists(file.path(root, "src", "jksum.c")) &&
        file.exists(file.path(root, "src", "beta_conditional.c")))
      return(root)
  }
  NULL
}

test_that("conditional count bootstrap enters the canonical scaled-row owner", {
  root <- locate_beta_conditional_count_sources()
  skip_if(is.null(root), "package sources unavailable")
  ingress <- paste(readLines(file.path(root, "src", "np.c"), warn = FALSE),
                   collapse = "\n")
  engine <- paste(readLines(file.path(root, "src", "jksum.c"), warn = FALSE),
                  collapse = "\n")
  sidecar <- paste(
    readLines(file.path(root, "src", "beta_conditional.c"), warn = FALSE),
    collapse = "\n"
  )

  expect_match(ingress, "count_status = np_conditional_count_levels(", fixed = TRUE)
  expect_match(ingress, "SEXP C_np_conditional_count_levels(", fixed = TRUE)
  expect_false(grepl("C_np_beta_conditional_bootstrap", ingress, fixed = TRUE))
  expect_false(grepl(
    "conditional_status = np_beta_conditional_lc_counts(",
    ingress, fixed = TRUE
  ))
  expect_match(
    engine,
    "int np_conditional_count_levels(const NPConditionalCountPlan * const plan,",
    fixed = TRUE
  )
  expect_gte(
    lengths(regmatches(engine, gregexpr("F77_CALL(dgemv)", engine, fixed = TRUE))),
    2L
  )
  expect_match(
    engine,
    "np_beta_scaled_row_context_fill(\n        &x_context",
    fixed = TRUE
  )
  expect_match(
    engine,
    "np_beta_scaled_row_context_fill(\n        &y_context",
    fixed = TRUE
  )

  callers <- lengths(regmatches(
    paste(ingress, engine, sep = "\n"),
    gregexpr("np_beta_conditional_lc_counts(",
             paste(ingress, engine, sep = "\n"), fixed = TRUE)
  ))
  expect_identical(callers, 0L)
  expect_match(sidecar, "np_beta_conditional_lc_counts(", fixed = TRUE)
})

test_that("canonical beta count levels match direct weighted conditional rows", {
  set.seed(20260801)
  xdat <- data.frame(x = sort(runif(18L, 0.02, 0.98)))
  ydat <- data.frame(y = sort(runif(18L, 0.02, 0.98)))
  exdat <- data.frame(x = c(0.04, 0.23, 0.52, 0.81, 0.96))
  eydat <- data.frame(y = c(0.05, 0.20, 0.56, 0.78, 0.95))
  counts <- cbind(
    rep.int(1, nrow(xdat)),
    c(2, 0, rep.int(1, nrow(xdat) - 2L)),
    c(0, 2, rep.int(1, nrow(xdat) - 2L))
  )

  for (cdf in c(FALSE, TRUE)) {
    constructor <- if (cdf) npcdistbw else npcdensbw
    for (placement in c("x", "y", "both")) {
      for (order in c(2L, 4L, 6L, 8L)) {
        args <- list(
          xdat = xdat, ydat = ydat,
          bws = c(0.18, 0.16), bandwidth.compute = FALSE,
          cxkertype = if (placement %in% c("x", "both")) "beta" else "gaussian",
          cxkerorder = order,
          cykertype = if (placement %in% c("y", "both")) "beta" else "gaussian",
          cykerorder = order
        )
        if (placement %in% c("x", "both"))
          args <- c(args, list(
            cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1
          ))
        if (placement %in% c("y", "both"))
          args <- c(args, list(
            cykerbound = "fixed", cykerlb = 0, cykerub = 1
          ))
        bw <- do.call(constructor, args)
        got <- np:::.np_conditional_count_levels(
          xdat = xdat, ydat = ydat, exdat = exdat, eydat = eydat,
          bws = bw, cdf = cdf, counts = counts
        )
        ops <- np:::.np_conditional_side_operators(
          xdat = xdat, ydat = ydat, exdat = exdat, eydat = eydat,
          bws = bw, cdf = cdf
        )
        expected <- t(ops$num %*% counts) /
          np:::NZD(t(ops$den %*% counts))
        expect_equal(got, expected, tolerance = 5e-12)
      }
    }
  }
})

test_that("canonical conditional count levels include categorical sides", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  n <- 18L
  xdat <- data.frame(
    x = seq(0.03, 0.97, length.out = n),
    u = factor(rep(c("a", "b", "c"), length.out = n)),
    o = ordered(rep(c("lo", "mid", "hi"), each = 6L),
                levels = c("lo", "mid", "hi"))
  )
  ydat <- data.frame(
    y = seq(0.04, 0.96, length.out = n),
    v = factor(rep(c("r", "s", "t"), each = 6L))
  )
  exdat <- xdat[c(1L, 5L, 9L, 13L, 17L), , drop = FALSE]
  eydat <- ydat[c(2L, 6L, 10L, 14L, 18L), , drop = FALSE]
  counts <- cbind(
    rep.int(1, n),
    c(2, 0, rep.int(1, n - 2L)),
    c(0, 2, rep.int(1, n - 2L))
  )

  for (cdf in c(FALSE, TRUE)) {
    constructor <- if (cdf) npcdistbw else npcdensbw
    estimator <- if (cdf) npcdist else npcdens
    result.name <- if (cdf) "condist" else "condens"
    bw <- do.call(constructor, list(
      xdat = xdat,
      ydat = ydat,
      bws = c(0.21, 0.22, 0.31, 0.18, 0.27),
      bandwidth.compute = FALSE,
      cxkertype = "beta",
      cxkerorder = 6L,
      cxkerbound = "fixed",
      cxkerlb = 0,
      cxkerub = 1,
      cykertype = "gaussian",
      cykerorder = 4L
    ))
    expected <- t(vapply(seq_len(ncol(counts)), function(j) {
      idx <- np:::.np_counts_to_indices(counts[, j])
      fit <- do.call(estimator, list(
        bws = bw,
        txdat = xdat[idx, , drop = FALSE],
        tydat = ydat[idx, , drop = FALSE],
        exdat = exdat,
        eydat = eydat,
        gradients = FALSE,
        proper = FALSE
      ))
      as.numeric(fit[[result.name]])
    }, numeric(nrow(exdat))))

    observed <- lapply(c(FALSE, TRUE), function(compress) {
      options(np.categorical.compress = compress)
      np:::.np_conditional_count_levels(
        xdat = xdat, ydat = ydat, exdat = exdat, eydat = eydat,
        bws = bw, cdf = cdf, counts = counts
      )
    })
    expect_equal(observed[[1L]], expected, tolerance = 5e-12,
                 info = if (cdf) "cdf dense" else "pdf dense")
    expect_equal(observed[[2L]], expected, tolerance = 5e-12,
                 info = if (cdf) "cdf compressed" else "pdf compressed")
    expect_equal(observed[[2L]], observed[[1L]], tolerance = 5e-14,
                 info = if (cdf) "cdf compression parity" else
                   "pdf compression parity")
  }
})

test_that("a beta response supports an all-categorical explanatory side", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  xdat <- data.frame(
    u = factor(rep(c("a", "b", "c"), 6L)),
    o = ordered(rep(c("lo", "mid", "hi"), each = 6L),
                levels = c("lo", "mid", "hi"))
  )
  ydat <- data.frame(y = seq(0.03, 0.97, length.out = nrow(xdat)))
  exdat <- xdat[c(2L, 7L, 12L, 17L), , drop = FALSE]
  eydat <- ydat[c(1L, 6L, 11L, 18L), , drop = FALSE]
  counts <- cbind(
    rep.int(1, nrow(xdat)),
    c(2, 0, rep.int(1, nrow(xdat) - 2L))
  )
  bw <- npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.24, 0.33, 0.17),
    bandwidth.compute = FALSE,
    cykertype = "beta",
    cykerorder = 8L,
    cykerbound = "fixed",
    cykerlb = 0,
    cykerub = 1
  )
  got <- np:::.np_conditional_count_levels(
    xdat = xdat, ydat = ydat, exdat = exdat, eydat = eydat,
    bws = bw, cdf = FALSE, counts = counts
  )
  expected <- t(vapply(seq_len(ncol(counts)), function(j) {
    idx <- np:::.np_counts_to_indices(counts[, j])
    as.numeric(npcdens(
      bws = bw,
      txdat = xdat[idx, , drop = FALSE],
      tydat = ydat[idx, , drop = FALSE],
      exdat = exdat,
      eydat = eydat,
      gradients = FALSE,
      proper = FALSE
    )$condens)
  }, numeric(nrow(exdat))))

  expect_equal(got, expected, tolerance = 5e-12)
})

test_that("a beta explanatory side supports an all-categorical response", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  xdat <- data.frame(x = seq(0.03, 0.97, length.out = 18L))
  ydat <- data.frame(
    u = factor(rep(c("a", "b", "c"), 6L)),
    o = ordered(rep(c("lo", "mid", "hi"), each = 6L),
                levels = c("lo", "mid", "hi"))
  )
  exdat <- xdat[c(2L, 7L, 12L, 17L), , drop = FALSE]
  eydat <- ydat[c(1L, 6L, 11L, 18L), , drop = FALSE]
  counts <- cbind(
    rep.int(1, nrow(xdat)),
    c(2, 0, rep.int(1, nrow(xdat) - 2L))
  )

  for (cdf in c(FALSE, TRUE)) {
    constructor <- if (cdf) npcdistbw else npcdensbw
    estimator <- if (cdf) npcdist else npcdens
    result.name <- if (cdf) "condist" else "condens"
    bw <- do.call(constructor, list(
      xdat = xdat,
      ydat = ydat,
      bws = c(0.26, 0.31, 0.18),
      bandwidth.compute = FALSE,
      cxkertype = "beta",
      cxkerorder = 6L,
      cxkerbound = "fixed",
      cxkerlb = 0,
      cxkerub = 1
    ))
    expected <- t(vapply(seq_len(ncol(counts)), function(j) {
      idx <- np:::.np_counts_to_indices(counts[, j])
      fit <- do.call(estimator, list(
        bws = bw,
        txdat = xdat[idx, , drop = FALSE],
        tydat = ydat[idx, , drop = FALSE],
        exdat = exdat,
        eydat = eydat,
        gradients = FALSE,
        proper = FALSE
      ))
      as.numeric(fit[[result.name]])
    }, numeric(nrow(exdat))))

    for (compress in c(FALSE, TRUE)) {
      options(np.categorical.compress = compress)
      got <- np:::.np_conditional_count_levels(
        xdat = xdat, ydat = ydat, exdat = exdat, eydat = eydat,
        bws = bw, cdf = cdf, counts = counts
      )
      expect_equal(got, expected, tolerance = 5e-12,
                   info = paste(if (cdf) "cdf" else "pdf",
                                "compression", compress))
    }
  }
})

test_that("nonfixed mixed conditional counts expand before NN bandwidths", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  n <- 18L
  xdat <- data.frame(
    x = seq(0.03, 0.97, length.out = n),
    u = factor(rep(c("a", "b", "c"), 6L))
  )
  ydat <- data.frame(
    y = seq(0.04, 0.96, length.out = n),
    v = factor(rep(c("r", "s", "t"), each = 6L))
  )
  exdat <- xdat[c(1L, 5L, 9L, 13L, 17L), , drop = FALSE]
  eydat <- ydat[c(2L, 6L, 10L, 14L, 18L), , drop = FALSE]
  counts <- cbind(
    rep.int(1, n),
    c(2, 0, rep.int(1, n - 2L)),
    c(0, 2, rep.int(1, n - 2L))
  )

  for (cdf in c(FALSE, TRUE)) {
    constructor <- if (cdf) npcdistbw else npcdensbw
    estimator <- if (cdf) npcdist else npcdens
    result.name <- if (cdf) "condist" else "condens"
    for (bwtype in c("generalized_nn", "adaptive_nn")) {
      bw <- do.call(constructor, list(
        xdat = xdat,
        ydat = ydat,
        bws = c(4, 0.22, 4, 0.18),
        bandwidth.compute = FALSE,
        bwtype = bwtype,
        cxkertype = "beta",
        cxkerorder = 6L,
        cxkerbound = "fixed",
        cxkerlb = 0,
        cxkerub = 1,
        cykertype = "gaussian",
        cykerorder = 4L
      ))
      expected <- t(vapply(seq_len(ncol(counts)), function(j) {
        idx <- np:::.np_counts_to_indices(counts[, j])
        fit <- do.call(estimator, list(
          bws = bw,
          txdat = xdat[idx, , drop = FALSE],
          tydat = ydat[idx, , drop = FALSE],
          exdat = exdat,
          eydat = eydat,
          gradients = FALSE,
          proper = FALSE
        ))
        as.numeric(fit[[result.name]])
      }, numeric(nrow(exdat))))

      for (compress in c(FALSE, TRUE)) {
        options(np.categorical.compress = compress)
        got <- np:::.np_conditional_count_levels(
          xdat = xdat, ydat = ydat, exdat = exdat, eydat = eydat,
          bws = bw, cdf = cdf, counts = counts
        )
        expect_equal(got, expected, tolerance = 5e-12,
                     info = paste(bwtype, if (cdf) "cdf" else "pdf",
                                  "compression", compress))
      }
    }
  }
})

test_that("nonfixed beta placements use literal resample NN radii", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  n <- 18L
  xdat <- data.frame(x = seq(0.03, 0.97, length.out = n))
  ydat <- data.frame(y = seq(0.04, 0.96, length.out = n))
  exdat <- xdat[c(1L, 5L, 9L, 13L, 17L), , drop = FALSE]
  eydat <- ydat[c(2L, 6L, 10L, 14L, 18L), , drop = FALSE]
  counts <- cbind(
    rep.int(1, n),
    c(2, 0, rep.int(1, n - 2L)),
    c(0, 2, rep.int(1, n - 2L))
  )

  for (cdf in c(FALSE, TRUE)) {
    constructor <- if (cdf) npcdistbw else npcdensbw
    estimator <- if (cdf) npcdist else npcdens
    result.name <- if (cdf) "condist" else "condens"
    for (bwtype in c("generalized_nn", "adaptive_nn")) {
      for (placement in c("x", "y", "both")) {
        for (order in c(2L, 8L)) {
          args <- list(
            xdat = xdat,
            ydat = ydat,
            bws = c(4, 4),
            bandwidth.compute = FALSE,
            bwtype = bwtype,
            cxkertype = if (placement %in% c("x", "both"))
              "beta" else "gaussian",
            cxkerorder = order,
            cykertype = if (placement %in% c("y", "both"))
              "beta" else "gaussian",
            cykerorder = order
          )
          if (placement %in% c("x", "both"))
            args <- c(args, list(
              cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1
            ))
          if (placement %in% c("y", "both"))
            args <- c(args, list(
              cykerbound = "fixed", cykerlb = 0, cykerub = 1
            ))
          bw <- do.call(constructor, args)
          expected <- t(vapply(seq_len(ncol(counts)), function(j) {
            idx <- np:::.np_counts_to_indices(counts[, j])
            fit <- do.call(estimator, list(
              bws = bw,
              txdat = xdat[idx, , drop = FALSE],
              tydat = ydat[idx, , drop = FALSE],
              exdat = exdat,
              eydat = eydat,
              gradients = FALSE,
              proper = FALSE
            ))
            as.numeric(fit[[result.name]])
          }, numeric(nrow(exdat))))
          got <- np:::.np_conditional_count_levels(
            xdat = xdat, ydat = ydat, exdat = exdat, eydat = eydat,
            bws = bw, cdf = cdf, counts = counts
          )
          expect_equal(
            got, expected, tolerance = 5e-12,
            info = paste(bwtype, placement, order,
                         if (cdf) "cdf" else "pdf")
          )
        }
      }
    }
  }

  expect_error(
    np:::.np_conditional_count_levels(
      xdat = xdat, ydat = ydat, exdat = exdat, eydat = eydat,
      bws = bw, cdf = TRUE, counts = counts,
      nn.sample.expanded = TRUE
    ),
    "expanded sample with unit counts"
  )
})
