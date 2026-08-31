h7c_explicit_bw <- function(x, y, bws, degree, bwtype = "fixed",
                            ckertype = "gaussian", ckerorder = 2L) {
  beta <- identical(ckertype, "beta")
  npregbw(
    xdat = x, ydat = y, bws = bws,
    bandwidth.compute = FALSE, bwscaling = FALSE,
    bwtype = bwtype, regtype = "lp", degree = degree,
    degree.select = "manual", ckertype = ckertype,
    ckerorder = ckerorder,
    ckerbound = if (beta) "fixed" else "none",
    ckerlb = if (beta) c(0, 0) else NULL,
    ckerub = if (beta) c(1, 1) else NULL
  )
}

h7c_fit <- function(bw, x, y, eval = NULL, se = FALSE,
                    order = c(1L, 1L), warn = FALSE) {
  call <- list(
    bws = bw, txdat = x, tydat = y,
    gradients = TRUE, se = se, gradient.order = order,
    warn.glp.gradient = warn
  )
  if (!is.null(eval))
    call$exdat <- eval
  do.call(npreg, call)
}

h7c_available_hc0_oracle <- function(bw, x, y, eval, order, column) {
  training <- suppressWarnings(h7c_fit(
    bw, x, y, se = FALSE, order = order
  ))
  residual <- y - training$mean
  n.eval <- if (is.null(eval)) nrow(x) else nrow(eval)
  influence <- vapply(seq_len(nrow(x)), function(donor) {
    unit <- numeric(nrow(x))
    unit[donor] <- 1
    suppressWarnings(h7c_fit(
      bw, x, unit, eval = eval, se = FALSE, order = order
    ))$grad[, column]
  }, numeric(n.eval))
  sqrt(rowSums(sweep(influence, 2L, residual, `*`)^2))
}

test_that("partial heterogeneous degrees preserve every available native column", {
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)

  n <- 13L
  x <- data.frame(
    x1 = seq(0.04, 0.96, length.out = n),
    x2 = (((seq_len(n) * 5L) %% n) + 0.5) / n,
    u = factor(rep(c("a", "b", "c"), length.out = n)),
    o = ordered(rep(c("low", "mid", "high"), length.out = n),
                levels = c("low", "mid", "high"))
  )
  y <- sin(2 * x$x1) + 0.3 * x$x2 + 0.35 * (x$u == "b") -
    0.2 * (x$o == "high") + seq_len(n) / 100
  external <- x[c(1L, 4L, 8L, 13L), , drop = FALSE]
  cases <- list(
    list(type = "fixed", bws = c(0.31, 0.34, 0.25, 0.3),
         degree = c(2L, 0L), eval = NULL, available = 1L,
         unavailable = 2L),
    list(type = "generalized_nn", bws = c(7, 7, 0.25, 0.3),
         degree = c(0L, 2L), eval = external, available = 2L,
         unavailable = 1L),
    list(type = "adaptive_nn", bws = c(7, 7, 0.25, 0.3),
         degree = c(2L, 0L), eval = external, available = 1L,
         unavailable = 2L)
  )

  for (case in cases) {
    bw <- h7c_explicit_bw(
      x, y, case$bws, case$degree, bwtype = case$type
    )
    no.se <- suppressWarnings(h7c_fit(
      bw, x, y, eval = case$eval, se = FALSE
    ))
    with.se <- suppressWarnings(h7c_fit(
      bw, x, y, eval = case$eval, se = TRUE
    ))
    oracle <- h7c_available_hc0_oracle(
      bw, x, y, case$eval, c(1L, 1L), case$available
    )

    expect_identical(with.se$mean, no.se$mean)
    expect_identical(with.se$grad, no.se$grad)
    expect_true(all(is.finite(with.se$grad[, case$available])))
    expect_true(all(is.finite(with.se$gerr[, case$available])))
    expect_equal(with.se$gerr[, case$available], oracle,
                 tolerance = 2e-11)
    expect_true(all(is.na(with.se$grad[, case$unavailable])))
    expect_true(all(is.na(with.se$gerr[, case$unavailable])))
    expect_true(all(is.finite(with.se$grad[, 3:4, drop = FALSE])))
    expect_true(all(is.finite(with.se$gerr[, 3:4, drop = FALSE])))
  }
})

test_that("partial availability warning and categorical-only remainder persist", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  n <- 15L
  x <- data.frame(
    x1 = seq(-1, 1, length.out = n),
    x2 = ((seq_len(n) * 4L) %% n) / n,
    u = factor(rep(c("a", "b", "c"), length.out = n))
  )
  y <- sin(x$x1) + 0.2 * x$x2 + 0.3 * (x$u == "b")
  bw <- h7c_explicit_bw(x, y, c(0.42, 0.38, 0.25), c(0L, 0L))

  expect_warning(
    fit <- h7c_fit(
      bw, x, y, se = TRUE, order = c(2L, 1L), warn = TRUE
    ),
    "unavailable"
  )
  expect_true(all(is.na(fit$grad[, 1:2, drop = FALSE])))
  expect_true(all(is.na(fit$gerr[, 1:2, drop = FALSE])))
  expect_true(all(is.finite(fit$grad[, 3L])))
  expect_true(all(is.finite(fit$gerr[, 3L])))
})

test_that("H7C uses the compiled coordinate mask without helper replacement", {
  source <- paste(
    readLines(test_path("..", "..", "R", "np.regression.R"), warn = FALSE),
    collapse = "\n"
  )
  start <- regexpr("if (gradients){", source, fixed = TRUE)[[1L]]
  tail <- substr(source, start, nchar(source))
  finish <- regexpr("if (compute.resid.from.fit)", tail, fixed = TRUE)[[1L]]
  active <- substr(tail, 1L, finish - 1L)

  expect_gt(start, 0L)
  expect_gt(finish, 0L)
  expect_match(source,
               "do.compiled.gradients <- isTRUE(gradients)",
               fixed = TRUE)
  expect_match(active, "glp.gradient.available", fixed = TRUE)
  expect_false(grepl("npreghat", active, fixed = TRUE))
  expect_false(grepl("matrix(NA_real_", active, fixed = TRUE))
})
