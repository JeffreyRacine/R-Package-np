test_that("operator weight policy owns beta and legacy complete units", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- data.frame(
    x1 = seq(0.05, 0.95, length.out = 12L),
    x2 = seq(0.08, 0.92, length.out = 12L)
  )
  eval <- dat[c(2L, 6L, 10L), , drop = FALSE]
  operators <- c("normal", "integral")

  for (family in c("gaussian", "beta")) {
    args <- list(
      dat = dat,
      bws = c(0.18, 0.23),
      bandwidth.compute = FALSE,
      ckertype = family,
      ckerorder = 4L
    )
    if (identical(family, "beta"))
      args <- c(args, list(
        ckerbound = "fixed", ckerlb = 0, ckerub = 1
      ))
    bw <- do.call(npudensbw, args)
    info <- np:::.np_operator_kernel_weight_scale(
      bws = bw,
      operator = operators,
      nvars = ncol(dat),
      where = "unit policy oracle"
    )

    expect_identical(info$bandwidth.divide, !identical(family, "beta"))
    expect_equal(
      info$scale,
      if (identical(family, "beta")) 1 else
        info$bws[["bw", exact = TRUE]][[1L]],
      tolerance = 0
    )

    explicit <- np:::.np_kernel_weights_direct(
      bws = info$bws,
      txdat = dat,
      exdat = eval,
      bandwidth.divide = !identical(family, "beta"),
      operator = operators
    )
    expected <- t(as.matrix(explicit)) / (nrow(dat) * info$scale)
    got <- np:::.np_ksum_unconditional_operator_fixed(
      xdat = dat, exdat = eval, bws = bw, operator = operators
    )
    expect_identical(got, expected)
  }
})

test_that("audited operator consumers use the shared unit policy", {
  candidates <- unique(c(
    normalizePath(file.path("..", ".."), mustWork = FALSE),
    normalizePath(".", mustWork = FALSE)
  ))
  root <- candidates[vapply(candidates, function(path) {
    file.exists(file.path(path, "R", "np.plot.helpers.R")) &&
      file.exists(file.path(path, "R", "np.cdhat.helpers.R"))
  }, logical(1))][1L]
  if (is.na(root)) root <- ""
  skip_if_not(file.exists(file.path(root, "R", "np.plot.helpers.R")))
  plot.source <- paste(readLines(
    file.path(root, "R", "np.plot.helpers.R"), warn = FALSE
  ), collapse = "\n")
  cdhat.source <- paste(readLines(
    file.path(root, "R", "np.cdhat.helpers.R"), warn = FALSE
  ), collapse = "\n")

  expect_match(
    plot.source,
    "bandwidth.divide = den.info$bandwidth.divide",
    fixed = TRUE
  )
  expect_match(
    plot.source,
    "bandwidth.divide = num.info$bandwidth.divide",
    fixed = TRUE
  )
  expect_gte(
    lengths(regmatches(
      cdhat.source,
      gregexpr("bandwidth.divide = op.info$bandwidth.divide",
               cdhat.source, fixed = TRUE)
    )),
    3L
  )
})
