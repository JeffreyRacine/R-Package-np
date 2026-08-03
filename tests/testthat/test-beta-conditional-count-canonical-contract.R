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

.beta_conditional_count_test_env <- new.env(parent = emptyenv())

.ensure_beta_conditional_count_pool <- function() {
  if (!isTRUE(.beta_conditional_count_test_env$started)) {
    skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
    .beta_conditional_count_test_env$started <- TRUE
    withr::defer({
      if (isTRUE(.beta_conditional_count_test_env$started)) {
        try(close_mpi_slaves(force = TRUE), silent = TRUE)
        .beta_conditional_count_test_env$started <- FALSE
      }
    }, envir = testthat::teardown_env())
  }
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
  .ensure_beta_conditional_count_pool()
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
        got <- npRmpi:::.np_beta_conditional_bootstrap_levels(
          xdat = xdat, ydat = ydat, exdat = exdat, eydat = eydat,
          bws = bw, cdf = cdf, counts = counts
        )
        ops <- npRmpi:::.np_conditional_side_operators(
          xdat = xdat, ydat = ydat, exdat = exdat, eydat = eydat,
          bws = bw, cdf = cdf
        )
        expected <- t(ops$num %*% counts) /
          npRmpi:::NZD(t(ops$den %*% counts))
        expect_equal(got, expected, tolerance = 5e-12)
      }
    }
  }
})
