.conditional_exact_lp_test_env <- new.env(parent = emptyenv())

.ensure_conditional_exact_lp_pool <- function() {
  if (!isTRUE(.conditional_exact_lp_test_env$started)) {
    skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
    .conditional_exact_lp_test_env$started <- TRUE
    withr::defer({
      if (isTRUE(.conditional_exact_lp_test_env$started)) {
        try(close_mpi_slaves(force = TRUE), silent = TRUE)
        .conditional_exact_lp_test_env$started <- FALSE
      }
    }, envir = testthat::teardown_env())
  }
}

test_that("nonfixed conditional exact bootstrap selects LP before kernel family", {
  .ensure_conditional_exact_lp_pool()
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  xdat <- data.frame(x = seq(0.05, 0.95, length.out = 14L))
  ydat <- data.frame(y = pmin(
    0.98,
    pmax(0.02, 0.18 + 0.63 * xdat$x + 0.04 * sin(seq_len(nrow(xdat))))
  ))
  exdat <- data.frame(x = c(0.08, 0.27, 0.51, 0.74, 0.93))
  eydat <- data.frame(y = c(0.10, 0.29, 0.53, 0.71, 0.91))
  counts <- cbind(
    c(2, 0, 1, 1, 0, 2, 1, 1, 0, 1, 1, 1, 1, 2),
    c(0, 2, 1, 0, 2, 1, 1, 0, 2, 1, 1, 1, 1, 1)
  )
  storage.mode(counts) <- "double"

  for (cdf in c(FALSE, TRUE)) {
    constructor <- if (cdf) npcdistbw else npcdensbw
    estimator <- if (cdf) npcdist else npcdens
    result.name <- if (cdf) "condist" else "condens"
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      for (bernstein in c(FALSE, TRUE)) {
        bw <- do.call(constructor, list(
          xdat = xdat,
          ydat = ydat,
          bwtype = bwtype,
          bws = if (identical(bwtype, "fixed")) c(0.2, 0.2) else c(4, 4),
          bandwidth.compute = FALSE,
          regtype = "lp",
          degree = 2L,
          bernstein.basis = bernstein,
          cxkertype = "beta",
          cykertype = "gaussian",
          cxkerorder = 4L,
          cykerorder = 6L,
          cxkerbound = "fixed",
          cxkerlb = 0,
          cxkerub = 1
        ))

        got <- npRmpi:::.np_inid_boot_from_ksum_conditional_exact(
          xdat = xdat,
          ydat = ydat,
          exdat = exdat,
          eydat = eydat,
          bws = bw,
          B = ncol(counts),
          cdf = cdf,
          counts = counts
        )

        oracle <- vapply(seq_len(ncol(counts)), function(j) {
          idx <- npRmpi:::.np_counts_to_indices(counts[, j])
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
        }, numeric(nrow(exdat)))
        center <- do.call(estimator, list(
          bws = bw,
          txdat = xdat,
          tydat = ydat,
          exdat = exdat,
          eydat = eydat,
          gradients = FALSE,
          proper = FALSE
        ))

        label <- paste(if (cdf) "cdf" else "pdf", bwtype,
                       if (bernstein) "bernstein" else "raw")
        expect_equal(got$t0, as.numeric(center[[result.name]]),
                     tolerance = 0, info = label)
        expect_equal(got$t, t(oracle), tolerance = 0, info = label)
      }
    }
  }
})

test_that("conditional exact LP dispatch does not depend on matched kernel families", {
  helper <- paste(
    deparse(body(npRmpi:::.np_inid_boot_from_ksum_conditional_exact)),
    collapse = "\n"
  )
  lp.branch <- regexpr("if \\(general\\.lp\\)", helper)
  beta.branch <- regexpr("if \\(any\\.beta\\)", helper)

  expect_gt(lp.branch[[1L]], 0L)
  expect_gt(beta.branch[[1L]], lp.branch[[1L]])
})
