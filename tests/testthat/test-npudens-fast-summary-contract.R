test_that("npudens exposes unconditional fast counts in bandwidth summary", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(42)
  n <- 100L
  y <- runif(n)

  fit <- npudens(
    ~ y,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1
  )

  expect_true(is.finite(as.numeric(fit$bws$num.feval[1L])))
  expect_true(is.finite(as.numeric(fit$bws$num.feval.fast[1L])))
  expect_gt(as.numeric(fit$bws$num.feval.fast[1L]), 0)
  expect_lte(as.numeric(fit$bws$num.feval.fast[1L]), as.numeric(fit$bws$num.feval[1L]))

  txt <- paste(capture.output(summary(fit$bws)), collapse = "\n")
  expect_true(grepl("Number of Function Evaluations:", txt, fixed = TRUE))

  fast <- as.numeric(fit$bws$num.feval.fast[1L])
  cache.hits <- if (!is.null(fit$bws$nn.cache) &&
                    "objective.hits" %in% names(fit$bws$nn.cache)) {
    as.numeric(fit$bws$nn.cache[["objective.hits"]])
  } else {
    0
  }

  if (is.finite(cache.hits) && cache.hits > 0)
    expect_true(grepl("Evaluation cache (Powell):", txt, fixed = TRUE))

  fast.extra <- fast - if (is.finite(cache.hits)) cache.hits else 0
  if (is.finite(fast.extra) && fast.extra > 0) {
    expect_true(grepl(
      sprintf("Fast CV route: %s of", format(round(fast.extra), big.mark = ",")),
      txt,
      fixed = TRUE
    ))
  } else {
    expect_false(grepl("Fast CV route:", txt, fixed = TRUE))
  }
})
