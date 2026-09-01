test_that("npconmode level batching matches independent density evaluations", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260811)
  n <- 48L
  x <- seq(-0.9, 0.9, length.out = n)
  y <- factor(sample(c("a", "b", "c"), n, replace = TRUE),
              levels = c("a", "b", "c"))
  xdat <- data.frame(x = x)
  ydat <- data.frame(y = y)
  bw <- npcdensbw(
    xdat = xdat, ydat = ydat, bws = c(0.12, 0.42),
    bandwidth.compute = FALSE, bwtype = "fixed",
    cxkertype = "epanechnikov", regtype = "ll"
  )
  efac <- factor(levels(y), levels = levels(y))

  batched <- getFromNamespace(".npConmodeEvaluateLevels", "npRmpi")(
    bws = bw, txdat = xdat, tydat = ydat, xeval = xdat,
    efac = efac, gradients = TRUE, gradient.level.index = 2L
  )

  legacy <- lapply(seq_along(efac), function(i) {
    npcdens(
      bws = bw, txdat = xdat, tydat = ydat, exdat = xdat,
      eydat = rep(efac[i], n), gradients = i == 2L
    )
  })
  legacy.p <- do.call(cbind, lapply(legacy, `[[`, "condens"))
  legacy.e <- do.call(cbind, lapply(legacy, `[[`, "conderr"))

  expect_equal(batched$probabilities, legacy.p, tolerance = 1e-13,
               ignore_attr = TRUE)
  expect_equal(batched$errors, legacy.e, tolerance = 1e-13,
               ignore_attr = TRUE)
  expect_equal(batched$gradients, legacy[[2L]]$congrad, tolerance = 1e-13,
               ignore_attr = TRUE)
})

test_that("npconmode batches ordinary levels around one gradient singleton", {
  calls <- list()
  local_mocked_bindings(
    npcdens = function(bws, txdat, tydat, exdat, eydat, gradients, ...) {
      calls[[length(calls) + 1L]] <<- list(
        levels = unique(as.character(eydat)),
        gradients = gradients
      )
      n <- nrow(exdat)
      list(
        condens = seq_len(n),
        conderr = seq_len(n),
        congrad = if (gradients) matrix(0, nrow = n, ncol = bws$xndim)
      )
    },
    .package = "npRmpi"
  )

  efac <- factor(c("a", "b", "c"), levels = c("a", "b", "c"))
  args <- list(
    bws = list(xndim = 1L, xnames = "x"),
    txdat = data.frame(x = 1:4),
    tydat = data.frame(y = efac[c(1L, 2L, 3L, 1L)]),
    xeval = data.frame(x = 1:4),
    efac = efac,
    gradient.level.index = 2L
  )
  do.call(getFromNamespace(".npConmodeEvaluateLevels", "npRmpi"),
          c(args, list(gradients = TRUE)))

  expect_length(calls, 2L)
  expect_identical(calls[[1L]], list(levels = c("a", "c"), gradients = FALSE))
  expect_identical(calls[[2L]], list(levels = "b", gradients = TRUE))

  calls <- list()
  do.call(getFromNamespace(".npConmodeEvaluateLevels", "npRmpi"),
          c(args, list(gradients = FALSE)))
  expect_length(calls, 1L)
  expect_identical(calls[[1L]],
                   list(levels = c("a", "b", "c"), gradients = FALSE))
})

test_that("npconmode level batching has a bounded temporary row plan", {
  planner <- getFromNamespace(".npConmodeLevelBlockWidth", "npRmpi")

  expect_identical(planner(2048L, 2L, 5L), 5L)
  expect_identical(planner(1000000L, 2L, 5L), 1L)
  expect_lte(planner(1L, 100L, 100000L), 65536L)
  expect_error(planner(0L, 2L, 5L), "invalid conditional-mode")
})
