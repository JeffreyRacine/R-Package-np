test_that("master-drawn wild gradient plots do not allocate unused task seeds", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_preserve_seed()
  withr::local_options(list(np.messages = FALSE, np.tree = FALSE))
  set.seed(719)
  n <- 40L
  x <- data.frame(x = runif(n), v = rnorm(n),
                  u = factor(rep(letters[1:2], n/2)))
  y <- sin(3*x$x) + .2*x$v + .2*(x$u == "b") + rnorm(n, sd = .2)
  counter <- new.env(parent = emptyenv())
  counter$seeds <- logical()
  trace(".npRmpi_bootstrap_chunk_tasks", where = asNamespace("npRmpi"),
    print = FALSE, tracer = bquote(assign("seeds",
      c(get("seeds", envir = .(counter)), isTRUE(with.seeds)),
      envir = .(counter))))
  on.exit(untrace(".npRmpi_bootstrap_chunk_tasks",
                 where = asNamespace("npRmpi")), add = TRUE)
  for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
    bw <- npregbw(xdat = x, ydat = y,
      bws = c(if (type == "fixed") c(.3, .6) else c(16,18), .3),
      bandwidth.compute = FALSE, bwtype = type, regtype = "lc")
    previous <- NULL
    for (chunk in c(3L, 5L)) {
      withr::local_options(list(np.plot.wild.chunk.size = chunk))
      set.seed(981)
      actual <- plot(bw, xdat = x, ydat = y, output = "data",
        errors = "bootstrap", plot.errors.boot.method = "wild",
        plot.errors.type = "pmzsd", B = 7L, neval = 5L, gradients = TRUE)
      values <- lapply(actual, function(z) z[c("grad", "gerr")])
      if (!is.null(previous)) expect_equal(values, previous, tolerance = 5e-10)
      previous <- values
    }
  }
  expect_gt(length(counter$seeds), 0L)
  expect_false(any(counter$seeds))
})
