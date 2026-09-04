test_that("all seven cell-search owners propagate an evaluator failure", {
  had.pool <- .mpi_pool_active()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  n <- 18L; i <- seq_len(n)
  d <- data.frame(x = sin(i * sqrt(2)) + i / 30,
                  z = cos(i * sqrt(3)) + i / 50)
  d$y <- 1 + d$x * (1 + d$z) + sin(i) / 10
  original <- getFromNamespace(".np_degree_search", "npRmpi")
  marker <- errorCondition("diagnostic owner evaluator failure",
    class = "np_test_owner_error", call = quote(fixed_degree_owner()), token = 808L)
  calls <- 0L; valid <- 0L
  caught <- NULL
  wrapped <- function(..., eval_fun) {
    calls <<- calls + 1L
    evaluate <- function(degree) {
      if (degree[[1L]] == 1L) stop(marker)
      out <- eval_fun(degree)
      valid <<- valid + 1L
      out
    }
    tryCatch(original(..., eval_fun = evaluate), error = function(e) {
      caught <<- e
      stop(e)
    })
  }
  # npindexbw uses a collective degree driver with more than one slave.
  # A master-only injected error would deliberately desynchronize that driver.
  # Install the same failure before native entry on every worker, preserving
  # and restoring each worker's own original binding without global objects.
  mpi.bcast.cmd(local({
    original <- getFromNamespace(".np_degree_search", "npRmpi")
    marker <- errorCondition("diagnostic owner evaluator failure",
      class = "np_test_owner_error", call = quote(fixed_degree_owner()), token = 808L)
    probe <- function(..., eval_fun) {
      evaluate <- function(degree) {
        if (degree[[1L]] == 1L) stop(marker)
        eval_fun(degree)
      }
      original(..., eval_fun = evaluate)
    }
    attr(probe, "np.test.original") <- original
    assignInNamespace(".np_degree_search", probe, "npRmpi")
  }), caller.execute = FALSE)
  on.exit({
    mpi.bcast.cmd(local({
      probe <- getFromNamespace(".np_degree_search", "npRmpi")
      assignInNamespace(".np_degree_search", attr(probe, "np.test.original"), "npRmpi")
    }), caller.execute = FALSE)
    if (!had.pool) close_mpi_slaves(force = TRUE)
  }, add = TRUE)
  local_mocked_bindings(.np_degree_search = wrapped, .package = "npRmpi")
  for (owner in c("npregbw", "npcdensbw", "npcdistbw", "npscoefbw",
                  "npplregbw", "npindexbw", "nplsqregbw")) {
    args <- list(formula = if (owner %in% c("npscoefbw", "npplregbw"))
                   y ~ x | z else y ~ x,
                 data = d, regtype = "lp", nomad = FALSE,
                 degree.select = "exhaustive", search.engine = "cell",
                 degree.min = 0L, degree.max = 2L, degree.start = 0L,
                 nmulti = 1L)
    if (owner %in% c("npindexbw", "npscoefbw")) args$optim.maxit <- 20L
    else {args$itmax <- 2L; args$powell.remin <- FALSE}
    if (owner == "npcdistbw") args$ngrid <- 5L
    if (owner == "npindexbw") args$formula <- y ~ x + z
    if (owner == "nplsqregbw") {
      args$bws <- args$formula
      args$formula <- NULL
      args$nomad <- "auto"
      args$scale <- rep(1, n)
    }
    npseed(934)
    before <- calls
    caught <- NULL
    condition <- tryCatch(do.call(getFromNamespace(owner, "npRmpi"), args), error = identity)
    # Preserve the exact condition at the shared driver. Existing outer SPMD
    # transport may wrap it; this tranche intentionally does not change that.
    expect_identical(caught, marker, info = owner)
    if (owner == "npindexbw" && mpi.comm.size(1L) > 2L) {
      expect_s3_class(condition, "error")
      expect_match(conditionMessage(condition), conditionMessage(marker), fixed = TRUE)
    } else expect_identical(condition, marker, info = owner)
    expect_equal(calls, before + 1L, info = owner)
  }
  expect_true(is.finite(npregbw(y ~ x, data = d, regtype = "ll",
    nomad = FALSE, nmulti = 1L, itmax = 1L)$fval))
})
