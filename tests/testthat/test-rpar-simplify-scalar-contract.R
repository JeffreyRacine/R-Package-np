test_that(".simplify enforces one scalar logical-control contract", {
  simplify_result <- getFromNamespace(".simplify", "npRmpi")
  answer <- list(a = 1L, b = 2L)

  expect_identical(simplify_result(2L, answer, TRUE), c(a = 1L, b = 2L))
  expect_identical(simplify_result(2L, answer, FALSE), answer)
  expect_identical(simplify_result(2L, answer, 1), c(a = 1L, b = 2L))
  expect_identical(simplify_result(2L, answer, 0), answer)

  for (value in list(c(TRUE, FALSE), logical(), NA)) {
    expect_error(
      simplify_result(2L, answer, value),
      "'simplify' must be TRUE or FALSE",
      fixed = TRUE
    )
  }
})

test_that("simplify validation preserves established earlier condition precedence", {
  expect_error(
    npRmpi:::mpi.remote.exec(1L, simplify = c(TRUE, FALSE)),
    "It seems no slaves running.",
    fixed = TRUE
  )
  expect_error(
    npRmpi:::mpi.parSapply(1L, "not_a_function", job.num = 2L,
                           simplify = c(TRUE, FALSE)),
    "was not found",
    fixed = TRUE
  )
})

test_that("invalid remote-exec simplify never dispatches and preserves the pool", {
  skip_on_cran()
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "installed npRmpi unavailable for subprocess MPI contract")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(np.messages = FALSE)",
      "trace <- tempfile('npRmpi-simplify-trace-', fileext='.tsv')",
      "options(npRmpi.transport.trace.file = trace)",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "condition_text <- function(expr) tryCatch({ force(expr); NA_character_ }, error=conditionMessage)",
      "for (value in list(c(TRUE, FALSE), logical(), NA)) {",
      "  before <- if (file.exists(trace)) readLines(trace, warn=FALSE) else character()",
      "  msg <- condition_text(npRmpi:::mpi.remote.exec(1L, simplify=value))",
      "  stopifnot(identical(msg, \"'simplify' must be TRUE or FALSE\"))",
      "  after <- if (file.exists(trace)) readLines(trace, warn=FALSE) else character()",
      "  stopifnot(identical(after, before))",
      "  ping <- npRmpi:::mpi.remote.exec(2L, simplify=TRUE)",
      "  stopifnot(is.data.frame(ping), identical(unname(unlist(ping)), 2L))",
      "}",
      "cat('SIMPLIFY_REMOTE_NO_DISPATCH_POOL_REUSE_OK\\n')"
    ),
    timeout = 45L,
    env = env
  )

  info <- paste(res$output, collapse = "\n")
  expect_true(res$status %in% c(0L, 137L), info = info)
  expect_true(any(grepl("SIMPLIFY_REMOTE_NO_DISPATCH_POOL_REUSE_OK",
                        res$output, fixed = TRUE)), info = info)
})

test_that("parallel simplify controls reject before work and preserve three-slave reuse", {
  skip_on_cran()
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "installed npRmpi unavailable for subprocess MPI contract")

  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "options(np.messages = FALSE)",
      "npRmpi.init(nslaves = 3L, quiet = TRUE)",
      "condition_text <- function(expr) tryCatch({ force(expr); NA_character_ }, error=conditionMessage)",
      "expected <- \"'simplify' must be TRUE or FALSE\"",
      "prefix <- tempfile('npRmpi-simplify-side-effect-')",
      "side_fun <- function(x, prefix) { file.create(paste0(prefix, '-', mpi.comm.rank())); x }",
      "msg <- condition_text(npRmpi:::mpi.parSapply(1:6, side_fun, prefix=prefix, simplify=c(TRUE, FALSE)))",
      "stopifnot(identical(msg, expected), !length(Sys.glob(paste0(prefix, '-*'))))",
      "stopifnot(identical(unname(npRmpi:::mpi.parSapply(1:6, identity)), 1:6))",
      "msg <- condition_text(npRmpi:::mpi.iparSapply(1:6, side_fun, prefix=prefix, simplify=logical()))",
      "stopifnot(identical(msg, expected), !length(Sys.glob(paste0(prefix, '-*'))))",
      "stopifnot(identical(unname(npRmpi:::mpi.iparSapply(1:6, identity)), 1:6))",
      "msg <- condition_text(npRmpi:::mpi.parReplicate(6L, { file.create(paste0(prefix, '-', mpi.comm.rank())); 1L }, simplify=NA))",
      "stopifnot(identical(msg, expected), !length(Sys.glob(paste0(prefix, '-*'))))",
      "stopifnot(identical(unname(npRmpi:::mpi.parReplicate(6L, 7L)), rep.int(7L, 6L)))",
      "msg <- condition_text(npRmpi:::mpi.iparReplicate(6L, { file.create(paste0(prefix, '-', mpi.comm.rank())); 1L }, simplify=c(FALSE, TRUE)))",
      "stopifnot(identical(msg, expected), !length(Sys.glob(paste0(prefix, '-*'))))",
      "stopifnot(identical(unname(npRmpi:::mpi.iparReplicate(6L, 8L)), rep.int(8L, 6L)))",
      "bad_rand <- function(n) { file.create(paste0(prefix, '-', mpi.comm.rank())); rep.int(1, n) }",
      "msg <- condition_text(npRmpi:::mpi.parSim(n=2L, rand.gen=bad_rand, statistic=sum, nsim=1L, simplify=NA))",
      "stopifnot(identical(msg, expected), !length(Sys.glob(paste0(prefix, '-*'))))",
      "ping <- npRmpi:::mpi.remote.exec(9L)",
      "stopifnot(is.data.frame(ping), identical(unname(unlist(ping)), rep.int(9L, 3L)))",
      "cat('SIMPLIFY_PARALLEL_NO_WORK_POOL_REUSE_OK\\n')"
    ),
    timeout = 60L,
    env = env
  )

  info <- paste(res$output, collapse = "\n")
  expect_true(res$status %in% c(0L, 137L), info = info)
  expect_true(any(grepl("SIMPLIFY_PARALLEL_NO_WORK_POOL_REUSE_OK",
                        res$output, fixed = TRUE)), info = info)
})
