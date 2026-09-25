test_that("public conditional-density MADS enters the prepared owner on every rank", {
  skip_on_cran()
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  probe <- function() {
    library(npRmpi)
    npRmpi.init(nslaves = 1L, quiet = TRUE)
    on.exit(npRmpi.quit(), add = TRUE)
    ns <- asNamespace("npRmpi")
    get(".npRmpi_bcast_cmd_expr", ns)(quote({
      assign(".prepared_search_entries", 0L, .GlobalEnv)
      trace("npRmpiPreparedObjectiveFixedNativeSearchConditionalDensity",
        where = asNamespace("npRmpi"), print = FALSE, tracer = quote({
          assign(".prepared_search_entries",
            get(".prepared_search_entries", .GlobalEnv) + 1L, .GlobalEnv)
        }))
    }), caller.execute = TRUE)
    options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE)
    set.seed(167)
    x <- data.frame(x = runif(24L, -1, 1))
    y <- data.frame(y = rnorm(24L))
    for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
      b <- npcdensbw(xdat = x, ydat = y,
        bws = if (type == "fixed") c(.6, .5) else c(12, 12),
        bwtype = type, bwmethod = "cv.ls", regtype = "lc",
        nomad = FALSE, bandwidth.compute = FALSE)
      fit <- npcdensbw(xdat = x, ydat = y, bws = b, nmulti = 1L,
        bwsolver = "mads", nomad.opts = list(MAX_BB_EVAL = 4L),
        powell.remin = FALSE)
      stopifnot(is.finite(fit$fval), abs(fit$fval) < 1e100,
        fit$num.feval == 4)
    }
    counts <- get(".npRmpi_bcast_cmd_expr", ns)(quote(
      npRmpi:::mpi.allgather.Robj(get(".prepared_search_entries", .GlobalEnv))
    ), caller.execute = TRUE)
    stopifnot(length(counts) == 2L, all(unlist(counts) == 3L))
  }
  result <- npRmpi_run_isolated_contract(
    lines = c(paste0("(", paste(deparse(probe), collapse = "\n"), ")()"),
      "cat('CONDITIONAL_PUBLIC_MADS_POOL_OK\\n')"),
    marker = "CONDITIONAL_PUBLIC_MADS_POOL_OK", timeout = 20L)
  expect_false(is.null(result))
  if (!is.null(result)) {
    expect_identical(result$status, 0L, info = paste(result$output, collapse = "\n"))
    expect_true(result$witnessed)
  }
})
