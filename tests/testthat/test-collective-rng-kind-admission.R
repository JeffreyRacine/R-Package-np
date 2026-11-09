test_that("common-kind admission is opt-in and examines the full certified tuple", {
  capture <- getFromNamespace(".np_rng_scope_capture", "npRmpi")
  e <- new.env(parent = environment(capture))
  environment(capture) <- e
  e$k <- 10403L
  e$mismatch <- FALSE
  e$calls <- 0L
  e$phases <- character()
  e$.np_rng_scope_capture_local <- function(uniform.only) list(seed = c(e$k, 7L))
  e$mpi.allreduce <- function(x, type, op) {
    e$calls <- e$calls + 1L
    expect_identical(x, c(e$k, -e$k))
    expect_identical(type, 1L)
    expect_identical(op, "max")
    if (e$mismatch) c(e$k + 1L, -e$k) else x
  }
  e$.npRmpi_rng_scope_agree <- function(failure, phase, deferred) {
    e$phases <- c(e$phases, phase)
    if (is.null(failure)) return(NULL)
    refusal <- structure(list(message = conditionMessage(failure), call = NULL),
                         class = c("np_rng_scope_refusal", "error", "condition"))
    if (deferred) refusal else stop(refusal)
  }
  capture(collective = TRUE)
  capture(same.kind = TRUE)
  expect_identical(e$calls, 0L)
  for (header in c(10403L, 10407L, 10103L, 403L)) {
    e$k <- header
    scope <- capture(collective = TRUE, same.kind = TRUE)
    expect_identical(scope$seed, c(header, 7L))
    expect_true(scope$collective)
  }
  expect_identical(e$calls, 4L)
  e$mismatch <- TRUE
  expect_error(capture(collective = TRUE, same.kind = TRUE),
               "different RNG kinds", class = "np_rng_scope_refusal")
  expect_s3_class(capture(collective = TRUE, same.kind = TRUE, deferred = TRUE),
                  "np_rng_scope_refusal")
  expect_identical(tail(e$phases, 2L), c("admission", "kind agreement"))
})

test_that("rank-partitioned public tests opt in while master-owned plans do not", {
  ns <- asNamespace("npRmpi")
  for (nm in c("npcmstest", "npqcmstest", "npdeneqtest",
               "npdeptest", "npsdeptest", "npsymtest", "npunitest")) {
    f <- get(nm, ns)
    calls <- list()
    walk <- function(x) {
      if (missing(x)) return(invisible(NULL))
      if (is.call(x) && identical(x[[1L]], as.name(".np_rng_scope_capture")))
        calls[[length(calls) + 1L]] <<- x
      if (is.call(x) || is.pairlist(x) || is.expression(x))
        for (a in as.list(x)) walk(a)
    }
    walk(body(f))
    expect_length(calls, 1L)
    expect_identical(calls[[1L]][["same.kind"]], TRUE, info = nm)
  }
  expect_false(grepl("same.kind = TRUE", paste(deparse(body(get("npsigtest.rbandwidth", ns))), collapse = " "), fixed = TRUE))
})
