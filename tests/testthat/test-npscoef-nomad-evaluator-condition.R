test_that("opted-in NOMAD preserves the original evaluator condition", {
  skip_on_cran()
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  original <- structure(list(
    message = paste(rep("original evaluator failure", 70L), collapse = ":"),
    call = quote(scoef_leaf(x = 1L)), token = "original",
    values = c(1.25, -7), detail = list(index = 2L)
  ), class = c("scoef_test_error", "error", "condition"))
  visits <- payloads <- 0L
  driver <- function(active = TRUE, condition = NULL) {
    visits <<- 0L
    payloads <<- 0L
    .np_nomad_search(
      engine = "nomad", baseline_record = NULL, start_degree = 1L,
      x0 = 0, bbin = 1L, lb = 0, ub = 4,
      eval_fun = function(point) {
        visits <<- visits + 1L
        if (visits == 2L && !is.null(condition)) stop(condition)
        list(objective = (point[1L] - 3)^2 + 1, degree = 1L, num.feval = 1L)
      },
      build_payload = function(point, best_record, solution, interrupted) {
        payloads <<- payloads + 1L
        list(payload = list(point = point, objective = best_record$objective))
      },
      native.r.bridge = TRUE, preserve.eval.error = active,
      nmulti = 2L, remin = TRUE, random.seed = 42L,
      nomad.opts = list(MAX_BB_EVAL = 20L)
    )
  }
  got <- tryCatch(driver(condition = original), error = identity)
  expect_identical(got, original)
  expect_identical(visits, 2L)
  expect_identical(payloads, 0L)
  # A fresh invocation has its own condition state.
  good <- driver()
  expect_identical(as.numeric(good$best_point), 3)
  expect_identical(payloads, 1L)
  typed <- original
  class(typed) <- c("np_nn_candidate_invalid", "error", "condition")
  expect_false(inherits(driver(condition = typed), "condition"))
  expect_gt(visits, 2L)
  expect_identical(payloads, 1L)
  expect_false(inherits(driver(active = FALSE, condition = original), "condition"))
  expect_gt(visits, 2L)
})

test_that("owned-pool error classification leaves direct invalidity policy intact", {
  original <- structure(list(message = "original", call = quote(moment()), tag = 7L),
                        class = c("scoef_test_error", "error", "condition"))
  typed <- original
  class(typed) <- c("np_nn_candidate_invalid", "error", "condition")
  classify <- .npscoefbw_nomad_unknown_nn_error
  for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
    bws <- list(type = type)
    expect_identical(classify(original, bws, preserve.eval.error = TRUE), original)
    expect_null(classify(typed, bws, preserve.eval.error = TRUE))
    expect_null(classify(typed, bws))
    if (type == "fixed") expect_null(classify(original, bws))
    else expect_identical(classify(original, bws), original)
  }
})

test_that("pool policy requires requested and actual root ownership", {
  ns <- asNamespace("npRmpi")
  active <- TRUE
  broadcast <- FALSE
  rank <- 0L
  size <- 2L
  sends <- 0L
  command <- NULL
  old <- options(npRmpi.local.regression.mode = FALSE)
  on.exit(options(old), add = TRUE)
  testthat::with_mocked_bindings({
    start <- function(request = TRUE)
      .npscoefbw_nomad_pool_start(list(), preserve.eval.error = request)
    good <- start()
    expect_true(good[["preserve.eval.error", exact = TRUE]])
    expect_identical(good$nslaves, 1L)
    expect_match(paste(deparse(command), collapse = ""), "PRESERVE_EVAL_ERROR = TRUE", fixed = TRUE)
    expect_false(start(FALSE)[["preserve.eval.error", exact = TRUE]])
    rank <- 1L
    expect_false(start()[["preserve.eval.error", exact = TRUE]])
    rank <- NA_integer_
    expect_false(start()[["preserve.eval.error", exact = TRUE]])
    rank <- 0L
    size <- 1L
    expect_false(start()[["preserve.eval.error", exact = TRUE]])
    size <- 2L
    prior <- sends
    active <- FALSE
    expect_null(start())
    active <- TRUE
    broadcast <- TRUE
    expect_null(start())
    broadcast <- FALSE
    options(npRmpi.local.regression.mode = TRUE)
    expect_null(start())
    expect_identical(sends, prior)
  }, .npRmpi_has_active_slave_pool = function(...) active,
     .npRmpi_autodispatch_called_from_bcast = function(...) broadcast,
     .npRmpi_autodispatch_remote_ref = function(...) "diagnostic-context",
     mpi.comm.size = function(...) size,
     mpi.comm.rank = function(...) rank,
     .npRmpi_bcast_cmd_expr = function(expr, ...) {
       sends <<- sends + 1L
       command <<- expr
       invisible(NULL)
     }, .package = "npRmpi")
})
