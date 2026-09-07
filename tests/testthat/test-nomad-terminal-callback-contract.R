t2_terminal_condition <- function(class = "t2_leaf_error") {
  structure(list(message = "terminal callback contract", call = quote(t2_leaf()),
                 token = list(value = 733L)),
            class = c(class, "error", "condition"))
}

t2_terminal_driver <- function(eval_fun, build_payload, ...) {
  getFromNamespace(".np_nomad_search", "npRmpi")(
    engine = "nomad", baseline_record = NULL, start_degree = 1L,
    x0 = 0, bbin = 1L, lb = 0, ub = 4,
    eval_fun = eval_fun, build_payload = build_payload,
    native.r.bridge = TRUE, preserve.eval.error = TRUE,
    random.seed = 42L, nomad.opts = list(MAX_BB_EVAL = 20L), ...)
}

test_that("terminal callbacks cannot evaluate or publish after an unexpected failure", {
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  for (stage in c("first", "interior", "last", "input", "result", "admissible")) {
    visits <- payloads <- requests <- 0L
    original <- t2_terminal_condition()
    failure.at <- switch(stage, first = 1L, interior = 2L, last = 4L, 2L)
    got <- testthat::with_mocked_bindings({
      tryCatch(t2_terminal_driver(
        function(point) {
          visits <<- visits + 1L
          if (visits == failure.at) {
            if (stage %in% c("first", "interior", "last")) stop(original)
            if (stage == "result")
              return(list(objective = list(new.env()), degree = 1L))
            if (stage == "admissible")
              return(list(objective = 1, degree = 1L, admissible = NA))
          }
          list(objective = point[1L]^2 + 1, degree = 1L)
        },
        function(...) { payloads <<- payloads + 1L; list(payload = TRUE) },
        nmulti = 2L, remin = TRUE), error = identity)
    }, .np_nomad_native_r_callback_search = function(eval.f, ...) {
      # Deliberately request more callbacks even after native-caught errors.
      for (i in 1:6) {
        requests <<- requests + 1L
        point <- if (stage == "input" && i == 2L) list(new.env()) else (i - 1) %% 5
        tryCatch(eval.f(point), error = function(e) NULL)
      }
      list(value = list(status = "error", native_status = 2L,
                        message = "generic native status"), output = character())
    }, .package = "npRmpi")
    expect_s3_class(got, "error")
    if (stage %in% c("first", "interior", "last")) expect_identical(got, original)
    expect_identical(visits, if (stage == "input") 1L else failure.at, info = stage)
    expect_identical(payloads, 0L, info = stage)
    expect_identical(requests, 6L, info = stage)
    expect_false(identical(conditionMessage(got), "generic native status"))
  }
})

test_that("rejection stays exploratory and required payload failures propagate", {
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  visits <- payloads <- 0L
  rejection <- t2_terminal_condition("np_nn_candidate_invalid")
  got <- t2_terminal_driver(function(point) {
    visits <<- visits + 1L
    if (visits == 2L) stop(rejection)
    list(objective = (point[1L] - 3)^2 + 1, degree = 1L)
  }, function(point, ...) { payloads <<- payloads + 1L; list(payload = point) },
  nmulti = 2L, remin = TRUE)
  expect_gt(visits, 2L)
  expect_identical(payloads, 1L)
  expect_equal(as.numeric(got$best_point), 3)
  # The same typed rejection is terminal in required payload/certification.
  got <- tryCatch(t2_terminal_driver(function(point)
    list(objective = point[1L]^2 + 1, degree = 1L),
    function(...) stop(rejection)), error = identity)
  expect_identical(got, rejection)
})

test_that("interrupt cannot publish an incumbent or enter another restart", {
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  for (interrupt.solve in 1:2) {
    visits <- payloads <- solves <- 0L
    original <- structure(list(message = "test interrupt", call = NULL),
                          class = c("interrupt", "condition"))
    got <- testthat::with_mocked_bindings({
      tryCatch(t2_terminal_driver(function(point) {
        visits <<- visits + 1L
        if (solves == interrupt.solve && point == 1) stop(original)
        list(objective = (point - 3)^2 + 1, degree = 1L)
      }, function(...) { payloads <<- payloads + 1L; list(payload = TRUE) },
      nmulti = 1L, remin = TRUE), interrupt = identity)
    }, .np_nomad_native_r_callback_search = function(eval.f, ...) {
      solves <<- solves + 1L
      for (i in 0:3) tryCatch(eval.f(i), interrupt = function(e) NULL)
      list(value = list(status = "ok", native_status = 0L, solution = 3,
                        objective = 1, message = "ok"), output = character())
    }, .package = "npRmpi")
    expect_identical(got, original)
    expect_identical(visits, (interrupt.solve - 1L) * 4L + 2L)
    expect_identical(payloads, 0L)
    expect_identical(solves, interrupt.solve)
  }
})

test_that("abort rendering cannot replace the original terminal condition", {
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  for (stage in c("callback", "payload")) {
    visits <- aborts <- 0L
    original <- t2_terminal_condition()
    got <- testthat::with_mocked_bindings({
      tryCatch(t2_terminal_driver(function(point) {
        visits <<- visits + 1L
        if (stage == "callback" && visits == 2L) stop(original)
        list(objective = point[1L]^2 + 1, degree = 1L)
      }, function(...) stop(original)), error = identity)
    }, .np_progress_abort = function(...) {
      aborts <<- aborts + 1L
      stop("diagnostic rendering failed")
    }, .package = "npRmpi")
    expect_identical(got, original)
    expect_identical(aborts, 1L)
  }
})

test_that("index service phase agreement preserves conditions and healthy values", {
  transaction <- getFromNamespace(".npindexbw_service_transaction", "npRmpi")
  original <- t2_terminal_condition()
  for (case in c("healthy", "skip", "native.unwind", "prepare.root", "returned.worker", "prepare.worker")) {
    rank <- if (case == "prepare.worker") 1L else 0L
    ctx <- list(rank = rank, size = 2L, comm = 1L, root = rank == 0L)
    replies <- switch(case, healthy = list(c(0L, 2L), 0L), skip = list(c(0L, 0L)),
      native.unwind = list(c(0L, 2L)),
      prepare.root = list(c(1L, 1L), 0L), returned.worker = list(c(0L, 2L), 1L, 1L),
      prepare.worker = list(c(1L, 1L), 1L))
    calls <- broadcasts <- 0L
    got <- testthat::with_mocked_bindings({
      tryCatch(transaction(ctx, function(guard) {
        if (startsWith(case, "prepare")) stop(original)
        if (case == "skip") return(7)
        guard$enter(native = TRUE)
        if (case == "native.unwind") stop(original)
        guard$phase <- "returned"
        7
      }), error = identity)
    }, mpi.allreduce = function(...) {
      calls <<- calls + 1L
      replies[[calls]]
    }, mpi.bcast.Robj = function(obj, rank, comm) {
      broadcasts <<- broadcasts + 1L
      # The existing sender returns NULL; only receivers return the object.
      if (ctx$rank == rank) NULL else original
    }, .package = "npRmpi")
    expect_identical(calls, length(replies), info = case)
    expect_identical(broadcasts, if (case %in% c("healthy", "skip", "native.unwind")) 0L else 1L)
    if (case %in% c("healthy", "skip")) expect_identical(got, 7)
    else if (case == "prepare.worker") expect_identical(got$service.error, original)
    else expect_identical(got, original)
  }
})
