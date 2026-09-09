b2_condition <- function(class = "b2_preparation_error") {
  structure(list(message = "B2 original local failure", call = quote(b2_leaf()),
                 token = list(value = 642L)),
            class = c(class, "error", "condition"))
}

test_that("B2 joins preparation and completion and preserves the first peer cause", {
  original <- b2_condition()
  for (case in c("healthy", "nested", "skip", "prepare.root", "prepare.peer",
                 "returned.root", "returned.peer", "native.unwind")) {
    replies <- switch(case,
      healthy =, nested = list(c(0L, 0L, 0L, 2L), c(0L, 0L)),
      skip = list(c(0L, 0L, 0L, 0L), c(0L, 0L)),
      prepare.root = list(c(1L, 0L, 0L, 1L), 0L),
      prepare.peer = list(c(1L, 0L, 0L, 1L), 1L),
      returned.root = list(c(0L, 0L, 0L, 2L), c(1L, 0L), 0L),
      returned.peer = list(c(0L, 0L, 0L, 2L), c(1L, 0L), 1L),
      native.unwind = list(c(0L, 0L, 0L, 2L)))
    counts <- new.env(parent = emptyenv())
    counts$reductions <- counts$broadcasts <- counts$evaluations <- 0L
    testthat::with_mocked_bindings({
      state <- .npscoefbw_collective_state(comm = 7L)
      evaluate <- function() {
        counts$evaluations <- counts$evaluations + 1L
        if (case == "prepare.root") stop(original)
        if (case == "skip") return(7)
        state$guard$enter(native = TRUE)
        if (case == "native.unwind") stop(original)
        state$guard$phase <- "returned"
        if (case == "returned.root") stop(original)
        7
      }
      got <- tryCatch(.npscoefbw_collective_transaction(state, function() {
        if (case == "nested")
          .npscoefbw_collective_transaction(state, evaluate) else evaluate()
      }), error = identity)
      expect_identical(counts$reductions, length(replies), info = case)
      expect_identical(counts$evaluations, 1L, info = case)
      expect_null(state$guard)
      if (case %in% c("healthy", "nested", "skip")) {
        expect_identical(got, 7, info = case)
        expect_null(state$terminal)
      } else {
        expect_identical(got, original, info = case)
        if (case == "native.unwind") {
          expect_null(state$terminal)
        } else {
          expect_identical(state$terminal, original)
          # The search latch prevents all later science and agreements.
          again <- tryCatch(.npscoefbw_collective_transaction(state, evaluate),
                            error = identity)
          expect_identical(again, original)
          expect_identical(counts$evaluations, 1L)
          expect_identical(counts$reductions, length(replies))
        }
      }
      expect_identical(counts$broadcasts,
        if (case %in% c("healthy", "nested", "skip", "native.unwind")) 0L else 1L)
    }, mpi.comm.rank = function(...) 0L, mpi.comm.size = function(...) 2L,
       mpi.allreduce = function(x, type, op, comm) {
         expect_identical(comm, 7L)
         counts$reductions <- counts$reductions + 1L
         replies[[counts$reductions]]
       }, mpi.bcast.Robj = function(obj, rank, comm) {
         expect_identical(comm, 7L)
         counts$broadcasts <- counts$broadcasts + 1L
         if (rank == 0L) NULL else original
       }, .npRmpi_autodispatch_in_context = function() FALSE,
       .package = "npRmpi")
  }
})

test_that("B2 typed preparation rejection meets peers at native readiness", {
  typed <- b2_condition("np_nn_candidate_invalid")
  original <- b2_condition()
  for (case in c("typed.local", "typed.peer", "direct.invalid", "error.wins")) {
    replies <- switch(case,
      typed.local = list(c(0L, 1L, 0L, 1L), 0L, c(0L, 0L)),
      typed.peer = list(c(0L, 1L, 0L, 1L), 1L, c(0L, 0L)),
      direct.invalid = list(c(0L, 0L, 1L, 1L), c(0L, 0L)),
      error.wins = list(c(1L, 1L, 0L, 0L), 1L))
    counts <- new.env(parent = emptyenv())
    counts$reductions <- counts$native <- 0L
    testthat::with_mocked_bindings({
      state <- .npscoefbw_collective_state()
      got <- tryCatch(.npscoefbw_collective_transaction(state, function() {
        # The generic callback consumes only an agreed exploratory rejection.
        tryCatch({
          entered <- if (case %in% c("typed.local", "error.wins"))
            state$guard$enter(native = FALSE, error = typed) else
            state$guard$enter(native = TRUE)
          if (!entered) return("existing invalid objective")
          counts$native <- counts$native + 1L
          stop("must not enter native")
        }, np_nn_candidate_invalid = function(e) "existing rejection record")
      }), error = identity)
      expect_identical(counts$native, 0L)
      expect_identical(counts$reductions, length(replies), info = case)
      if (case == "error.wins") expect_identical(got, original)
      else expect_identical(got, if (case == "direct.invalid")
        "existing invalid objective" else "existing rejection record")
    }, mpi.comm.rank = function(...) 0L, mpi.comm.size = function(...) 2L,
       mpi.allreduce = function(...) {
         counts$reductions <- counts$reductions + 1L
         replies[[counts$reductions]]
       }, mpi.bcast.Robj = function(obj, rank, comm) {
         if (rank == 0L) NULL else if (case == "error.wins") original else typed
       }, .npRmpi_autodispatch_in_context = function() FALSE,
       .package = "npRmpi")
  }
})

test_that("B2 same-class errors after native return are not candidate rejection", {
  typed <- b2_condition("np_nn_candidate_invalid")
  for (raw.recovery in c(FALSE, TRUE)) for (consumed in c(FALSE, TRUE)) {
    counts <- new.env(parent = emptyenv())
    counts$reductions <- 0L
    replies <- list(c(0L, 0L, 0L, 2L), c(1L, 0L), 0L)
    testthat::with_mocked_bindings({
      state <- .npscoefbw_collective_state()
      got <- tryCatch(.npscoefbw_collective_transaction(state, function() {
        state$guard$enter(native = TRUE)
        state$guard$phase <- "returned"
        if (consumed) {
          # The direct-family catch retained the actual returned-phase error;
          # later generic evaluation may have consumed its exploratory class.
          state$guard$error <- typed
          return(7)
        }
        stop(typed)
      }, allow.typed.rejection = raw.recovery), error = identity)
      expect_identical(got, typed)
      expect_identical(state$terminal, typed)
      expect_identical(counts$reductions, 3L)
    }, mpi.comm.rank = function(...) 0L, mpi.comm.size = function(...) 2L,
       mpi.allreduce = function(...) {
         counts$reductions <- counts$reductions + 1L
         replies[[counts$reductions]]
       }, mpi.bcast.Robj = function(...) NULL,
       .npRmpi_autodispatch_in_context = function() FALSE,
       .package = "npRmpi")
  }
})

test_that("opted-in whole callbacks latch result and recording failures", {
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  for (stage in c("evaluate", "result", "typed.record")) {
    counts <- new.env(parent = emptyenv())
    counts$visits <- counts$transactions <- counts$records <- counts$payloads <- 0L
    original <- b2_condition(if (stage == "typed.record")
      "np_nn_candidate_invalid" else "b2_preparation_error")
    got <- testthat::with_mocked_bindings({
      tryCatch(.np_nomad_search(
        engine = "nomad+powell", baseline_record = NULL, start_degree = 1L,
        x0 = 0, bbin = 1L, lb = 0, ub = 4,
        eval_fun = function(point) {
          counts$visits <- counts$visits + 1L
          if (counts$visits == 2L) {
            if (stage == "evaluate") stop(original)
            if (stage == "result")
              return(list(objective = list(new.env()), degree = 1L))
          }
          list(objective = point[1L]^2 + 1, degree = 1L)
        },
        build_payload = function(...) {
          counts$payloads <- counts$payloads + 1L
          list(payload = TRUE)
        }, native.r.bridge = TRUE, preserve.eval.error = TRUE,
        .native.callback.transaction = function(evaluate, point) {
          counts$transactions <- counts$transactions + 1L
          evaluate()
        }, nmulti = 2L, remin = TRUE,
        random.seed = 42L, nomad.opts = list(MAX_BB_EVAL = 20L)), error = identity)
    }, .np_nomad_native_r_callback_search = function(eval.f, ...) {
      # Simulate an optimizer requesting callbacks after its first caught error.
      for (i in 0:5) tryCatch(eval.f(i %% 5), error = function(e) NULL)
      list(value = list(status = "error", native_status = 2L,
                        message = "native marker"), output = character())
    }, .np_degree_progress_step = function(state, ...) {
      counts$records <- counts$records + 1L
      if (stage == "typed.record" && counts$visits == 2L) stop(original)
      state
    }, .package = "npRmpi")
    expect_s3_class(got, "error")
    if (stage != "result") expect_identical(got, original)
    expect_identical(counts$visits, 2L, info = stage)
    expect_identical(counts$transactions, 2L, info = stage)
    expect_identical(counts$payloads, 0L, info = stage)
  }
})

test_that("B2 readiness restores actual native locality before scientific entry", {
  for (case in c("default", "local.success", "local.failure", "option.only")) {
    counts <- new.env(parent = emptyenv())
    counts$mode <- case %in% c("local.success", "local.failure")
    previous <- counts$mode
    counts$set <- logical()
    counts$reductions <- 0L
    original <- b2_condition()
    old <- options(npRmpi.local.regression.mode = case != "default")
    env <- new.env(parent = environment(.npscoefbw_collective_transaction))
    env$.Call <- function(name, active, PACKAGE) {
      expect_identical(name, "C_np_set_local_regression_mode")
      before <- counts$mode
      counts$set <- c(counts$set, active)
      counts$mode <- active
      before
    }
    env$mpi.allreduce <- function(x, type, op, comm) {
      expect_false(counts$mode)
      counts$reductions <- counts$reductions + 1L
      if (case == "local.failure") {
        if (op == "sum") c(1L, 0L, 0L, 1L) else 1L
      } else if (length(x) == 4L) c(0L, 0L, 0L, 2L) else c(0L, 0L)
    }
    env$mpi.bcast.Robj <- function(...) {
      expect_false(counts$mode)
      original
    }
    env$.npRmpi_raise_completed_failure <- function(result) stop(result$condition)
    transaction <- .npscoefbw_collective_transaction
    environment(transaction) <- env
    state <- list2env(list(comm = 1L, rank = 0L, size = 2L, guard = NULL,
                          terminal = NULL), parent = emptyenv())
    got <- tryCatch(transaction(state, function() {
      on.exit(counts$mode <- FALSE)
      entered <- tryCatch(state$guard$enter(native = TRUE), error = identity)
      expect_identical(counts$mode, previous)
      expect_identical(getOption("npRmpi.local.regression.mode"), case != "default")
      if (inherits(entered, "condition")) stop(entered)
      state$guard$phase <- "returned"
      7
    }), error = identity)
    expect_identical(counts$set, if (case == "default") logical() else c(FALSE, previous))
    expect_false(counts$mode)
    expect_identical(got, if (case == "local.failure") original else 7)
    options(old)
  }
})

test_that("B2 hooks stay private and materialize native arguments before entry", {
  expect_null(formals(.np_nomad_search)$.native.callback.transaction)
  kernel <- gsub("[[:space:]]+", " ", paste(deparse(body(npksum.default)), collapse = " "))
  expect_match(kernel, "args <- list(...)", fixed = TRUE)
  expect_match(kernel, "internal.entry.guard$enter(native = TRUE)", fixed = TRUE)
  expect_lt(regexpr("args <- list(...)", kernel, fixed = TRUE)[1L],
            regexpr("internal.entry.guard$enter(native = TRUE)", kernel, fixed = TRUE)[1L])
  expect_match(kernel, ".npRmpi_with_local_regression", fixed = TRUE)
  expect_match(kernel, "retain <- function(expr) tryCatch(force(expr)", fixed = TRUE)
  expect_match(kernel, "as.call(list(retain, reentry.call))", fixed = TRUE)
  expect_match(kernel, "eval.native <- if (is.null(internal.entry.guard)) .Call", fixed = TRUE)
  generic <- gsub("[[:space:]]+", " ", paste(deparse(body(.np_nomad_search)), collapse = " "))
  expect_match(generic, "value <- as.numeric(wrapped_eval(point)[1L])", fixed = TRUE)
  expect_match(generic, "function() native.eval.original(point), point", fixed = TRUE)
})
