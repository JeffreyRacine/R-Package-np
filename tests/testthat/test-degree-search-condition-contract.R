degree_search_contract_call <- function(eval_fun) {
  getFromNamespace(".np_degree_search", "npRmpi")(
    method = "exhaustive",
    candidates = list(0:1),
    baseline_degree = 0L,
    start_degree = 0L,
    restarts = 0L,
    max_cycles = 2L,
    verify = FALSE,
    eval_fun = eval_fun,
    direction = "min",
    trace_level = "none"
  )
}

test_that("cell degree search propagates unexpected conditions without a survivor", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  search <- getFromNamespace(".np_degree_search", "npRmpi")
  for (method in c("coordinate", "exhaustive")) for (bad in 0:2) {
    visits <- integer()
    unexpected <- errorCondition("diagnostic unexpected evaluator failure",
      class = "np_test_configuration_error", call = quote(candidate_site()),
      evidence = list(token = 934L, degree = bad))
    condition <- tryCatch(search(method = method, candidates = list(0:2),
      baseline_degree = 0L, start_degree = 0L, verify = TRUE,
      eval_fun = function(degree) {
        visits <<- c(visits, degree)
        if (degree == bad) stop(unexpected)
        list(objective = 2 - degree, payload = list(degree = degree))
      }), error = identity)
    expect_identical(condition, unexpected)
    expect_identical(tail(visits, 1L), as.integer(bad))
    expect_false(is.list(condition) && isTRUE(condition$certified))
  }
  # Identical text from distinct sites does not justify swallowing either error.
  first <- errorCondition("same text", call = quote(first_site()), token = 7L)
  expect_identical(tryCatch(degree_search_contract_call(function(degree)
    stop(if (degree == 0L) first else errorCondition("same text",
      call = quote(second_site())))), error = identity), first)
})

test_that("only typed mathematical admission failures are pruned", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  search <- getFromNamespace(".np_degree_search", "npRmpi")
  for (method in c("coordinate", "exhaustive")) for (kind in c("degree", "nn")) {
    rejected <- if (kind == "degree")
      errorCondition("inadmissible degree", class = "np_degree_candidate_invalid",
                     call = NULL, degree = 1L)
    else npRmpi:::.np_nn_candidate_invalid_condition("invalid raw NN objective",
                     owner = "contract", point = 1, raw.objective = Inf)
    result <- search(method = method, candidates = list(0:2),
      baseline_degree = 0L, start_degree = 0L, verify = TRUE,
      eval_fun = function(degree) {
        if (degree == 1L) stop(rejected)
        list(objective = 2 - degree, payload = list(degree = degree))
      })
    expect_identical(result$best_payload$degree, 2L)
    expect_true(result$completed && result$certified)
    expect_equal(sum(result$trace$status == "error"), 1L)
    expect_identical(tryCatch(degree_search_contract_call(function(degree)
      stop(rejected)), error = identity), rejected)
  }
  expect_error(degree_search_contract_call(function(degree)
    list(objective = Inf, payload = NULL)),
    "automatic degree search failed to obtain any admissible fitted model", fixed = TRUE)
})

test_that("LP admission errors are typed but malformed inputs are not", {
  admit <- getFromNamespace("npValidateLpBasisAdmission", "npRmpi")
  degrees <- list(glp = c(2L, 2L), additive = c(2L, 3L), tensor = c(1L, 2L))
  for (basis in names(degrees)) for (nobs in 5:7) {
    out <- tryCatch(admit(basis, degrees[[basis]], nobs, where = "contract"), error = identity)
    if (nobs == 7L) expect_identical(out$status, "valid") else {
      expect_s3_class(out, "np_degree_candidate_invalid")
      expect_identical(out$degree, degrees[[basis]])
      expect_identical(out$admission$status, "rank_inadmissible")
      expect_null(conditionCall(out))
    }
  }
  capacity <- tryCatch(admit("tensor", rep(10L, 5L), 200000L), error = identity)
  expect_s3_class(capacity, "np_degree_candidate_invalid")
  expect_identical(capacity$admission$status, "capacity_exceeded")
  for (args in list(list("invalid", 1L, 10L), list("glp", -1L, 10L),
                    list("glp", .5, 10L), list("glp", 1L, NA_real_))) {
    condition <- tryCatch(do.call(admit, args), error = identity)
    expect_s3_class(condition, "error")
    expect_false(inherits(condition, "np_degree_candidate_invalid"))
  }
})

test_that("exceptional search unwind finishes progress and permits reuse", {
  old <- options(np.messages = TRUE)
  on.exit(options(old), add = TRUE)
  begin <- abort <- end <- 0L
  local_mocked_bindings(
    .np_degree_progress_begin = function(...) { begin <<- begin + 1L; list(id = 1L) },
    .np_degree_progress_step = function(state, ...) state,
    .np_progress_abort = function(...) { abort <<- abort + 1L; invisible(NULL) },
    .np_degree_progress_end = function(state, ...) { end <<- end + 1L; state },
    .package = "npRmpi")
  expected <- errorCondition("failed after progress began", class = "np_test_error")
  condition <- tryCatch(degree_search_contract_call(function(degree) {
    if (degree == 1L) stop(expected)
    list(objective = 2, payload = list(degree = degree))
  }), error = identity)
  expect_identical(condition, expected)
  expect_identical(c(begin, abort, end), c(1L, 1L, 0L))
  expect_identical(degree_search_contract_call(function(degree)
    list(objective = 2 - degree, payload = list(degree = degree)))$best_payload$degree, 1L)
  expect_identical(c(begin, abort, end), c(2L, 1L, 1L))
})
test_that("real admission pruning and intentional interrupts retain their contracts", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  search <- getFromNamespace(".np_degree_search", "npRmpi")
  admitted <- search(method = "exhaustive", candidates = list(0:3),
    baseline_degree = 0L, start_degree = 0L, eval_fun = function(degree) {
      npRmpi:::npValidateLpBasisAdmission("glp", degree, 3L)
      list(objective = -degree, payload = list(degree = degree))
    })
  expect_identical(admitted$best_payload$degree, 1L)
  expect_equal(sum(admitted$trace$status == "error"), 2L)
  interrupted <- degree_search_contract_call(function(degree) {
    if (degree == 1L)
      stop(structure(list(message = "interrupt"), class = c("interrupt", "condition")))
    list(objective = 1, payload = list(degree = degree))
  })
  expect_true(interrupted$interrupted)
  expect_false(interrupted$completed || interrupted$certified)
  expect_identical(interrupted$best_payload$degree, 0L)
})

test_that("cell degree search success and accounting remain unchanged", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  result <- degree_search_contract_call(function(degree) {
    objective <- if (degree[[1L]] == 0L) 2 else 1
    list(objective = objective, payload = list(degree = degree))
  })

  expect_identical(result$best_payload$degree, 1L)
  expect_identical(result$n.unique, 2L)
  expect_identical(result$n.cached, result$n.visits - result$n.unique)
})

test_that("automatic degree search exposes bounded-kernel contracts", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI slave pool is unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  i <- seq_len(18L)
  dat <- data.frame(x = ((i * sqrt(2)) %% 1))
  dat$y <- 0.2 + 0.6 * dat$x + sin(i) / 50

  expect_error(
    npreg(y ~ x, data = dat, nomad = "auto", ckertype = "beta"),
    'beta kernels require ckerbound = "fixed" or "range" with finite lower and upper bounds',
    fixed = TRUE
  )
  expect_error(
    npreg(
      y ~ x, data = dat, nomad = "auto", ckertype = "beta",
      ckerbound = "fixed"
    ),
    "'ckerbound' requires 'ckerlb' and 'ckerub'.",
    fixed = TRUE
  )
  for (kernel in c("gaussian", "epanechnikov", "uniform")) {
    expect_error(
      npreg(
        y ~ x, data = dat, nomad = "auto", ckertype = kernel,
        ckerbound = "fixed"
      ),
      "'ckerbound' requires 'ckerlb' and 'ckerub'.",
      fixed = TRUE
    )
  }
  expect_error(
    npcdens(y ~ x, data = dat, nomad = "auto", cxkertype = "beta"),
    paste0(
      "conditional density explanatory beta kernel requires ",
      'cxkerbound = "fixed" or "range" with finite cxkerlb and cxkerub'
    ),
    fixed = TRUE
  )
  expect_error(
    do.call(npcdist, list(
      formula = y ~ x, data = dat, nomad = "auto", cykertype = "beta"
    )),
    paste0(
      "conditional distribution dependent beta kernel requires ",
      'cykerbound = "fixed" or "range" with finite cykerlb and cykerub'
    ),
    fixed = TRUE
  )
  expect_error(
    npcdens(
      y ~ x, data = dat, nomad = "auto", cxkertype = "gaussian",
      cxkerbound = "fixed"
    ),
    "'cxkerbound' requires 'cxkerlb' and 'cxkerub'.",
    fixed = TRUE
  )
  expect_error(
    do.call(npcdist, list(
      formula = y ~ x, data = dat, nomad = "auto",
      cykertype = "gaussian", cykerbound = "fixed"
    )),
    "'cykerbound' requires 'cykerlb' and 'cykerub'.",
    fixed = TRUE
  )
})
