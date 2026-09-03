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

test_that("cell degree search preserves only uniform evaluator conditions", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  shared <- errorCondition(
    "shared configuration failure",
    class = "np_test_configuration_error",
    call = quote(validate_configuration())
  )
  uniform <- tryCatch(
    degree_search_contract_call(function(degree) stop(shared)),
    error = identity
  )
  expect_s3_class(uniform, "np_test_configuration_error")
  expect_identical(conditionMessage(uniform), "shared configuration failure")
  expect_identical(conditionCall(uniform), quote(validate_configuration()))

  heterogeneous <- tryCatch(
    degree_search_contract_call(function(degree) {
      stop(if (degree[[1L]] == 0L) "first failure" else "second failure")
    }),
    error = identity
  )
  expect_s3_class(heterogeneous, "simpleError")
  expect_identical(
    conditionMessage(heterogeneous),
    "automatic degree search failed to obtain any admissible fitted model"
  )

  different_calls <- tryCatch(
    degree_search_contract_call(function(degree) {
      call <- if (degree[[1L]] == 0L) quote(first_site()) else quote(second_site())
      stop(errorCondition("shared text", call = call))
    }),
    error = identity
  )
  expect_s3_class(different_calls, "simpleError")
  expect_identical(
    conditionMessage(different_calls),
    "automatic degree search failed to obtain any admissible fitted model"
  )

  mixed <- tryCatch(
    degree_search_contract_call(function(degree) {
      if (degree[[1L]] == 0L)
        stop("candidate-local failure")
      list(objective = Inf, payload = NULL)
    }),
    error = identity
  )
  expect_s3_class(mixed, "simpleError")
  expect_identical(
    conditionMessage(mixed),
    "automatic degree search failed to obtain any admissible fitted model"
  )
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
