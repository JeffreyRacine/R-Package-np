test_that("single-index service owns collective RNG before coordinator phases", {
  skip_on_cran()
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "installed npRmpi unavailable for service ownership test")
  res <- npRmpi_run_rscript_subprocess(lines = c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(nslaves = 2L, quiet = TRUE)",
    "mpi.bcast.cmd({",
    "  options(np.messages = FALSE, np.tree = FALSE)",
    "  assign('.rng.phases', character(), .GlobalEnv); assign('.service.evals', 0L, .GlobalEnv)",
    "  trace('.npRmpi_rng_scope_agree', where = asNamespace('npRmpi'),",
    "    print = FALSE, tracer = quote(assign('.rng.phases', c(get('.rng.phases', .GlobalEnv), phase), .GlobalEnv)))",
    "  trace('.npindexbw_eval_objective_service_traced', where = asNamespace('npRmpi'),",
    "    print = FALSE, tracer = quote(assign('.service.evals', get('.service.evals', .GlobalEnv) + 1L, .GlobalEnv)))",
    "}, caller.execute = TRUE)",
    "set.seed(42); x1 <- runif(30, -1, 1); x2 <- runif(30, -1, 1)",
    "y <- (x1 - x2)^3 + rnorm(30, sd = .2); saved <- .Random.seed",
    "for (type in c('fixed', 'generalized_nn', 'adaptive_nn')) {",
    "  fit <- npindex(y ~ x1 + x2, nomad = TRUE, bwtype = type,",
    "    nmulti = 2L, random.seed = 73L, degree.max = 1L,",
    "    degree.verify = TRUE, nomad.opts = list(MAX_BB_EVAL = 8L))",
    "  stopifnot(inherits(fit, 'singleindex'), identical(.Random.seed, saved))",
    "}",
    "witness <- get('.npRmpi_bcast_cmd_expr', asNamespace('npRmpi'))(",
    "  quote(get('mpi.gather.Robj', asNamespace('npRmpi'))(",
    "    list(phases = get('.rng.phases', .GlobalEnv), evals = get('.service.evals', .GlobalEnv)), simplify = FALSE)),",
    "  caller.execute = TRUE)",
    "stopifnot(length(witness) == 3L, all(vapply(witness, function(x)",
    "  x$evals > 0L && all(c('admission','activation') %in% x$phases), TRUE)))",
    "cat('SERVICE_RNG_OWNERSHIP_OK\\n')"
  ), timeout = 100L, env = env)
  info <- paste(res$output, collapse = "\n")
  expect_identical(res$status, 0L, info = info)
  expect_true(any(grepl("SERVICE_RNG_OWNERSHIP_OK", res$output, fixed = TRUE)), info = info)
})
