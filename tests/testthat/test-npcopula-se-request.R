test_that("copula extraction reads stored errors without refitting", {
  stored <- structure(list(copula = c(.2, .8),
                           copulaerr = c(NA_real_, 0)), class = "npcopula")
  expect_identical(se(stored), stored$copulaerr)
  stored$se <- TRUE
  expect_identical(se(stored), stored$copulaerr)
  omitted <- stored
  omitted$se <- FALSE
  expect_error(se(omitted), "npcopula(bws = omitted$bws, se = TRUE",
               fixed = TRUE)
  omitted$se <- NULL
  omitted$copulaerr <- NULL
  expect_error(se(omitted), "without repeating bandwidth search", fixed = TRUE)
  expect_error(predict(omitted, se.fit = TRUE),
               "were not computed", fixed = TRUE)
  for (output in c("object", "data"))
    expect_error(predict(omitted, se.fit = TRUE, output = output),
                 "were not computed", fixed = TRUE)
  expect_identical(predict(omitted), omitted$copula)
})

test_that("copula SE fitting request is independent of point evaluation", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  dat <- data.frame(x = seq(-1, 1, length.out = 24),
                    y = sin(seq_len(24)) + seq_len(24)/12)
  ns <- asNamespace(getNamespaceName(environment(npcopula)))
  old.target <- get(".npcopula_asymptotic_se", ns)
  eval.grid <- get(".npcopula_eval_xgrid", ns)
  for (target in c("density", "distribution")) {
    make.bw <- if (target == "density") npudensbw else npudistbw
    bw <- make.bw(dat = dat, bws = c(.4, .6), bandwidth.compute = FALSE)
    for (u in list(NULL, data.frame(x = c(.3, .7), y = c(.35, .65)))) {
      off <- npcopula(bws = bw, data = dat, u = u, n.quasi.inv = 24)
      full <- npcopula(bws = bw, data = dat, u = u, n.quasi.inv = 24,
                       se = TRUE)
      expect_identical(off$se, FALSE)
      expect_null(off$copulaerr)
      expect_identical(full$se, TRUE)
      expect_identical(as.data.frame(off), as.data.frame(full))
      expect_identical(se(full), old.target(full, dat, eval.grid(full)))
      expect_error(se(off), "without repeating bandwidth search", fixed = TRUE)
      expect_identical(predict(full, se.fit = TRUE)$se.fit, se(full))
    }
    grid <- data.frame(x = c(.3, .7), y = c(.35, .65))
    pred <- predict(off, u = grid, se.fit = TRUE, n.quasi.inv = 24)
    expected <- npcopula(bws = bw, data = dat, u = grid,
                         n.quasi.inv = 24, se = TRUE)
    expect_identical(pred$fit, fitted(expected))
    expect_identical(pred$se.fit, se(expected))
  }
  for (method in c("formula", "default")) {
    expect_identical(formals(getS3method("npcopula", method))$se, FALSE)
    public.formals <- names(formals(getS3method("npcopula", method)))
    expect_gt(match("se", public.formals), match("...", public.formals))
  }
  expect_error(npcopula(se = NA), "'se' must", fixed = TRUE)
  expect_error(npcopula(s.e = stop("must not evaluate")),
               "did you mean 'se'", fixed = TRUE)
  formula.fit <- npcopula(~ x + y, data = dat, evaluation = "sample",
                         bandwidth.compute = FALSE, se = TRUE)
  formula.reference <- npcopula(bws = formula.fit$bws, data = dat, se = TRUE)
  expect_identical(fitted(formula.fit), fitted(formula.reference))
  expect_identical(se(formula.fit), se(formula.reference))
  expect_error(predict(formula.fit, se.fit = TRUE, se = FALSE),
               "conflicting 'se' and 'se.fit'", fixed = TRUE)

  # Test-only closure: uncertainty recovery must not execute on the OFF path.
  guarded <- get(".npcopula_eval", ns)
  environment(guarded) <- new.env(parent = environment(guarded))
  environment(guarded)$.npcopula_asymptotic_se <- function(...)
    stop("uncertainty owner entered", call. = FALSE)
  args <- list(bws = bw, data = dat, target = "distribution",
               evaluation = "sample", neval = 2,
               n.quasi.inv = 24, er.quasi.inv = 1)
  expect_identical(do.call(guarded, c(args, list(se = FALSE)))$se, FALSE)
  expect_error(do.call(guarded, c(args, list(se = TRUE))),
               "uncertainty owner entered", fixed = TRUE)

  # The new fitting-time demand announces its real stage before entering it.
  events <- list()
  observed <- get(".npcopula_eval", ns)
  environment(observed) <- new.env(parent = environment(observed))
  record <- function(label) {
    events[[length(events) + 1L]] <<- list(label = label, at = proc.time()[3L])
  }
  environment(observed)$.npcopula_progress_step <- function(state, done, detail) {
    record(detail)
    get(".npcopula_progress_step", ns)(state, done, detail)
  }
  environment(observed)$.npcopula_asymptotic_se <- function(...) {
    record("uncertainty-entry")
    result <- old.target(...)
    record("uncertainty-complete")
    result
  }
  result <- do.call(observed, c(args, list(se = TRUE)))
  labels <- vapply(events, function(event) event$label, "")
  expect_lt(match("asymptotic standard errors", labels),
            match("uncertainty-entry", labels))
  expect_lt(match("uncertainty-entry", labels),
            match("uncertainty-complete", labels))
  expect_identical(se(result), old.target(result, dat, eval.grid(result)))
})
