test_that("regression percentile tasks replay the coordinator generator", {
  skip_if(getRversion() < "3.6.0")
  kind <- RNGkind(); had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({do.call(RNGkind, as.list(kind)); if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)}, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  state <- .np_reg_percentile_rng_state(91)
  contract <- list(runtime = paste(R.version$major, R.version$minor, sep = "."),
    header = state[1L], length = length(state), kind = RNGkind())
  expected <- .np_reg_percentile_with_seed(91, runif(37))
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(17)
  worker <- .Random.seed
  expect_identical(.np_reg_percentile_with_state(state, contract, runif(37)), expected)
  expect_identical(.Random.seed, worker)
  expect_error(.np_reg_percentile_with_state(NULL, contract, stop("entered")), "invalid regression percentile")
  expect_error(.np_reg_percentile_with_state(state, contract, stop("primary")), "primary")
  expect_identical(.Random.seed, worker)
  state[20L] <- NA_integer_
  expect_silent(.np_reg_percentile_rng_validate(state, contract))
})


test_that("task RNG cleanup survives throwing warning handlers", {
  skip_if(getRversion() < "3.6.0")
  kind <- RNGkind(); had <- exists('.Random.seed', .GlobalEnv, inherits = FALSE)
  old <- if (had) get('.Random.seed', .GlobalEnv) else NULL
  on.exit({suppressWarnings(do.call(RNGkind, as.list(kind)))
    if (had) assign('.Random.seed', old, .GlobalEnv) else
      if (exists('.Random.seed', .GlobalEnv, inherits = FALSE)) rm('.Random.seed', envir = .GlobalEnv)})
  RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
  state <- .np_reg_percentile_rng_state(91)
  contract <- list(runtime = paste(R.version$major, R.version$minor, sep = '.'),
    header = state[1L], length = length(state), kind = RNGkind())
  bad <- state; bad[1L] <- 99999L; invalid <- contract; invalid$header <- bad[1L]
  expect_error(.np_reg_percentile_rng_validate(bad, invalid), 'invalid regression percentile')
  invalid <- contract; invalid$kind[1L] <- 'Super-Duper'
  expect_error(.np_reg_percentile_rng_validate(state, invalid), 'invalid regression percentile')
  for (absent in c(FALSE, TRUE)) {
    suppressWarnings(RNGkind("L'Ecuyer-CMRG", 'Inversion', 'Rounding'))
    set.seed(42)
    if (absent) rm('.Random.seed', envir = .GlobalEnv)
    entry <- if (!absent) .Random.seed else NULL
    expected.kind <- c("L'Ecuyer-CMRG", 'Inversion', 'Rounding')
    err <- tryCatch(withCallingHandlers(
      .np_reg_percentile_with_state(state, contract, {runif(3); warning('body')}),
      warning = function(w) stop('handler primary')), error = identity)
    expect_match(conditionMessage(err), if (absent) 'caller seed is absent' else 'handler primary')
    expect_identical(if (exists('.Random.seed', .GlobalEnv, inherits = FALSE)) .Random.seed else NULL, entry)
    if (!absent) expect_identical(RNGkind(), expected.kind)
  }
})


test_that("task RNG validation respects the runtime kind layout", {
  # R before 3.6 has two kind fields and the Rounding sampler is implicit.
  validate <- .np_reg_percentile_rng_validate
  environment(validate) <- list2env(list(getRversion = function() numeric_version('3.5.0')),
    parent = environment(validate))
  state <- c(403L, rep.int(1L, 625L))
  contract <- list(runtime = paste(R.version$major, R.version$minor, sep = '.'),
    header = state[1L], length = length(state), kind = c('Mersenne-Twister', 'Inversion'))
  expect_silent(validate(state, contract))
  contract$kind <- c(contract$kind, 'Rounding')
  expect_error(validate(state, contract), 'invalid regression percentile')
})
